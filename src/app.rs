//! Main application state and UI

use eframe::egui;
use std::sync::mpsc::{channel, Receiver};
use std::thread;

use crate::analysis::{
    align_all_sequences, parse_sequence_set, parse_template_fasta, reverse_complement,
    run_pairwise_screening, AlignmentResults, AnalysisMethod, AnalysisParams, ProgressUpdate,
    ScreeningResults, SequenceSet, TemplateSequence, ThreadCount,
};

/// Application state
pub struct OligoscreenApp {
    // Input state — two separate inputs
    template_input: String,
    sequence_set_input: String,
    template: Option<TemplateSequence>,
    sequence_set: Option<SequenceSet>,
    template_error: Option<String>,
    sequence_set_error: Option<String>,

    // Alignment state
    alignment_results: Option<AlignmentResults>,
    is_aligning: bool,
    alignment_progress: Option<ProgressUpdate>,
    alignment_progress_rx: Option<Receiver<ProgressUpdate>>,
    alignment_results_rx: Option<Receiver<AlignmentResults>>,

    // Analysis parameters
    params: AnalysisParams,
    method_selection: MethodSelection,
    thread_selection: ThreadSelection,
    manual_thread_count: usize,
    incremental_limit_ambiguities: bool,
    incremental_max_ambiguities: u32,

    // Analysis state (screening)
    is_analyzing: bool,
    analysis_progress: Option<ProgressUpdate>,
    progress_rx: Option<Receiver<ProgressUpdate>>,
    results_rx: Option<Receiver<ScreeningResults>>,

    // Results state
    results: Option<ScreeningResults>,
    selected_length: Option<u32>,
    selected_position: Option<usize>,
    show_detail_window: bool,

    // Detail window display options
    detail_show_reverse_complement: bool,
    detail_show_codon_spacing: bool,

    // View state
    current_tab: Tab,

    // Save/Load
    save_error: Option<String>,
    load_error: Option<String>,
    pending_save: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Tab {
    Input,
    Analysis,
    Results,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum MethodSelection {
    NoAmbiguities,
    FixedAmbiguities,
    Incremental,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ThreadSelection {
    Auto,
    Manual,
}

impl Default for OligoscreenApp {
    fn default() -> Self {
        let available_threads = std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1);
        Self {
            template_input: String::new(),
            sequence_set_input: String::new(),
            template: None,
            sequence_set: None,
            template_error: None,
            sequence_set_error: None,
            alignment_results: None,
            is_aligning: false,
            alignment_progress: None,
            alignment_progress_rx: None,
            alignment_results_rx: None,
            params: AnalysisParams::default(),
            method_selection: MethodSelection::NoAmbiguities,
            thread_selection: ThreadSelection::Auto,
            manual_thread_count: available_threads,
            incremental_limit_ambiguities: false,
            incremental_max_ambiguities: 3,
            is_analyzing: false,
            analysis_progress: None,
            progress_rx: None,
            results_rx: None,
            results: None,
            selected_length: None,
            selected_position: None,
            show_detail_window: false,
            detail_show_reverse_complement: false,
            detail_show_codon_spacing: true,
            current_tab: Tab::Input,
            save_error: None,
            load_error: None,
            pending_save: false,
        }
    }
}

impl OligoscreenApp {
    pub fn new(_cc: &eframe::CreationContext<'_>) -> Self {
        Self::default()
    }

    fn parse_template(&mut self) {
        self.template_error = None;
        self.template = None;

        if self.template_input.trim().is_empty() {
            return;
        }

        match parse_template_fasta(&self.template_input) {
            Ok(t) => self.template = Some(t),
            Err(e) => self.template_error = Some(e),
        }
    }

    fn parse_sequence_set(&mut self) {
        self.sequence_set_error = None;
        self.sequence_set = None;

        if self.sequence_set_input.trim().is_empty() {
            return;
        }

        match parse_sequence_set(&self.sequence_set_input) {
            Ok(s) => self.sequence_set = Some(s),
            Err(e) => self.sequence_set_error = Some(e),
        }
    }

    fn start_alignment(&mut self) {
        let Some(template) = &self.template else {
            return;
        };
        let Some(seq_set) = &self.sequence_set else {
            return;
        };

        let template_clone = template.clone();
        let seq_set_clone = seq_set.clone();
        let scoring = self.params.alignment_scoring.clone();

        let (progress_tx, progress_rx) = channel();
        let (results_tx, results_rx) = channel();

        self.alignment_progress_rx = Some(progress_rx);
        self.alignment_results_rx = Some(results_rx);
        self.is_aligning = true;
        self.alignment_progress = None;
        self.alignment_results = None;
        self.results = None;

        thread::spawn(move || {
            let results =
                align_all_sequences(&template_clone, &seq_set_clone, &scoring, Some(&progress_tx));
            let _ = results_tx.send(results);
        });
    }

    fn start_analysis(&mut self) {
        let Some(alignment_results) = &self.alignment_results else {
            return;
        };

        // Update method from selection
        self.params.method = match self.method_selection {
            MethodSelection::NoAmbiguities => AnalysisMethod::NoAmbiguities,
            MethodSelection::FixedAmbiguities => {
                AnalysisMethod::FixedAmbiguities(self.params.method.get_fixed_ambiguities())
            }
            MethodSelection::Incremental => {
                let max_amb = if self.incremental_limit_ambiguities {
                    Some(self.incremental_max_ambiguities)
                } else {
                    None
                };
                AnalysisMethod::Incremental(self.params.method.get_incremental_pct(), max_amb)
            }
        };

        // Update thread count from selection
        self.params.thread_count = match self.thread_selection {
            ThreadSelection::Auto => ThreadCount::Auto,
            ThreadSelection::Manual => ThreadCount::Fixed(self.manual_thread_count),
        };

        let alignment_clone = alignment_results.clone();
        let params_clone = self.params.clone();

        let (progress_tx, progress_rx) = channel();
        let (results_tx, results_rx) = channel();

        self.progress_rx = Some(progress_rx);
        self.results_rx = Some(results_rx);
        self.is_analyzing = true;
        self.analysis_progress = None;

        thread::spawn(move || {
            let results =
                run_pairwise_screening(&alignment_clone, &params_clone, Some(progress_tx));
            let _ = results_tx.send(results);
        });
    }

    fn check_progress(&mut self) {
        // Check alignment progress
        if let Some(rx) = &self.alignment_progress_rx {
            while let Ok(progress) = rx.try_recv() {
                self.alignment_progress = Some(progress);
            }
        }

        if let Some(rx) = &self.alignment_results_rx {
            if let Ok(results) = rx.try_recv() {
                self.alignment_results = Some(results);
                self.is_aligning = false;
                self.alignment_progress_rx = None;
                self.alignment_results_rx = None;
            }
        }

        // Check analysis progress
        if let Some(rx) = &self.progress_rx {
            while let Ok(progress) = rx.try_recv() {
                self.analysis_progress = Some(progress);
            }
        }

        if let Some(rx) = &self.results_rx {
            if let Ok(results) = rx.try_recv() {
                self.selected_length = results.results_by_length.keys().min().copied();
                self.results = Some(results);
                self.is_analyzing = false;
                self.progress_rx = None;
                self.results_rx = None;
                self.current_tab = Tab::Results;
            }
        }
    }

    fn save_results(&mut self) {
        let Some(results) = &self.results else {
            self.save_error = Some("No results to save".to_string());
            return;
        };

        if let Some(path) = rfd::FileDialog::new()
            .add_filter("JSON", &["json"])
            .set_file_name("screening_results.json")
            .save_file()
        {
            match serde_json::to_string_pretty(results) {
                Ok(json) => {
                    if let Err(e) = std::fs::write(&path, json) {
                        self.save_error = Some(format!("Failed to write file: {}", e));
                    } else {
                        self.save_error = None;
                    }
                }
                Err(e) => {
                    self.save_error = Some(format!("Failed to serialize: {}", e));
                }
            }
        }
    }

    fn load_results(&mut self) {
        if let Some(path) = rfd::FileDialog::new()
            .add_filter("JSON", &["json"])
            .pick_file()
        {
            match std::fs::read_to_string(&path) {
                Ok(json) => match serde_json::from_str::<ScreeningResults>(&json) {
                    Ok(results) => {
                        self.selected_length = results.results_by_length.keys().min().copied();
                        self.results = Some(results);
                        self.load_error = None;
                        self.current_tab = Tab::Results;
                    }
                    Err(e) => {
                        self.load_error = Some(format!("Failed to parse: {}", e));
                    }
                },
                Err(e) => {
                    self.load_error = Some(format!("Failed to read file: {}", e));
                }
            }
        }
    }

    fn load_template_file(&mut self) {
        if let Some(path) = rfd::FileDialog::new()
            .add_filter("FASTA", &["fasta", "fa", "fna", "fas", "txt"])
            .pick_file()
        {
            match std::fs::read_to_string(&path) {
                Ok(content) => {
                    self.template_input = content;
                    self.parse_template();
                }
                Err(e) => {
                    self.template_error = Some(format!("Failed to read file: {}", e));
                }
            }
        }
    }

    fn load_sequence_set_file(&mut self) {
        if let Some(path) = rfd::FileDialog::new()
            .add_filter("FASTA", &["fasta", "fa", "fna", "fas", "txt"])
            .pick_file()
        {
            match std::fs::read_to_string(&path) {
                Ok(content) => {
                    self.sequence_set_input = content;
                    self.parse_sequence_set();
                }
                Err(e) => {
                    self.sequence_set_error = Some(format!("Failed to read file: {}", e));
                }
            }
        }
    }
}

impl AnalysisMethod {
    fn get_fixed_ambiguities(&self) -> u32 {
        match self {
            AnalysisMethod::FixedAmbiguities(n) => *n,
            _ => 1,
        }
    }

    fn get_incremental_pct(&self) -> u32 {
        match self {
            AnalysisMethod::Incremental(pct, _) => *pct,
            _ => 50,
        }
    }

    fn get_incremental_max_amb(&self) -> Option<u32> {
        match self {
            AnalysisMethod::Incremental(_, max_amb) => *max_amb,
            _ => None,
        }
    }
}

impl eframe::App for OligoscreenApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        if self.is_aligning || self.is_analyzing {
            self.check_progress();
            ctx.request_repaint();
        }

        // Handle deferred save
        if self.pending_save {
            self.pending_save = false;
            self.save_results();
        }

        // Top menu bar
        egui::TopBottomPanel::top("menu_bar").show(ctx, |ui| {
            egui::menu::bar(ui, |ui| {
                ui.menu_button("File", |ui| {
                    if ui.button("Load Template...").clicked() {
                        self.load_template_file();
                        ui.close_menu();
                    }
                    if ui.button("Load Sequences...").clicked() {
                        self.load_sequence_set_file();
                        ui.close_menu();
                    }
                    ui.separator();
                    if ui.button("Load Results...").clicked() {
                        self.load_results();
                        ui.close_menu();
                    }
                    if ui.button("Save Results...").clicked() {
                        self.save_results();
                        ui.close_menu();
                    }
                });
            });
        });

        // Tab bar
        egui::TopBottomPanel::top("tabs").show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.selectable_value(&mut self.current_tab, Tab::Input, "Input Data");
                ui.selectable_value(&mut self.current_tab, Tab::Analysis, "Analysis Setup");
                ui.selectable_value(&mut self.current_tab, Tab::Results, "Results");
            });
        });

        // Status bar
        egui::TopBottomPanel::bottom("status").show(ctx, |ui| {
            ui.horizontal(|ui| {
                if self.is_aligning {
                    ui.spinner();
                    if let Some(ref progress) = self.alignment_progress {
                        ui.label(&progress.message);
                    } else {
                        ui.label("Starting alignment...");
                    }
                } else if self.is_analyzing {
                    ui.spinner();
                    if let Some(ref progress) = self.analysis_progress {
                        ui.label(&progress.message);
                    } else {
                        ui.label("Starting analysis...");
                    }
                } else if let Some(ref results) = self.results {
                    ui.label(format!(
                        "Results: {} sequences, {} bp template",
                        results.total_sequences, results.alignment_length
                    ));
                } else if let Some(ref ar) = self.alignment_results {
                    ui.label(format!(
                        "Aligned: {}/{} sequences valid, {} bp template",
                        ar.valid_count,
                        ar.total_count,
                        ar.template.sequence.len()
                    ));
                } else {
                    ui.label("Load template and sequences to begin");
                }
            });
        });

        // Main content
        egui::CentralPanel::default().show(ctx, |ui| {
            match self.current_tab {
                Tab::Input => self.show_input_tab(ui),
                Tab::Analysis => self.show_analysis_tab(ui),
                Tab::Results => self.show_results_tab(ui),
            }
        });

        // Detail window
        if self.show_detail_window {
            self.show_variant_detail_window(ctx);
        }
    }
}

impl OligoscreenApp {
    fn show_input_tab(&mut self, ui: &mut egui::Ui) {
        ui.heading("Input Data");
        ui.separator();

        let available_height = ui.available_height();
        let half_height = (available_height - 120.0) / 2.0;

        // Template section
        ui.group(|ui| {
            ui.horizontal(|ui| {
                ui.strong("Template Sequence");
                ui.separator();
                if ui.button("Load File").clicked() {
                    self.load_template_file();
                }
                if ui.button("Clear").clicked() {
                    self.template_input.clear();
                    self.template = None;
                    self.template_error = None;
                    self.alignment_results = None;
                }
            });

            ui.label("Single FASTA sequence (e.g., RefSeq of the target):");

            egui::ScrollArea::vertical()
                .id_salt("template_scroll")
                .max_height(half_height.max(80.0))
                .show(ui, |ui| {
                    let response = ui.add(
                        egui::TextEdit::multiline(&mut self.template_input)
                            .font(egui::TextStyle::Monospace)
                            .desired_width(f32::INFINITY)
                            .desired_rows(6),
                    );
                    if response.changed() {
                        self.parse_template();
                        self.alignment_results = None;
                    }
                });

            if let Some(ref error) = self.template_error {
                ui.colored_label(egui::Color32::RED, format!("Error: {}", error));
            }
            if let Some(ref t) = self.template {
                ui.colored_label(
                    egui::Color32::GREEN,
                    format!("{}: {} bp", t.name, t.sequence.len()),
                );
            }
        });

        ui.add_space(5.0);

        // Sequence set section
        ui.group(|ui| {
            ui.horizontal(|ui| {
                ui.strong("Sequence Set");
                ui.separator();
                if ui.button("Load File").clicked() {
                    self.load_sequence_set_file();
                }
                if ui.button("Clear").clicked() {
                    self.sequence_set_input.clear();
                    self.sequence_set = None;
                    self.sequence_set_error = None;
                    self.alignment_results = None;
                }
            });

            ui.label("Multiple FASTA sequences (unaligned):");

            egui::ScrollArea::vertical()
                .id_salt("seqset_scroll")
                .max_height(half_height.max(80.0))
                .show(ui, |ui| {
                    let response = ui.add(
                        egui::TextEdit::multiline(&mut self.sequence_set_input)
                            .font(egui::TextStyle::Monospace)
                            .desired_width(f32::INFINITY)
                            .desired_rows(6),
                    );
                    if response.changed() {
                        self.parse_sequence_set();
                        self.alignment_results = None;
                    }
                });

            if let Some(ref error) = self.sequence_set_error {
                ui.colored_label(egui::Color32::RED, format!("Error: {}", error));
            }
            if let Some(ref s) = self.sequence_set {
                ui.colored_label(egui::Color32::GREEN, format!("{} sequences loaded", s.len()));
            }
        });

        ui.add_space(10.0);

        // Align button
        ui.horizontal(|ui| {
            let can_align = self.template.is_some()
                && self.sequence_set.is_some()
                && !self.is_aligning
                && !self.is_analyzing;

            if ui
                .add_enabled(can_align, egui::Button::new("Align Sequences"))
                .clicked()
            {
                self.start_alignment();
            }

            if self.is_aligning {
                ui.spinner();
                if let Some(ref progress) = self.alignment_progress {
                    ui.label(&progress.message);
                }
            }

            if let Some(ref ar) = self.alignment_results {
                ui.colored_label(
                    egui::Color32::GREEN,
                    format!("Aligned: {}/{} valid", ar.valid_count, ar.total_count),
                );
            }
        });
    }

    fn show_analysis_tab(&mut self, ui: &mut egui::Ui) {
        ui.heading("Analysis Setup");
        ui.separator();

        if self.alignment_results.is_none() {
            ui.colored_label(
                egui::Color32::YELLOW,
                "Please load and align sequences in the Input tab first.",
            );
            return;
        }

        egui::ScrollArea::vertical().show(ui, |ui| {
            // Analysis method selection
            ui.group(|ui| {
                ui.heading("Analysis Method");

                ui.radio_value(
                    &mut self.method_selection,
                    MethodSelection::NoAmbiguities,
                    "No Ambiguities - Find all unique exact variants",
                );

                ui.radio_value(
                    &mut self.method_selection,
                    MethodSelection::FixedAmbiguities,
                    "Fixed Ambiguities - Use up to N ambiguity codes per variant",
                );

                if self.method_selection == MethodSelection::FixedAmbiguities {
                    ui.horizontal(|ui| {
                        ui.add_space(20.0);
                        ui.label("Max ambiguities:");
                        let mut n = self.params.method.get_fixed_ambiguities();
                        if ui
                            .add(egui::DragValue::new(&mut n).range(0..=20))
                            .changed()
                        {
                            self.params.method = AnalysisMethod::FixedAmbiguities(n);
                        }
                    });
                }

                ui.radio_value(
                    &mut self.method_selection,
                    MethodSelection::Incremental,
                    "Incremental - Find variants covering X% of remaining sequences",
                );

                if self.method_selection == MethodSelection::Incremental {
                    ui.horizontal(|ui| {
                        ui.add_space(20.0);
                        ui.label("Target coverage per step (%):");
                        let mut pct = self.params.method.get_incremental_pct();
                        let max_amb = self.params.method.get_incremental_max_amb();
                        if ui
                            .add(egui::DragValue::new(&mut pct).range(1..=100))
                            .changed()
                        {
                            self.params.method = AnalysisMethod::Incremental(pct, max_amb);
                        }
                    });
                    ui.horizontal(|ui| {
                        ui.add_space(20.0);
                        ui.checkbox(
                            &mut self.incremental_limit_ambiguities,
                            "Limit ambiguities:",
                        );
                        ui.add_enabled(
                            self.incremental_limit_ambiguities,
                            egui::DragValue::new(&mut self.incremental_max_ambiguities)
                                .range(0..=20),
                        );
                        ui.label("max");
                    });
                }
            });

            ui.add_space(10.0);

            // Global options
            ui.group(|ui| {
                ui.heading("Global Options");
                ui.checkbox(
                    &mut self.params.exclude_n,
                    "Exclude N (any base) as ambiguity code",
                );
            });

            ui.add_space(10.0);

            // Oligo length range
            ui.group(|ui| {
                ui.heading("Oligo Length Range");
                ui.horizontal(|ui| {
                    ui.label("Minimum length:");
                    ui.add(
                        egui::DragValue::new(&mut self.params.min_oligo_length).range(3..=100),
                    );
                    ui.add_space(20.0);
                    ui.label("Maximum length:");
                    ui.add(
                        egui::DragValue::new(&mut self.params.max_oligo_length).range(3..=100),
                    );
                });

                if self.params.min_oligo_length > self.params.max_oligo_length {
                    self.params.max_oligo_length = self.params.min_oligo_length;
                }

                let range = self.params.max_oligo_length - self.params.min_oligo_length + 1;
                if range > 20 {
                    ui.colored_label(
                        egui::Color32::YELLOW,
                        format!(
                            "Warning: Large length range ({}) may take significant time",
                            range
                        ),
                    );
                }
            });

            ui.add_space(10.0);

            // Resolution
            ui.group(|ui| {
                ui.heading("Analysis Resolution");
                ui.horizontal(|ui| {
                    ui.label("Step size (bases):");
                    ui.add(egui::DragValue::new(&mut self.params.resolution).range(1..=100));
                });
                ui.label("Lower values = more positions analyzed, higher resolution");
            });

            ui.add_space(10.0);

            // Coverage threshold
            ui.group(|ui| {
                ui.heading("Coverage Threshold");
                ui.horizontal(|ui| {
                    ui.label("Target coverage (%):");
                    ui.add(
                        egui::DragValue::new(&mut self.params.coverage_threshold)
                            .range(1.0..=100.0),
                    );
                });
                ui.label("Number of variants needed to reach this coverage will be reported");
            });

            ui.add_space(10.0);

            // Alignment scoring
            ui.group(|ui| {
                ui.heading("Alignment Scoring");
                ui.horizontal(|ui| {
                    ui.label("Match:");
                    ui.add(
                        egui::DragValue::new(&mut self.params.alignment_scoring.match_score)
                            .range(0..=10),
                    );
                    ui.add_space(10.0);
                    ui.label("Mismatch:");
                    ui.add(
                        egui::DragValue::new(&mut self.params.alignment_scoring.mismatch_score)
                            .range(-10..=0),
                    );
                });
                ui.horizontal(|ui| {
                    ui.label("Gap open:");
                    ui.add(
                        egui::DragValue::new(&mut self.params.alignment_scoring.gap_open)
                            .range(-20..=0),
                    );
                    ui.add_space(10.0);
                    ui.label("Gap extend:");
                    ui.add(
                        egui::DragValue::new(&mut self.params.alignment_scoring.gap_extend)
                            .range(-10..=0),
                    );
                });
                ui.horizontal(|ui| {
                    ui.label("Min alignment quality:");
                    ui.add(
                        egui::DragValue::new(
                            &mut self.params.alignment_scoring.min_score_fraction,
                        )
                        .range(0.0..=1.0)
                        .speed(0.05),
                    );
                });
                ui.label("Sequences below this quality threshold are excluded");
            });

            ui.add_space(10.0);

            // Thread count
            ui.group(|ui| {
                ui.heading("Parallelization");

                let available_threads = std::thread::available_parallelism()
                    .map(|n| n.get())
                    .unwrap_or(1);

                ui.radio_value(
                    &mut self.thread_selection,
                    ThreadSelection::Auto,
                    format!("Auto ({} threads)", available_threads),
                );

                ui.horizontal(|ui| {
                    ui.radio_value(
                        &mut self.thread_selection,
                        ThreadSelection::Manual,
                        "Manual:",
                    );
                    let enabled = self.thread_selection == ThreadSelection::Manual;
                    ui.add_enabled(
                        enabled,
                        egui::DragValue::new(&mut self.manual_thread_count)
                            .range(1..=available_threads.max(32)),
                    );
                    ui.label("threads");
                });
            });

            ui.add_space(20.0);

            // Run button
            ui.horizontal(|ui| {
                let can_run = self.alignment_results.is_some()
                    && !self.is_analyzing
                    && !self.is_aligning;
                if ui
                    .add_enabled(can_run, egui::Button::new("Run Analysis"))
                    .clicked()
                {
                    self.start_analysis();
                }

                if self.is_analyzing {
                    ui.spinner();
                    if let Some(ref progress) = self.analysis_progress {
                        ui.label(&progress.message);
                    }
                }
            });
        });
    }

    fn show_results_tab(&mut self, ui: &mut egui::Ui) {
        if self.results.is_none() {
            ui.heading("Results");
            ui.separator();
            ui.label("No results yet. Run an analysis from the Analysis Setup tab.");
            return;
        }

        let has_results = self.results.is_some();

        ui.horizontal(|ui| {
            ui.heading("Results");
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                if ui.button("Save Results").clicked() && has_results {
                    self.pending_save = true;
                }
            });
        });
        ui.separator();

        // Get lengths list and template sequence
        let (lengths, template_seq) = {
            let results = self.results.as_ref().unwrap();
            let mut lengths: Vec<u32> = results.results_by_length.keys().copied().collect();
            lengths.sort();
            (lengths, results.consensus_sequence.clone())
        };

        // Length selector
        ui.horizontal(|ui| {
            ui.label("Oligo Length:");

            egui::ComboBox::from_id_salt("length_selector")
                .selected_text(
                    self.selected_length
                        .map(|l| format!("{} bp", l))
                        .unwrap_or_else(|| "Select...".to_string()),
                )
                .show_ui(ui, |ui| {
                    for length in &lengths {
                        ui.selectable_value(
                            &mut self.selected_length,
                            Some(*length),
                            format!("{} bp", length),
                        );
                    }
                });
        });

        ui.add_space(10.0);

        // Results display
        if let Some(length) = self.selected_length {
            let length_result = self
                .results
                .as_ref()
                .and_then(|r| r.results_by_length.get(&length))
                .cloned();

            if let Some(length_result) = length_result {
                let coverage_threshold = self
                    .results
                    .as_ref()
                    .map(|r| r.params.coverage_threshold)
                    .unwrap_or(95.0);
                self.show_length_results(ui, &length_result, &template_seq, coverage_threshold);
            }
        }

        // Error messages
        if let Some(ref error) = self.save_error {
            ui.colored_label(egui::Color32::RED, error);
        }
        if let Some(ref error) = self.load_error {
            ui.colored_label(egui::Color32::RED, error);
        }
    }

    fn show_length_results(
        &mut self,
        ui: &mut egui::Ui,
        length_result: &crate::analysis::LengthResult,
        template_seq: &str,
        coverage_threshold: f64,
    ) {
        let positions = &length_result.positions;
        if positions.is_empty() {
            ui.label("No positions analyzed for this length.");
            return;
        }

        let oligo_len = length_result.oligo_length as usize;

        // Summary stats
        let non_skipped: Vec<_> = positions.iter().filter(|p| !p.analysis.skipped).collect();
        let min_variants = non_skipped
            .iter()
            .map(|p| p.variants_needed)
            .min()
            .unwrap_or(0);
        let max_variants = non_skipped
            .iter()
            .map(|p| p.variants_needed)
            .max()
            .unwrap_or(0);
        let avg_variants: f64 = if non_skipped.is_empty() {
            0.0
        } else {
            non_skipped
                .iter()
                .map(|p| p.variants_needed)
                .sum::<usize>() as f64
                / non_skipped.len() as f64
        };

        ui.group(|ui| {
            ui.horizontal(|ui| {
                ui.label(format!("Positions analyzed: {}", positions.len()));
                ui.separator();
                ui.label(format!("Non-skipped: {}", non_skipped.len()));
                ui.separator();
                ui.label(format!("Min variants: {}", min_variants));
                ui.separator();
                ui.label(format!("Max variants: {}", max_variants));
                ui.separator();
                ui.label(format!("Avg variants: {:.1}", avg_variants));
            });
        });

        ui.add_space(5.0);

        // Legend
        ui.horizontal(|ui| {
            ui.label("Legend:");
            ui.add_space(5.0);
            let legend_colors = [
                (egui::Color32::from_rgb(0, 180, 0), "1 (best)"),
                (egui::Color32::from_rgb(150, 180, 0), "2-5"),
                (egui::Color32::from_rgb(255, 165, 0), "6-10"),
                (egui::Color32::from_rgb(220, 50, 50), ">10"),
                (egui::Color32::DARK_GRAY, "skipped"),
            ];
            for (color, label) in legend_colors {
                let (rect, _) =
                    ui.allocate_exact_size(egui::vec2(12.0, 12.0), egui::Sense::hover());
                ui.painter().rect_filled(rect, 2.0, color);
                ui.label(label);
                ui.add_space(5.0);
            }
        });

        ui.add_space(5.0);

        ui.label(format!(
            "Variants needed to reach {:.0}% coverage:",
            coverage_threshold
        ));

        // Vertical scrollable table
        egui::ScrollArea::vertical()
            .id_salt("results_table")
            .show(ui, |ui| {
                egui::Grid::new("position_grid")
                    .striped(true)
                    .min_col_width(40.0)
                    .spacing([8.0, 4.0])
                    .show(ui, |ui| {
                        // Header row
                        ui.strong("Pos");
                        ui.strong("Template Window");
                        ui.strong("Variants");
                        ui.strong("Coverage");
                        ui.strong("Seqs");
                        ui.strong("");
                        ui.end_row();

                        for pos_result in positions {
                            let pos = pos_result.position;
                            let end = (pos + oligo_len).min(template_seq.len());
                            let template_window = if end <= template_seq.len() && pos < end {
                                &template_seq[pos..end]
                            } else {
                                "—"
                            };

                            // Position (1-indexed)
                            ui.label(format!("{}", pos + 1));

                            // Template subsequence (monospace)
                            ui.add(
                                egui::Label::new(
                                    egui::RichText::new(template_window)
                                        .monospace()
                                        .size(11.0),
                                )
                                .wrap_mode(egui::TextWrapMode::Extend),
                            );

                            // Variants needed
                            if pos_result.analysis.skipped {
                                ui.colored_label(egui::Color32::DARK_GRAY, "skipped");
                            } else {
                                ui.label(format!("{}", pos_result.variants_needed));
                            }

                            // Coverage
                            if !pos_result.analysis.skipped {
                                ui.label(format!(
                                    "{:.1}%",
                                    pos_result.analysis.coverage_at_threshold
                                ));
                            } else {
                                ui.label("—");
                            }

                            // Sequences analyzed
                            if !pos_result.analysis.skipped {
                                ui.label(format!("{}", pos_result.analysis.sequences_analyzed));
                            } else {
                                ui.label("—");
                            }

                            // Color indicator + click
                            let color = if pos_result.analysis.skipped {
                                egui::Color32::DARK_GRAY
                            } else {
                                variant_count_color(pos_result.variants_needed)
                            };

                            let (rect, response) = ui.allocate_exact_size(
                                egui::vec2(20.0, 16.0),
                                egui::Sense::click(),
                            );
                            ui.painter().rect_filled(rect, 2.0, color);

                            if response.clicked() && !pos_result.analysis.skipped {
                                self.selected_position = Some(pos_result.position);
                                self.show_detail_window = true;
                            }

                            response.on_hover_ui(|ui| {
                                ui.label(format!("Position: {}", pos + 1));
                                if pos_result.analysis.skipped {
                                    ui.label(format!(
                                        "Skipped: {}",
                                        pos_result
                                            .analysis
                                            .skip_reason
                                            .as_deref()
                                            .unwrap_or("Unknown")
                                    ));
                                } else {
                                    ui.label(format!(
                                        "Variants needed: {}",
                                        pos_result.variants_needed
                                    ));
                                    ui.label(format!(
                                        "Coverage: {:.1}%",
                                        pos_result.analysis.coverage_at_threshold
                                    ));
                                    ui.label("Click for details");
                                }
                            });

                            ui.end_row();
                        }
                    });
            });
    }

    fn show_variant_detail_window(&mut self, ctx: &egui::Context) {
        let Some(ref results) = self.results else {
            self.show_detail_window = false;
            return;
        };

        let Some(length) = self.selected_length else {
            self.show_detail_window = false;
            return;
        };

        let Some(position) = self.selected_position else {
            self.show_detail_window = false;
            return;
        };

        let Some(length_result) = results.results_by_length.get(&length) else {
            self.show_detail_window = false;
            return;
        };

        let Some(pos_result) = length_result
            .positions
            .iter()
            .find(|p| p.position == position)
        else {
            self.show_detail_window = false;
            return;
        };

        let pos_result = pos_result.clone();
        let coverage_threshold = results.params.coverage_threshold;

        let show_reverse_complement = self.detail_show_reverse_complement;
        let show_codon_spacing = self.detail_show_codon_spacing;

        egui::Window::new(format!("Position {} Details", position + 1))
            .open(&mut self.show_detail_window)
            .default_width(600.0)
            .default_height(400.0)
            .show(ctx, |ui| {
                ui.horizontal(|ui| {
                    ui.label(format!("Position: {}", position + 1));
                    ui.separator();
                    ui.label(format!("Oligo length: {} bp", length));
                });

                ui.separator();

                if pos_result.analysis.skipped {
                    ui.colored_label(
                        egui::Color32::YELLOW,
                        format!(
                            "This window was skipped: {}",
                            pos_result
                                .analysis
                                .skip_reason
                                .as_deref()
                                .unwrap_or("Unknown reason")
                        ),
                    );
                    return;
                }

                ui.label(format!(
                    "Total sequences: {}",
                    pos_result.analysis.total_sequences
                ));
                ui.label(format!(
                    "Sequences analyzed: {}",
                    pos_result.analysis.sequences_analyzed
                ));
                ui.label(format!(
                    "Variants needed for {:.0}% coverage: {}",
                    coverage_threshold, pos_result.variants_needed
                ));
                ui.label(format!(
                    "Coverage at threshold: {:.1}%",
                    pos_result.analysis.coverage_at_threshold
                ));

                ui.separator();

                // Display options
                ui.horizontal(|ui| {
                    ui.heading("Variants");
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        ui.checkbox(&mut self.detail_show_codon_spacing, "Codon spacing");
                        ui.checkbox(
                            &mut self.detail_show_reverse_complement,
                            "Reverse complement",
                        );
                    });
                });

                egui::ScrollArea::vertical()
                    .max_height(300.0)
                    .show(ui, |ui| {
                        egui::Grid::new("variants_grid")
                            .striped(true)
                            .min_col_width(50.0)
                            .show(ui, |ui| {
                                ui.strong("#");
                                ui.strong("Sequence");
                                ui.strong("Count");
                                ui.strong("Percentage");
                                ui.strong("Cumulative");
                                ui.end_row();

                                let mut cumulative = 0.0;
                                for (i, variant) in
                                    pos_result.analysis.variants.iter().enumerate()
                                {
                                    cumulative += variant.percentage;

                                    let is_threshold = i + 1 == pos_result.variants_needed;

                                    if is_threshold {
                                        ui.colored_label(
                                            egui::Color32::GREEN,
                                            format!("{}", i + 1),
                                        );
                                    } else {
                                        ui.label(format!("{}", i + 1));
                                    }

                                    let display_seq = format_sequence_for_display(
                                        &variant.sequence,
                                        show_reverse_complement,
                                        show_codon_spacing,
                                    );

                                    ui.add(
                                        egui::Label::new(
                                            egui::RichText::new(&display_seq)
                                                .monospace()
                                                .size(11.0),
                                        )
                                        .wrap_mode(egui::TextWrapMode::Extend),
                                    );

                                    ui.label(format!("{}", variant.count));
                                    ui.label(format!("{:.1}%", variant.percentage));

                                    if is_threshold {
                                        ui.colored_label(
                                            egui::Color32::GREEN,
                                            format!("{:.1}%", cumulative),
                                        );
                                    } else {
                                        ui.label(format!("{:.1}%", cumulative));
                                    }

                                    ui.end_row();
                                }
                            });
                    });
            });
    }
}

/// Format a sequence for display with optional transformations
fn format_sequence_for_display(seq: &str, reverse_comp: bool, codon_spacing: bool) -> String {
    let mut result = if reverse_comp {
        reverse_complement(seq)
    } else {
        seq.to_string()
    };

    if codon_spacing {
        result = add_codon_spacing(&result);
    }

    result
}

/// Add spaces every 3 characters (codon format)
fn add_codon_spacing(seq: &str) -> String {
    seq.chars()
        .enumerate()
        .flat_map(|(i, c)| {
            if i > 0 && i % 3 == 0 {
                vec![' ', c]
            } else {
                vec![c]
            }
        })
        .collect()
}

/// Get color for variant count (green = good, red = bad)
fn variant_count_color(count: usize) -> egui::Color32 {
    if count == 0 {
        return egui::Color32::DARK_GRAY;
    }

    if count == 1 {
        egui::Color32::from_rgb(0, 180, 0) // Green - best
    } else if count <= 5 {
        egui::Color32::from_rgb(150, 180, 0) // Yellow-green
    } else if count <= 10 {
        egui::Color32::from_rgb(255, 165, 0) // Orange
    } else {
        egui::Color32::from_rgb(220, 50, 50) // Red - worst
    }
}
