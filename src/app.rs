//! Main application state and UI

use eframe::egui;
use std::sync::mpsc::{channel, Receiver};
use std::thread;

use crate::analysis::{
    parse_reference_fasta, parse_template_fasta, reverse_complement, run_screening,
    AnalysisMethod, AnalysisParams, ProgressUpdate, ReferenceData, ScreeningResults, TemplateData,
    ThreadCount,
};

/// Application state
pub struct OligoscreenApp {
    // Input tab state - template
    template_input: String,
    template_data: Option<TemplateData>,
    template_error: Option<String>,

    // Input tab state - references
    reference_input: String,
    reference_data: Option<ReferenceData>,
    reference_error: Option<String>,

    // Analysis parameters
    params: AnalysisParams,
    method_selection: MethodSelection,
    thread_selection: ThreadSelection,
    manual_thread_count: usize,

    // Incremental method options
    incremental_limit_ambiguities: bool,
    incremental_max_ambiguities: u32,

    // Analysis state
    is_analyzing: bool,
    analysis_progress: Option<ProgressUpdate>,
    progress_rx: Option<Receiver<ProgressUpdate>>,
    results_rx: Option<Receiver<ScreeningResults>>,

    // Results state
    results: Option<ScreeningResults>,
    selected_position: Option<usize>,
    selected_length_for_detail: Option<u32>,
    show_detail_window: bool,

    // Detail window display options
    detail_show_reverse_complement: bool,
    detail_show_codon_spacing: bool,

    // View state
    current_tab: Tab,
    zoom_level: f32,

    // Save/Load
    save_error: Option<String>,
    load_error: Option<String>,

    // Deferred actions
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
            template_data: None,
            template_error: None,
            reference_input: String::new(),
            reference_data: None,
            reference_error: None,
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
            selected_position: None,
            selected_length_for_detail: None,
            show_detail_window: false,
            detail_show_reverse_complement: false,
            detail_show_codon_spacing: true,
            current_tab: Tab::Input,
            zoom_level: 1.0,
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

    fn parse_template_input(&mut self) {
        self.template_error = None;
        self.template_data = None;

        if self.template_input.trim().is_empty() {
            return;
        }

        match parse_template_fasta(&self.template_input) {
            Ok(data) => {
                self.template_data = Some(data);
            }
            Err(e) => {
                self.template_error = Some(e);
            }
        }
    }

    fn parse_reference_input(&mut self) {
        self.reference_error = None;
        self.reference_data = None;

        if self.reference_input.trim().is_empty() {
            return;
        }

        match parse_reference_fasta(&self.reference_input) {
            Ok(data) => {
                self.reference_data = Some(data);
            }
            Err(e) => {
                self.reference_error = Some(e);
            }
        }
    }

    fn start_analysis(&mut self) {
        let Some(template) = &self.template_data else {
            return;
        };
        let Some(references) = &self.reference_data else {
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

        let template_clone = template.clone();
        let references_clone = references.clone();
        let params_clone = self.params.clone();

        let (progress_tx, progress_rx) = channel();
        let (results_tx, results_rx) = channel();

        self.progress_rx = Some(progress_rx);
        self.results_rx = Some(results_rx);
        self.is_analyzing = true;
        self.analysis_progress = None;

        thread::spawn(move || {
            let results = run_screening(
                &template_clone,
                &references_clone,
                &params_clone,
                Some(progress_tx),
            );
            let _ = results_tx.send(results);
        });
    }

    fn check_analysis_progress(&mut self) {
        if let Some(rx) = &self.progress_rx {
            while let Ok(progress) = rx.try_recv() {
                self.analysis_progress = Some(progress);
            }
        }

        if let Some(rx) = &self.results_rx {
            if let Ok(results) = rx.try_recv() {
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
                    self.parse_template_input();
                }
                Err(e) => {
                    self.template_error = Some(format!("Failed to read file: {}", e));
                }
            }
        }
    }

    fn load_reference_file(&mut self) {
        if let Some(path) = rfd::FileDialog::new()
            .add_filter("FASTA", &["fasta", "fa", "fna", "fas", "txt"])
            .pick_file()
        {
            match std::fs::read_to_string(&path) {
                Ok(content) => {
                    self.reference_input = content;
                    self.parse_reference_input();
                }
                Err(e) => {
                    self.reference_error = Some(format!("Failed to read file: {}", e));
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
        if self.is_analyzing {
            self.check_analysis_progress();
            ctx.request_repaint();
        }

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
                    if ui.button("Load References...").clicked() {
                        self.load_reference_file();
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
                if self.is_analyzing {
                    ui.spinner();
                    if let Some(ref progress) = self.analysis_progress {
                        ui.label(&progress.message);
                    } else {
                        ui.label("Starting analysis...");
                    }
                } else if let Some(ref results) = self.results {
                    ui.label(format!(
                        "Results: {} references, {} bp template",
                        results.total_sequences, results.template_length
                    ));
                } else {
                    let mut parts = Vec::new();
                    if let Some(ref t) = self.template_data {
                        parts.push(format!("Template: {} bp", t.sequence.len()));
                    }
                    if let Some(ref r) = self.reference_data {
                        parts.push(format!("References: {} sequences", r.len()));
                    }
                    if parts.is_empty() {
                        ui.label("Load template and reference sequences to begin");
                    } else {
                        ui.label(parts.join(" | "));
                    }
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

        // Use available height for two panels
        let available_height = ui.available_height();
        let panel_height = (available_height / 2.0 - 60.0).max(120.0);

        // --- Template Sequence ---
        ui.group(|ui| {
            ui.horizontal(|ui| {
                ui.heading("Template Sequence");
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    if ui.button("Clear").clicked() {
                        self.template_input.clear();
                        self.template_data = None;
                        self.template_error = None;
                    }
                    if ui.button("Load File").clicked() {
                        self.load_template_file();
                    }
                    if ui.button("Load Example").clicked() {
                        self.template_input = EXAMPLE_TEMPLATE.to_string();
                        self.parse_template_input();
                    }
                });
            });

            ui.label("Single sequence in FASTA format (A, C, G, T only):");

            egui::ScrollArea::vertical()
                .id_salt("template_scroll")
                .max_height(panel_height)
                .show(ui, |ui| {
                    let response = ui.add(
                        egui::TextEdit::multiline(&mut self.template_input)
                            .font(egui::TextStyle::Monospace)
                            .desired_width(f32::INFINITY)
                            .desired_rows(6),
                    );
                    if response.changed() {
                        self.parse_template_input();
                    }
                });

            if let Some(ref error) = self.template_error {
                ui.colored_label(egui::Color32::RED, format!("Error: {}", error));
            }
            if let Some(ref data) = self.template_data {
                ui.colored_label(
                    egui::Color32::from_rgb(100, 200, 100),
                    format!("Template: {} ({} bp)", data.name, data.sequence.len()),
                );
            }
        });

        ui.add_space(5.0);

        // --- Reference Sequences ---
        ui.group(|ui| {
            ui.horizontal(|ui| {
                ui.heading("Reference Sequences");
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    if ui.button("Clear").clicked() {
                        self.reference_input.clear();
                        self.reference_data = None;
                        self.reference_error = None;
                    }
                    if ui.button("Load File").clicked() {
                        self.load_reference_file();
                    }
                    if ui.button("Load Example").clicked() {
                        self.reference_input = EXAMPLE_REFERENCES.to_string();
                        self.parse_reference_input();
                    }
                });
            });

            ui.label("Multiple sequences in FASTA format (unaligned):");

            egui::ScrollArea::vertical()
                .id_salt("reference_scroll")
                .max_height(panel_height)
                .show(ui, |ui| {
                    let response = ui.add(
                        egui::TextEdit::multiline(&mut self.reference_input)
                            .font(egui::TextStyle::Monospace)
                            .desired_width(f32::INFINITY)
                            .desired_rows(6),
                    );
                    if response.changed() {
                        self.parse_reference_input();
                    }
                });

            if let Some(ref error) = self.reference_error {
                ui.colored_label(egui::Color32::RED, format!("Error: {}", error));
            }
            if let Some(ref data) = self.reference_data {
                let min_len = data.sequences.iter().map(|s| s.len()).min().unwrap_or(0);
                let max_len = data.sequences.iter().map(|s| s.len()).max().unwrap_or(0);
                ui.colored_label(
                    egui::Color32::from_rgb(100, 200, 100),
                    format!(
                        "References: {} sequences ({}-{} bp)",
                        data.len(),
                        min_len,
                        max_len
                    ),
                );
            }
        });
    }

    fn show_analysis_tab(&mut self, ui: &mut egui::Ui) {
        ui.heading("Analysis Setup");
        ui.separator();

        // Check if data is loaded
        let has_template = self.template_data.is_some();
        let has_references = self.reference_data.is_some();

        if !has_template || !has_references {
            ui.colored_label(
                egui::Color32::YELLOW,
                if !has_template && !has_references {
                    "Please load template and reference sequences in the Input tab."
                } else if !has_template {
                    "Please load a template sequence in the Input tab."
                } else {
                    "Please load reference sequences in the Input tab."
                },
            );
            return;
        }

        egui::ScrollArea::vertical().show(ui, |ui| {
            // Pairwise Aligner Settings
            ui.group(|ui| {
                ui.heading("Pairwise Aligner Settings");

                ui.horizontal(|ui| {
                    ui.label("Match score:");
                    ui.add(egui::DragValue::new(&mut self.params.pairwise.match_score).range(-10..=10));
                    ui.add_space(20.0);
                    ui.label("Mismatch score:");
                    ui.add(egui::DragValue::new(&mut self.params.pairwise.mismatch_score).range(-10..=10));
                });

                ui.horizontal(|ui| {
                    ui.label("Gap open penalty:");
                    ui.add(egui::DragValue::new(&mut self.params.pairwise.gap_open_penalty).range(-20..=0));
                    ui.add_space(20.0);
                    ui.label("Gap extend penalty:");
                    ui.add(egui::DragValue::new(&mut self.params.pairwise.gap_extend_penalty).range(-20..=0));
                });

                ui.horizontal(|ui| {
                    ui.label("Maximum allowed mismatches:");
                    ui.add(egui::DragValue::new(&mut self.params.pairwise.max_mismatches).range(0..=50));
                });
                ui.label("Matches exceeding this mismatch count are recorded as 'no match'.");
            });

            ui.add_space(10.0);

            // Analysis method selection
            ui.group(|ui| {
                ui.heading("Analysis Method");

                ui.radio_value(
                    &mut self.method_selection,
                    MethodSelection::NoAmbiguities,
                    "No Ambiguities - Find all unique exact variants",
                );

                ui.horizontal(|ui| {
                    ui.radio_value(
                        &mut self.method_selection,
                        MethodSelection::FixedAmbiguities,
                        "Fixed Ambiguities - Use up to N ambiguity codes per variant",
                    );
                });

                if self.method_selection == MethodSelection::FixedAmbiguities {
                    ui.horizontal(|ui| {
                        ui.add_space(20.0);
                        ui.label("Max ambiguities:");
                        let mut n = self.params.method.get_fixed_ambiguities();
                        if ui.add(egui::DragValue::new(&mut n).range(0..=20)).changed() {
                            self.params.method = AnalysisMethod::FixedAmbiguities(n);
                        }
                    });
                }

                ui.horizontal(|ui| {
                    ui.radio_value(
                        &mut self.method_selection,
                        MethodSelection::Incremental,
                        "Incremental - Find variants covering X% of remaining sequences",
                    );
                });

                if self.method_selection == MethodSelection::Incremental {
                    ui.horizontal(|ui| {
                        ui.add_space(20.0);
                        ui.label("Target coverage per step (%):");
                        let mut pct = self.params.method.get_incremental_pct();
                        let max_amb = self.params.method.get_incremental_max_amb();
                        if ui.add(egui::DragValue::new(&mut pct).range(1..=100)).changed() {
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
                    if self.incremental_limit_ambiguities {
                        ui.horizontal(|ui| {
                            ui.add_space(20.0);
                            ui.label(
                                "If target % cannot be reached, accepts best variant within limit.",
                            );
                        });
                    }
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

            // Thread count
            ui.group(|ui| {
                ui.heading("Parallelization");

                let available_threads = std::thread::available_parallelism()
                    .map(|n| n.get())
                    .unwrap_or(1);

                ui.horizontal(|ui| {
                    ui.radio_value(
                        &mut self.thread_selection,
                        ThreadSelection::Auto,
                        format!("Auto ({} threads)", available_threads),
                    );
                });

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

                ui.label("More threads = faster analysis but higher CPU usage");
            });

            ui.add_space(20.0);

            // Run button
            ui.horizontal(|ui| {
                let can_run = has_template && has_references && !self.is_analyzing;
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

        // Extract data we need
        let (lengths, template_seq, total_seqs, coverage_threshold) = {
            let results = self.results.as_ref().unwrap();
            let mut lengths: Vec<u32> = results.results_by_length.keys().copied().collect();
            lengths.sort();
            (
                lengths,
                results.template_sequence.clone(),
                results.total_sequences,
                results.params.coverage_threshold,
            )
        };

        if lengths.is_empty() {
            ui.label("No length results available.");
            return;
        }

        // Controls
        ui.horizontal(|ui| {
            ui.label("Zoom:");
            ui.add(egui::Slider::new(&mut self.zoom_level, 0.5..=3.0));
            ui.add_space(20.0);
            ui.label(format!(
                "{} reference sequences | Template: {} bp",
                total_seqs,
                template_seq.len()
            ));
        });

        ui.add_space(5.0);

        // Heatmap display
        self.show_heatmap(ui, &lengths, &template_seq, coverage_threshold);

        // Error messages
        if let Some(ref error) = self.save_error {
            ui.colored_label(egui::Color32::RED, error);
        }
        if let Some(ref error) = self.load_error {
            ui.colored_label(egui::Color32::RED, error);
        }
    }

    fn show_heatmap(
        &mut self,
        ui: &mut egui::Ui,
        lengths: &[u32],
        template_seq: &str,
        coverage_threshold: f64,
    ) {
        let results = self.results.as_ref().unwrap();

        // Get positions from the first length result
        let first_length_result = results.results_by_length.get(&lengths[0]);
        let positions: Vec<usize> = first_length_result
            .map(|lr| lr.positions.iter().map(|p| p.position).collect())
            .unwrap_or_default();

        if positions.is_empty() {
            ui.label("No positions analyzed.");
            return;
        }

        // Cell dimensions
        let cell_w = (14.0 * self.zoom_level).max(6.0);
        let cell_h = (18.0 * self.zoom_level).max(10.0);
        let label_width = 50.0; // Width for row labels
        let header_height = 30.0 * self.zoom_level; // Height for template sequence row
        let pos_label_height = 16.0 * self.zoom_level;

        let num_cols = positions.len();
        let num_rows = lengths.len();

        // Summary stats per length
        ui.group(|ui| {
            ui.horizontal_wrapped(|ui| {
                for &length in lengths {
                    if let Some(lr) = results.results_by_length.get(&length) {
                        let non_skipped: Vec<_> =
                            lr.positions.iter().filter(|p| !p.analysis.skipped).collect();
                        if !non_skipped.is_empty() {
                            let avg: f64 =
                                non_skipped.iter().map(|p| p.variants_needed).sum::<usize>()
                                    as f64
                                    / non_skipped.len() as f64;
                            let min = non_skipped
                                .iter()
                                .map(|p| p.variants_needed)
                                .min()
                                .unwrap_or(0);
                            let max = non_skipped
                                .iter()
                                .map(|p| p.variants_needed)
                                .max()
                                .unwrap_or(0);
                            ui.label(format!(
                                "{}bp: {}-{} (avg {:.1})",
                                length, min, max, avg
                            ));
                            ui.separator();
                        }
                    }
                }
            });
        });

        ui.add_space(5.0);

        ui.label(format!(
            "Variants needed to reach {:.0}% coverage (click cell for details):",
            coverage_threshold
        ));

        // Build heatmap data: lookup by (length, position) -> variants_needed
        let heatmap_data: std::collections::HashMap<(u32, usize), &crate::analysis::PositionResult> =
            {
                let mut map = std::collections::HashMap::new();
                for &length in lengths {
                    if let Some(lr) = results.results_by_length.get(&length) {
                        for pr in &lr.positions {
                            map.insert((length, pr.position), pr);
                        }
                    }
                }
                map
            };

        // Total width/height for the heatmap area
        let total_width = label_width + (num_cols as f32 * cell_w);
        let total_height =
            pos_label_height + header_height + (num_rows as f32 * cell_h) + 30.0; // +30 for legend

        egui::ScrollArea::horizontal()
            .id_salt("heatmap_scroll")
            .show(ui, |ui| {
                let (response, painter) = ui.allocate_painter(
                    egui::vec2(total_width, total_height),
                    egui::Sense::click_and_drag(),
                );
                let origin = response.rect.min;

                // --- Position numbers row ---
                let show_every_n = if cell_w < 12.0 {
                    (12.0 / cell_w).ceil() as usize
                } else {
                    1
                };

                for (col, &pos) in positions.iter().enumerate() {
                    if col % show_every_n != 0 {
                        continue;
                    }
                    let x = origin.x + label_width + (col as f32 * cell_w) + cell_w / 2.0;
                    let y = origin.y + pos_label_height / 2.0;
                    painter.text(
                        egui::pos2(x, y),
                        egui::Align2::CENTER_CENTER,
                        format!("{}", pos + 1),
                        egui::FontId::proportional(8.0 * self.zoom_level.min(1.5)),
                        egui::Color32::GRAY,
                    );
                }

                // --- Template sequence row ---
                let seq_y_start = origin.y + pos_label_height;
                for (col, &pos) in positions.iter().enumerate() {
                    if pos < template_seq.len() {
                        let base = &template_seq[pos..pos + 1];
                        let x = origin.x + label_width + (col as f32 * cell_w) + cell_w / 2.0;
                        let y = seq_y_start + header_height / 2.0;

                        let color = base_color(base.chars().next().unwrap_or('N'));
                        painter.text(
                            egui::pos2(x, y),
                            egui::Align2::CENTER_CENTER,
                            base,
                            egui::FontId::monospace(10.0 * self.zoom_level.min(1.5)),
                            color,
                        );
                    }
                }

                // --- Row labels (oligo lengths) ---
                let grid_y_start = seq_y_start + header_height;
                for (row, &length) in lengths.iter().enumerate() {
                    let y = grid_y_start + (row as f32 * cell_h) + cell_h / 2.0;
                    painter.text(
                        egui::pos2(origin.x + label_width - 5.0, y),
                        egui::Align2::RIGHT_CENTER,
                        format!("{} bp", length),
                        egui::FontId::proportional(10.0 * self.zoom_level.min(1.5)),
                        egui::Color32::LIGHT_GRAY,
                    );
                }

                // --- Heatmap cells ---
                let mut hovered_cell: Option<(u32, usize)> = None;
                let mut clicked_cell: Option<(u32, usize)> = None;

                for (row, &length) in lengths.iter().enumerate() {
                    for (col, &pos) in positions.iter().enumerate() {
                        let cell_x = origin.x + label_width + (col as f32 * cell_w);
                        let cell_y = grid_y_start + (row as f32 * cell_h);
                        let cell_rect = egui::Rect::from_min_size(
                            egui::pos2(cell_x, cell_y),
                            egui::vec2(cell_w - 1.0, cell_h - 1.0),
                        );

                        let color = if let Some(pr) = heatmap_data.get(&(length, pos)) {
                            if pr.analysis.skipped {
                                egui::Color32::from_rgb(40, 40, 40)
                            } else {
                                let no_match_frac = if pr.analysis.total_sequences > 0 {
                                    pr.analysis.no_match_count as f64
                                        / pr.analysis.total_sequences as f64
                                } else {
                                    0.0
                                };
                                position_color(pr.variants_needed, no_match_frac)
                            }
                        } else {
                            egui::Color32::from_rgb(30, 30, 30)
                        };

                        painter.rect_filled(cell_rect, 1.0, color);

                        // Check hover/click using the response's pointer
                        if let Some(pointer_pos) = response.hover_pos() {
                            if cell_rect.contains(pointer_pos) {
                                hovered_cell = Some((length, pos));
                                painter.rect_stroke(
                                    cell_rect,
                                    1.0,
                                    egui::Stroke::new(1.5, egui::Color32::WHITE),
                                    egui::StrokeKind::Outside,
                                );
                            }
                        }

                        if response.clicked() {
                            if let Some(pointer_pos) = ui.ctx().pointer_latest_pos() {
                                if cell_rect.contains(pointer_pos) {
                                    clicked_cell = Some((length, pos));
                                }
                            }
                        }
                    }
                }

                // Handle tooltip
                if let Some((length, pos)) = hovered_cell {
                    if let Some(pr) = heatmap_data.get(&(length, pos)) {
                        let tooltip_text = if pr.analysis.skipped {
                            format!(
                                "Position: {}, Length: {} bp\nSkipped: {}",
                                pos + 1,
                                length,
                                pr.analysis
                                    .skip_reason
                                    .as_deref()
                                    .unwrap_or("Unknown")
                            )
                        } else {
                            format!(
                                "Position: {}, Length: {} bp\nVariants needed: {}\nCoverage: {:.1}%\nMatched: {}/{}\nNo match: {}",
                                pos + 1,
                                length,
                                pr.variants_needed,
                                pr.analysis.coverage_at_threshold,
                                pr.analysis.sequences_analyzed,
                                pr.analysis.total_sequences,
                                pr.analysis.no_match_count,
                            )
                        };
                        response.clone().on_hover_text(tooltip_text);
                    }
                }

                // Handle click
                if let Some((length, pos)) = clicked_cell {
                    self.selected_position = Some(pos);
                    self.selected_length_for_detail = Some(length);
                    self.show_detail_window = true;
                }
            });

        // Legend
        ui.add_space(5.0);
        ui.horizontal(|ui| {
            ui.label("Legend:");
            ui.add_space(10.0);

            let legend_colors = [
                (egui::Color32::from_rgb(0, 180, 0), "1 (best)"),
                (egui::Color32::from_rgb(150, 180, 0), "2-5"),
                (egui::Color32::from_rgb(255, 165, 0), "6-10"),
                (egui::Color32::from_rgb(220, 50, 50), ">10"),
                (egui::Color32::from_rgb(40, 40, 40), "skipped/no data"),
            ];

            for (color, label) in legend_colors {
                let (rect, _) =
                    ui.allocate_exact_size(egui::vec2(15.0, 15.0), egui::Sense::hover());
                ui.painter().rect_filled(rect, 2.0, color);
                ui.label(label);
                ui.add_space(10.0);
            }
        });
    }

    fn show_variant_detail_window(&mut self, ctx: &egui::Context) {
        let Some(ref results) = self.results else {
            self.show_detail_window = false;
            return;
        };

        let Some(length) = self.selected_length_for_detail else {
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

        // Extract template oligo for display
        let template_oligo = if position + length as usize <= results.template_sequence.len() {
            &results.template_sequence[position..position + length as usize]
        } else {
            ""
        };
        let template_oligo = template_oligo.to_string();

        let show_reverse_complement = self.detail_show_reverse_complement;
        let show_codon_spacing = self.detail_show_codon_spacing;

        egui::Window::new(format!("Position {} Details", position + 1))
            .open(&mut self.show_detail_window)
            .default_width(650.0)
            .default_height(450.0)
            .show(ctx, |ui| {
                ui.horizontal(|ui| {
                    ui.label(format!("Position: {}", position + 1));
                    ui.separator();
                    ui.label(format!("Oligo length: {} bp", length));
                });

                // Template oligo display
                if !template_oligo.is_empty() {
                    let display_template = format_sequence_for_display(
                        &template_oligo,
                        show_reverse_complement,
                        show_codon_spacing,
                    );
                    ui.horizontal(|ui| {
                        ui.label("Template oligo:");
                        ui.add(
                            egui::Label::new(
                                egui::RichText::new(&display_template)
                                    .monospace()
                                    .size(11.0)
                                    .color(egui::Color32::from_rgb(100, 180, 255)),
                            )
                            .wrap_mode(egui::TextWrapMode::Extend),
                        );
                    });
                }

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
                    "Total references: {}",
                    pos_result.analysis.total_sequences
                ));
                ui.label(format!(
                    "Matched: {}",
                    pos_result.analysis.sequences_analyzed
                ));
                if pos_result.analysis.no_match_count > 0 {
                    ui.colored_label(
                        egui::Color32::from_rgb(255, 180, 100),
                        format!(
                            "No match: {}/{} ({:.1}%)",
                            pos_result.analysis.no_match_count,
                            pos_result.analysis.total_sequences,
                            (pos_result.analysis.no_match_count as f64
                                / pos_result.analysis.total_sequences as f64)
                                * 100.0
                        ),
                    );
                }
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

                                // No match row
                                if pos_result.analysis.no_match_count > 0 {
                                    ui.label("");
                                    ui.colored_label(
                                        egui::Color32::from_rgb(255, 180, 100),
                                        "No match",
                                    );
                                    ui.colored_label(
                                        egui::Color32::from_rgb(255, 180, 100),
                                        format!("{}", pos_result.analysis.no_match_count),
                                    );
                                    let no_match_pct = (pos_result.analysis.no_match_count as f64
                                        / pos_result.analysis.total_sequences as f64)
                                        * 100.0;
                                    ui.colored_label(
                                        egui::Color32::from_rgb(255, 180, 100),
                                        format!("{:.1}%", no_match_pct),
                                    );
                                    ui.label("");
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

/// Get color for a position based on variant count and no-match fraction.
///
/// Green only when both variant count is low AND no-match fraction is low.
/// High no-match fraction shifts the color toward red regardless of variant count.
fn position_color(variant_count: usize, no_match_fraction: f64) -> egui::Color32 {
    if variant_count == 0 {
        return egui::Color32::from_rgb(40, 40, 40);
    }

    // Base color from variant count (same scale as before)
    let (base_r, base_g, base_b) = if variant_count == 1 {
        (0u8, 180u8, 0u8) // Green - best
    } else if variant_count <= 5 {
        (150, 180, 0) // Yellow-green
    } else if variant_count <= 10 {
        (255, 165, 0) // Orange
    } else {
        (220, 50, 50) // Red - worst
    };

    // Blend toward red based on no-match fraction
    // At 0% no-match: use base color as-is
    // At 100% no-match: fully red
    let nm = no_match_fraction.clamp(0.0, 1.0);
    let red_r: u8 = 220;
    let red_g: u8 = 50;
    let red_b: u8 = 50;

    let r = (base_r as f64 * (1.0 - nm) + red_r as f64 * nm) as u8;
    let g = (base_g as f64 * (1.0 - nm) + red_g as f64 * nm) as u8;
    let b = (base_b as f64 * (1.0 - nm) + red_b as f64 * nm) as u8;

    egui::Color32::from_rgb(r, g, b)
}

/// Color for DNA base letters in the template display
fn base_color(base: char) -> egui::Color32 {
    match base {
        'A' => egui::Color32::from_rgb(100, 200, 100), // Green
        'T' => egui::Color32::from_rgb(220, 80, 80),   // Red
        'G' => egui::Color32::from_rgb(255, 200, 60),   // Yellow/gold
        'C' => egui::Color32::from_rgb(100, 150, 255),  // Blue
        _ => egui::Color32::GRAY,
    }
}

const EXAMPLE_TEMPLATE: &str = r#">Template
TATGGTACGTCATGTTCTAGAAATGGGCTGT
"#;

const EXAMPLE_REFERENCES: &str = r#">Ref1
TATGGTACGTCATGTTCTAGAAATGGGCTGT
>Ref2
AATATGGTACGTCATGTTCTAGAAATGGGCTGT
>Ref3
TATGGTTCGTCATGTTCTAGAAATGGGCTGTTTT
>Ref4
GTATGGTACGTCATGTTCTAGAAATGGGCTGT
"#;
