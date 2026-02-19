//! Oligoscreen Tool - Primer Site Screening Tool
//!
//! A Rust application for screening large sets of aligned DNA sequences
//! to find suitable primer sites with low variability.

mod analysis;
mod app;

use app::OligoscreenApp;

fn main() -> eframe::Result<()> {
    let native_options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1200.0, 800.0])
            .with_min_inner_size([900.0, 600.0])
            .with_title("Oligoscreen Tool"),
        ..Default::default()
    };

    eframe::run_native(
        "Oligoscreen Tool",
        native_options,
        Box::new(|cc| Ok(Box::new(OligoscreenApp::new(cc)))),
    )
}
