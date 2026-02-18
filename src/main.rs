mod analysis;
mod app;

use app::OligoscreenApp;

fn main() -> eframe::Result<()> {
    let native_options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1100.0, 800.0])
            .with_min_inner_size([800.0, 600.0]),
        ..Default::default()
    };

    eframe::run_native(
        "Oligoscreen Pairwise",
        native_options,
        Box::new(|cc| Ok(Box::new(OligoscreenApp::new(cc)))),
    )
}
