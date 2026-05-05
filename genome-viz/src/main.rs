//! Terminal UI application for visualizing DNA sequences.

use std::{
    io::{self, stdout},
    path::PathBuf,
    time::Duration,
};

use anyhow::{Context, Result};
use crossterm::{
    ExecutableCommand,
    event::{self, Event, KeyEventKind},
    terminal::{EnterAlternateScreen, LeaveAlternateScreen, disable_raw_mode, enable_raw_mode},
};
use ratatui::{Terminal, backend::CrosstermBackend};

use crate::app::AppState;
use crate::display::draw_ui;

mod app;
mod complement;
mod display;
mod reader;

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: {} <fasta_file> [index_file]", args[0]);
        std::process::exit(1);
    }

    let fasta_path = PathBuf::from(&args[1]);
    let index_path = args.get(2).map(PathBuf::from);

    // Install panic hook to restore terminal on crash.
    let original_hook = std::panic::take_hook();
    std::panic::set_hook(Box::new(move |info| {
        let _ = disable_raw_mode();
        let _ = stdout().execute(LeaveAlternateScreen);
        original_hook(info);
    }));

    // Enter alternate screen and raw mode.
    stdout().execute(EnterAlternateScreen)?;
    enable_raw_mode()?;

    let mut terminal = Terminal::new(CrosstermBackend::new(stdout()))?;
    let mut app = AppState::new(&fasta_path, index_path.as_deref())
        .context("failed to initialize application")?;

    // Set initial size.
    let size = terminal.size()?;
    app.resize(size.width, size.height);

    let result = run_app(&mut terminal, &mut app);

    // Restore terminal.
    disable_raw_mode()?;
    stdout().execute(LeaveAlternateScreen)?;

    result
}

fn run_app(
    terminal: &mut Terminal<CrosstermBackend<io::Stdout>>,
    app: &mut AppState,
) -> Result<()> {
    while app.running {
        // Fetch current window for rendering.
        let bases = app.fetch_current_window();

        terminal.draw(|f| {
            draw_ui(f, app, &bases);
        })?;

        // Poll for events with a timeout so resize is handled promptly.
        if event::poll(Duration::from_millis(100))? {
            match event::read()? {
                Event::Key(key) if key.kind == KeyEventKind::Press => {
                    app.on_key(key);
                }
                Event::Resize(width, height) => {
                    app.resize(width, height);
                }
                _ => {}
            }
        }
    }

    Ok(())
}
