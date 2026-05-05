//! Application state and key-event handling.

use std::path::{Path, PathBuf};

use anyhow::{Context, Result};

use crate::reader::{FastaReader, IndexRecord, default_index_path, read_index};

/// Main application state.
pub struct AppState {
    pub records: Vec<IndexRecord>,
    pub current_index: usize,
    pub current_pos: u64,
    pub scroll_step: u64,
    pub page_size: u64,
    pub view_width: usize,
    pub running: bool,
    pub error_message: Option<String>,
    reader: FastaReader,
}

impl AppState {
    /// Initialize the application: load index, open FASTA, set defaults.
    pub fn new(fasta_path: &Path, index_path: Option<&Path>) -> Result<Self> {
        let index_path = index_path
            .map(PathBuf::from)
            .unwrap_or_else(|| default_index_path(fasta_path));

        let records = read_index(&index_path)
            .with_context(|| format!("failed to read index {}", index_path.display()))?;

        if records.is_empty() {
            anyhow::bail!("index file is empty");
        }

        let reader = FastaReader::open(fasta_path)?;

        Ok(Self {
            records,
            current_index: 0,
            current_pos: 1,
            scroll_step: 10,
            page_size: 20,
            view_width: 40,
            running: true,
            error_message: None,
            reader,
        })
    }

    /// Resize-aware recalculation of page and view dimensions.
    pub fn resize(&mut self, width: u16, height: u16) {
        // Reserve 2 lines for header and footer
        let body_height = height.saturating_sub(2).max(1);
        self.page_size = body_height as u64;

        // view_width is kept for potential future multi-base-per-line support.
        // Currently each TUI line shows one base, so this is mainly padding calc.
        self.view_width = (width as usize).saturating_sub(22).max(1);
    }

    /// Handle a key press and update state.
    pub fn on_key(&mut self, key: crossterm::event::KeyEvent) {
        use crossterm::event::{KeyCode, KeyModifiers};

        match key.code {
            KeyCode::Up | KeyCode::Char('k') => self.scroll_up(self.scroll_step),
            KeyCode::Down | KeyCode::Char('j') => self.scroll_down(self.scroll_step),
            KeyCode::PageUp => self.scroll_up(self.page_size),
            KeyCode::PageDown => self.scroll_down(self.page_size),
            KeyCode::Char('n') => self.next_chromosome(),
            KeyCode::Char('p') => self.prev_chromosome(),
            KeyCode::Char('m') => self.jump_to_mitochondria(),
            KeyCode::Char('g') => self.jump_to_start(),
            KeyCode::Char('G') => self.jump_to_end(),
            KeyCode::Char('q') | KeyCode::Char('Q') => self.running = false,
            KeyCode::Char('c') if key.modifiers.contains(KeyModifiers::CONTROL) => {
                self.running = false;
            }
            _ => {}
        }
    }

    /// Scroll up by the given amount, clamped to the start of the chromosome.
    pub fn scroll_up(&mut self, amount: u64) {
        self.current_pos = self.current_pos.saturating_sub(amount).max(1);
    }

    /// Scroll down by the given amount, clamped to the chromosome length.
    pub fn scroll_down(&mut self, amount: u64) {
        let record = &self.records[self.current_index];
        self.current_pos = (self.current_pos + amount).min(record.length);
    }

    /// Move to the next chromosome, wrapping around.
    pub fn next_chromosome(&mut self) {
        self.current_index = (self.current_index + 1) % self.records.len();
        self.current_pos = 1;
    }

    /// Move to the previous chromosome, wrapping around.
    pub fn prev_chromosome(&mut self) {
        if self.current_index == 0 {
            self.current_index = self.records.len() - 1;
        } else {
            self.current_index -= 1;
        }
        self.current_pos = 1;
    }

    /// Jump to the mitochondria contig if present.
    pub fn jump_to_mitochondria(&mut self) {
        if let Some(idx) = self
            .records
            .iter()
            .position(|r| matches!(r.name.as_str(), "chrM" | "MT" | "M"))
        {
            self.current_index = idx;
            self.current_pos = 1;
        }
    }

    /// Jump to the start of the current chromosome.
    pub fn jump_to_start(&mut self) {
        self.current_pos = 1;
    }

    /// Jump to the end of the current chromosome.
    pub fn jump_to_end(&mut self) {
        let record = &self.records[self.current_index];
        self.current_pos = record.length;
    }

    /// Fetch the current window of bases for display.
    pub fn fetch_current_window(&mut self) -> Vec<u8> {
        let record = &self.records[self.current_index];
        match self
            .reader
            .fetch_bases(record, self.current_pos, self.page_size)
        {
            Ok(bases) => {
                self.error_message = None;
                bases
            }
            Err(e) => {
                self.error_message = Some(format!("{e}"));
                Vec::new()
            }
        }
    }

    /// Return the length of the current contig.
    pub fn current_length(&self) -> u64 {
        self.records[self.current_index].length
    }
}
