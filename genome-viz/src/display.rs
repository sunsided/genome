//! TUI rendering for the DNA ladder visualizer.

use ratatui::{
    Frame,
    layout::{Alignment, Constraint, Direction, Layout, Rect},
    style::{Color, Modifier, Style},
    text::{Line, Span, Text},
    widgets::Paragraph,
};

use crate::complement::complement;

/// Render the full application UI into the given frame.
pub fn draw_ui(f: &mut Frame, app: &crate::app::AppState, bases: &[u8]) {
    let chunks = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(1),
            Constraint::Min(1),
            Constraint::Length(1),
        ])
        .split(f.area());

    let visible = bases.len() as u64;
    draw_header(f, app, visible, chunks[0]);
    draw_body(f, app, bases, chunks[1]);
    draw_footer(f, chunks[2]);
}

fn draw_header(f: &mut Frame, app: &crate::app::AppState, visible: u64, area: Rect) {
    let record = &app.records[app.current_index];
    let text = format!(
        " {}:{} | pos {} | {}/{} lines visible ",
        record.name,
        format_number(record.length),
        format_number(app.current_pos),
        visible,
        app.page_size
    );
    let paragraph = Paragraph::new(text)
        .style(Style::default().fg(Color::Cyan))
        .alignment(Alignment::Center);
    f.render_widget(paragraph, area);
}

fn draw_body(f: &mut Frame, app: &crate::app::AppState, bases: &[u8], area: Rect) {
    let pos_width = digit_count(app.current_length());

    let mut lines = Vec::with_capacity(bases.len());
    for (i, &base) in bases.iter().enumerate() {
        let pos = app.current_pos + i as u64;
        let comp = complement(base);

        let base_st = base_style(base);
        let comp_st = base_style(comp);

        let pos_str = format!("{:>1$}", pos, pos_width);
        let line = Line::from(vec![
            Span::styled(format!("{}  ", pos_str), Style::default().fg(Color::Gray)),
            Span::raw("|-"),
            Span::styled(String::from_utf8_lossy(&[base]).into_owned(), base_st),
            Span::raw(" "),
            Span::styled(String::from_utf8_lossy(&[comp]).into_owned(), comp_st),
            Span::raw("-|"),
        ]);
        lines.push(line);
    }

    if let Some(ref err) = app.error_message {
        lines.push(Line::from(""));
        lines.push(Line::from(Span::styled(
            format!("Error: {}", err),
            Style::default().fg(Color::Red),
        )));
    } else if bases.is_empty() && app.error_message.is_none() {
        lines.push(Line::from(Span::styled(
            "No data to display.",
            Style::default().fg(Color::DarkGray),
        )));
    }

    let paragraph = Paragraph::new(Text::from(lines));
    f.render_widget(paragraph, area);
}

fn draw_footer(f: &mut Frame, area: Rect) {
    let text = " ↑↓/kj:scroll pgup/pgdn:page [n]ext/[p]rev chr [1]chr1 [s]kip-N [m]ito [g]start [G]end [q]uit ";
    let paragraph = Paragraph::new(text)
        .style(Style::default().fg(Color::Gray))
        .alignment(Alignment::Center);
    f.render_widget(paragraph, area);
}

/// Return the number of digits in a positive integer.
fn digit_count(n: u64) -> usize {
    if n == 0 { 1 } else { n.ilog10() as usize + 1 }
}

/// Format a number with comma separators.
fn format_number(n: u64) -> String {
    let s = n.to_string();
    let mut result = String::with_capacity(s.len() + s.len() / 3);
    let first = s.len() % 3;
    let first = if first == 0 { 3 } else { first };

    for (i, c) in s.chars().enumerate() {
        if i > 0 && i == first {
            result.push(',');
        } else if i > first && (i - first) % 3 == 0 {
            result.push(',');
        }
        result.push(c);
    }
    result
}

/// Return the style for a given DNA base.
fn base_style(base: u8) -> Style {
    let is_lower = base.is_ascii_lowercase();
    let upper = base.to_ascii_uppercase();

    let color = match upper {
        b'A' => Color::Green,
        b'C' => Color::Blue,
        b'G' => Color::Yellow,
        b'T' => Color::Red,
        b'N' => Color::DarkGray,
        _ => Color::White,
    };

    let style = Style::default().fg(color);
    if is_lower {
        style.add_modifier(Modifier::DIM)
    } else {
        style
    }
}
