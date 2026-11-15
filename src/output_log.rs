//! The output log, `<output>/output_log.txt`, holds everything CluE prints
//! to the terminal during a run, stdout and stderr, in the order printed.
//!
//! The prints are spread over modules that do not know the output directory,
//! so the log is a process-wide state, `LOG`. The macros `clue_print!`,
//! `clue_println!`, and `clue_eprintln!` print their text as `print!`,
//! `println!`, and `eprintln!` do, and pass the same text to the log:
//! * before the output directory is known, the text is buffered;
//! * `open()` writes the buffer to the log file and appends every later text
//!   to the file directly (unbuffered, so the log of an interrupted run is
//!   complete up to the interruption);
//! * `discard()` drops the buffer and stops logging (output not saved);
//! * `finish()` closes the file and starts a new, empty buffer, so that a
//!   second run in the same process starts a new log.
//!
//! A failing write prints a warning to stderr and stops the log; it never
//! stops a calculation. The log is one per process: runs that execute at the
//! same time in one process write into whichever log was opened last.
use crate::clue_errors::CluEError;

use std::fs::{File,OpenOptions};
use std::io::Write;
use std::sync::{Mutex,MutexGuard};

/// `LogState` is the state of the output log.
/// * `Buffering(text)`: the log is not open yet; `text` is what was printed
///   so far.
/// * `Writing(file, path)`: the text is appended to `file`, opened from
///   `path`.
/// * `Off`: nothing is recorded.
#[derive(Debug)]
pub enum LogState{
  Buffering(String),
  Writing(File,String),
  Off,
}

impl LogState{
  //----------------------------------------------------------------------------
  /// This function records `text`: it appends it to the buffer or to the
  /// file, or does nothing if the log is off. If writing fails, it prints a
  /// warning to stderr and switches the log off.
  pub fn record(&mut self, text: &str){
    match self{
      LogState::Buffering(buffer) => buffer.push_str(text),
      LogState::Writing(file, _) => {
        if let Err(err) = file.write_all(text.as_bytes()){
          self.stop(&err);
        }
      },
      LogState::Off => (),
    }
  }
  //----------------------------------------------------------------------------
  /// This function opens the file `path` for appending, writes the buffered
  /// text to it, and switches to `Writing`. If the log is already writing to
  /// `path`, it does nothing; if it is writing to another file, that file is
  /// closed first. If writing the buffer fails, it prints a warning to
  /// stderr and switches the log off.
  pub fn open(&mut self, path: &str) -> Result<(),CluEError>{
    if let LogState::Writing(_, current) = self{
      if current == path{
        return Ok(());
      }
    }

    let file = OpenOptions::new().create(true).append(true).open(path)
      .map_err(|err| CluEError::CannotOpenFile(path.to_string(),
            err.to_string()))?;

    let buffer = match std::mem::replace(self, LogState::Off){
      LogState::Buffering(buffer) => buffer,
      LogState::Writing(..) | LogState::Off => String::new(),
    };

    *self = LogState::Writing(file, path.to_string());
    self.record(&buffer);
    Ok(())
  }
  //----------------------------------------------------------------------------
  /// This function drops the buffered text and switches the log off.
  pub fn discard(&mut self){
    *self = LogState::Off;
  }
  //----------------------------------------------------------------------------
  /// This function closes the log file, if one is open, and switches to an
  /// empty buffer.
  pub fn finish(&mut self){
    *self = LogState::Buffering(String::new());
  }
  //----------------------------------------------------------------------------
  // This function prints the warning for the failed write `err` to stderr
  // and switches the log off.
  fn stop(&mut self, err: &std::io::Error){
    if let LogState::Writing(_, path) = self{
      eprintln!("Warning: cannot write to \"{}\": {}; the output log is not \
continued.", path, err);
    }
    *self = LogState::Off;
  }
}

/// The output log of the process.
static LOG: Mutex<LogState> = Mutex::new(LogState::Buffering(String::new()));

//------------------------------------------------------------------------------
// This function locks `LOG`. A panic of another thread while it held the
// lock leaves a valid state, so a poisoned lock is used as it is.
fn lock() -> MutexGuard<'static, LogState>{
  LOG.lock().unwrap_or_else(|poisoned| poisoned.into_inner())
}
//------------------------------------------------------------------------------
/// This function prints `text` to stderr if `to_stderr`, and to stdout
/// otherwise, and records it in the output log. The lock is held for both,
/// so the texts of different threads are not mixed.
/// It is called by the macros `clue_print!`, `clue_println!`, and
/// `clue_eprintln!`.
pub fn print_and_record(text: &str, to_stderr: bool){
  let mut log = lock();
  if to_stderr{
    eprint!("{}", text);
  }else{
    print!("{}", text);
  }
  log.record(text);
}
//------------------------------------------------------------------------------
/// This function opens the output log at `path` (`LogState::open`).
pub fn open(path: &str) -> Result<(),CluEError>{
  lock().open(path)
}
//------------------------------------------------------------------------------
/// This function drops the buffered text and switches the output log off
/// (`LogState::discard`).
pub fn discard(){
  lock().discard();
}
//------------------------------------------------------------------------------
/// This function closes the output log and starts an empty buffer
/// (`LogState::finish`).
pub fn finish(){
  lock().finish();
}
//------------------------------------------------------------------------------
/// This macro prints to stdout as `print!` does, and records the text in the
/// output log.
#[macro_export]
macro_rules! clue_print{
  ($($arg:tt)*) => {
    $crate::output_log::print_and_record(&format!($($arg)*), false)
  };
}
//------------------------------------------------------------------------------
/// This macro prints to stdout as `println!` does, and records the text in
/// the output log.
#[macro_export]
macro_rules! clue_println{
  () => {
    $crate::output_log::print_and_record("\n", false)
  };
  ($($arg:tt)*) => {
    $crate::output_log::print_and_record(
        &format!("{}\n", format_args!($($arg)*)), false)
  };
}
//------------------------------------------------------------------------------
/// This macro prints to stderr as `eprintln!` does, and records the text in
/// the output log.
#[macro_export]
macro_rules! clue_eprintln{
  () => {
    $crate::output_log::print_and_record("\n", true)
  };
  ($($arg:tt)*) => {
    $crate::output_log::print_and_record(
        &format!("{}\n", format_args!($($arg)*)), true)
  };
}
//------------------------------------------------------------------------------


#[cfg(test)]
mod tests{
  use super::*;

  //----------------------------------------------------------------------------
  // This function returns a path for test `name` under the system temp
  // directory; nothing exists at the path.
  fn temp_file(name: &str) -> String{
    let path = std::env::temp_dir().join(format!(
          "clue_oxide_test_output_log_{}_{}.txt", name, std::process::id()));
    let _ = std::fs::remove_file(&path);
    path.to_str().unwrap().to_string()
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_buffered_text_is_written_first_on_open(){
    let path = temp_file("open");
    let mut log = LogState::Buffering(String::new());
    log.record("title\n");
    log.record("note\n");
    log.open(&path).unwrap();
    assert!(matches!(&log, LogState::Writing(_, p) if *p == path));
    log.record("progress\n");
    assert_eq!(std::fs::read_to_string(&path).unwrap(),
        "title\nnote\nprogress\n");

    // Opening the same path again does nothing.
    log.open(&path).unwrap();
    log.record("end\n");
    assert_eq!(std::fs::read_to_string(&path).unwrap(),
        "title\nnote\nprogress\nend\n");
    std::fs::remove_file(&path).unwrap();
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_log_is_appended(){
    let path = temp_file("append");
    std::fs::write(&path, "earlier run\n").unwrap();
    let mut log = LogState::Buffering("title\n".to_string());
    log.open(&path).unwrap();
    log.record("done\n");
    assert_eq!(std::fs::read_to_string(&path).unwrap(),
        "earlier run\ntitle\ndone\n");
    std::fs::remove_file(&path).unwrap();
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_discard_drops_the_buffer(){
    let mut log = LogState::Buffering("title\n".to_string());
    log.discard();
    assert!(matches!(log, LogState::Off));
    log.record("progress\n");
    assert!(matches!(log, LogState::Off));
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_finish_returns_to_an_empty_buffer(){
    let path = temp_file("finish");
    let mut log = LogState::Buffering("first run\n".to_string());
    log.open(&path).unwrap();
    log.finish();
    assert!(matches!(&log, LogState::Buffering(b) if b.is_empty()));
    log.record("second run\n");
    assert!(matches!(&log, LogState::Buffering(b) if b == "second run\n"));
    assert_eq!(std::fs::read_to_string(&path).unwrap(), "first run\n");

    // After discard(), finish() also starts an empty buffer.
    log.discard();
    log.finish();
    assert!(matches!(&log, LogState::Buffering(b) if b.is_empty()));
    std::fs::remove_file(&path).unwrap();
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_open_switches_to_another_file(){
    let path_a = temp_file("switch_a");
    let path_b = temp_file("switch_b");
    let mut log = LogState::Buffering("a\n".to_string());
    log.open(&path_a).unwrap();
    log.open(&path_b).unwrap();
    log.record("b\n");
    assert_eq!(std::fs::read_to_string(&path_a).unwrap(), "a\n");
    assert_eq!(std::fs::read_to_string(&path_b).unwrap(), "b\n");
    std::fs::remove_file(&path_a).unwrap();
    std::fs::remove_file(&path_b).unwrap();
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_open_error(){
    let path = format!("{}/no_such_dir/output_log.txt",
        temp_file("missing_dir"));
    let mut log = LogState::Buffering("title\n".to_string());
    let err = log.open(&path).unwrap_err();
    assert!(matches!(&err, CluEError::CannotOpenFile(p,_) if *p == path));
    // The buffer is kept.
    assert!(matches!(&log, LogState::Buffering(b) if b == "title\n"));
  }
  //----------------------------------------------------------------------------
  // /dev/full accepts the open, and every write fails with "no space left
  // on device".
  #[cfg(target_os = "linux")]
  #[test]
  fn test_failing_write_switches_off(){
    let mut log = LogState::Buffering(String::new());
    log.open("/dev/full").unwrap();
    assert!(matches!(log, LogState::Writing(..)));
    log.record("progress\n");
    assert!(matches!(log, LogState::Off));

    // Writing the buffer on open fails the same way.
    let mut log = LogState::Buffering("title\n".to_string());
    log.open("/dev/full").unwrap();
    assert!(matches!(log, LogState::Off));
  }
  //----------------------------------------------------------------------------
  // Every form of the macros expands and runs. The global log is shared
  // with the other tests, which may run at the same time, so its contents
  // are not checked here; tests/test_output_log.rs checks them in a
  // separate process.
  #[test]
  fn test_macro_forms(){
    crate::clue_print!("");
    crate::clue_println!();
    crate::clue_println!("{} {}", "output_log", "test");
    crate::clue_eprintln!("{}", "output_log test");
  }
}
