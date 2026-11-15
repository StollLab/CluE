// Integration test of the output log, <output>/output_log.txt
// (plan_multistructures.md, 1.G item 12). It runs the clue_oxide binary as
// a subprocess, so that no other test shares the process-wide log.
use std::path::{Path,PathBuf};
use std::process::{Command,Output};

const ASSET: &str = "assets/1omp_K26R1_0.003228966616048703.toml";
const STRUCTURE: &str = "assets/1omp_K26R1_0.003228966616048703.pdb";
const LOG_FILE: &str = "output_log.txt";

//------------------------------------------------------------------------------
// This function returns an empty directory for test `name` under the system
// temp directory.
fn temp_dir(name: &str) -> PathBuf{
  let dir = std::env::temp_dir().join(format!(
        "clue_oxide_int_test_output_log_{}_{}", name, std::process::id()));
  let _ = std::fs::remove_dir_all(&dir);
  std::fs::create_dir_all(&dir).unwrap();
  dir
}
//------------------------------------------------------------------------------
// This function writes, in `work`, the 1omp test input with a radius of 8 Å
// and 11 time points, the absolute structure path, `output.directory =
// output_dir`, and `rng_seed = seed`, and returns the path of the input file.
fn write_config(work: &Path, output_dir: &str, seed: u64) -> PathBuf{
  let root = Path::new(env!("CARGO_MANIFEST_DIR"));
  let structure = root.join(STRUCTURE);
  let toml = std::fs::read_to_string(root.join(ASSET)).unwrap();
  for old in [STRUCTURE, "radius = 300", "number_timepoints = [301]",
    "[output]"]
  {
    assert!(toml.contains(old), "{}", old);
  }
  let toml = toml
    .replace(STRUCTURE, structure.to_str().unwrap())
    .replace("radius = 300", "radius = 8")
    .replace("number_timepoints = [301]", "number_timepoints = [11]")
    .replace("[output]", &format!("[output]\n  directory = '{}'",
          output_dir));
  let toml = format!("{}\n[run]\nnumber_runs = 1\nrng_seed = {}\n",
      toml, seed);
  let path = work.join(format!("input_{}.toml", seed));
  std::fs::write(&path, toml).unwrap();
  path
}
//------------------------------------------------------------------------------
// This function runs the binary on `config` with the working directory
// `work`.
fn run_clue(config: &Path, work: &Path) -> Output{
  Command::new(env!("CARGO_BIN_EXE_clue_oxide"))
    .arg(config)
    .current_dir(work)
    .output()
    .unwrap()
}
//------------------------------------------------------------------------------
fn lines(bytes: &[u8]) -> Vec::<String>{
  String::from_utf8(bytes.to_vec()).unwrap()
    .lines().map(|l| l.to_string()).collect()
}
//------------------------------------------------------------------------------
// This function checks that the log `section` consists of exactly the
// lines of the captured stdout and stderr of `output`, with the stdout lines
// in their stdout order and the stderr lines in their stderr order. (The two
// pipes are read separately, so the order between them is not known.)
fn assert_section_is_terminal_output(section: &str, output: &Output){
  let out = lines(&output.stdout);
  let err = lines(&output.stderr);
  let (mut i, mut j) = (0, 0);
  for line in section.lines(){
    if i < out.len() && line == out[i]{
      i += 1;
    }else if j < err.len() && line == err[j]{
      j += 1;
    }else{
      panic!("log line {:?} is not the next stdout line {:?} or the next \
stderr line {:?}", line, out.get(i), err.get(j));
    }
  }
  assert_eq!((i, j), (out.len(), err.len()),
      "not every printed line is in the log");
}
//------------------------------------------------------------------------------
// This function checks that a log section starts with the title, and ends
// with the success message of main.
fn assert_section_frame(section: &str){
  let first = section.lines().next().unwrap();
  assert!(first.len() == 78 && first.chars().all(|c| c == '#'), "{}", first);
  assert!(section.contains("\nCLuE Copyright"));
  let last = section.lines().last().unwrap();
  assert!(last.starts_with("CluE successfully completed in ")
      && last.ends_with(" s."), "{}", last);
}
//------------------------------------------------------------------------------
#[test]
fn int_test_output_log(){
  let work = temp_dir("work");
  let out = work.join("out");
  let out_str = out.to_str().unwrap();
  let log_path = out.join(LOG_FILE);

  // A new output: the log holds exactly what was printed.
  let config = write_config(&work, out_str, 5);
  let first = run_clue(&config, &work);
  assert!(first.status.success(), "{}",
      String::from_utf8_lossy(&first.stderr));
  let log_1 = std::fs::read_to_string(&log_path).unwrap();
  assert_section_is_terminal_output(&log_1, &first);
  assert_section_frame(&log_1);
  assert!(log_1.contains("System instance: 1/1."));
  assert!(log_1.contains("Found "));

  // A second run into the complete output appends its own section; the
  // first section is unchanged.
  let second = run_clue(&config, &work);
  assert!(second.status.success());
  let log_2 = std::fs::read_to_string(&log_path).unwrap();
  let Some(section_2) = log_2.strip_prefix(&log_1) else{
    panic!("the first section changed");
  };
  assert_section_is_terminal_output(section_2, &second);
  assert_section_frame(section_2);
  assert!(section_2.contains("already complete"));

  // A run that stops with OutputSeedMismatch does not touch the log.
  let other_seed = write_config(&work, out_str, 6);
  let mismatch = run_clue(&other_seed, &work);
  assert!(!mismatch.status.success());
  assert!(String::from_utf8_lossy(&mismatch.stderr)
      .contains(&format!("CluE Error: the output in \"{}\" was computed with \
run.rng_seed = 5, but run.rng_seed is 6", out_str)), "{}",
      String::from_utf8_lossy(&mismatch.stderr));
  assert_eq!(std::fs::read_to_string(&log_path).unwrap(), log_2);

  std::fs::remove_dir_all(&work).unwrap();
}
//------------------------------------------------------------------------------
#[test]
fn int_test_no_output_log_without_output(){
  // The working directory holds only the input file; with output.directory
  // = "", the run writes nothing, so no log either.
  let work = temp_dir("not_saved");
  let config = write_config(&work, "", 5);
  let output = run_clue(&config, &work);
  assert!(output.status.success(), "{}",
      String::from_utf8_lossy(&output.stderr));
  let entries = std::fs::read_dir(&work).unwrap()
    .map(|e| e.unwrap().file_name().to_string_lossy().to_string())
    .collect::<Vec::<String>>();
  assert_eq!(entries, vec!["input_5.toml"]);
  std::fs::remove_dir_all(&work).unwrap();
}
//------------------------------------------------------------------------------
