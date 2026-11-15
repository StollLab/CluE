//! CluE Oxide (Cluster Evolution Oxide) is a spin dynamics simulation program 
//! for electron spin decoherence.
//!
//! CluE Oxide implements Yang an Liu's cluster correlation expansion for an
//! input structure, usually a PDB file.  
//! 
//! Documentation on CluE Oxide the program is provided in the 
//![manual](https://github.com/jahnsam/clue_oxide/blob/main/manual/CluE_Oxide.pdf).
pub mod config;
pub mod clue_errors;
pub mod cluster;
pub mod cluster_methods;
pub mod elements;
pub mod info;
pub mod isotopes;
pub mod integration_grid;
pub mod physical_constants;
pub mod quantum;
pub mod structure;
pub mod signal;
pub mod space_3d;
pub mod symmetric_list_2d;
pub mod math;
pub mod misc;
pub mod kmeans;
pub mod io;
pub mod output_manifest;
pub mod output_log;



use crate::config::Config;
use crate::config::command_line_input::CommandLineInput;
use crate::clue_errors::CluEError;
use crate::structure::Structure;
use crate::quantum::tensors::HamiltonianTensors;
use crate::cluster::build_adjacency_list::build_adjacency_list;
use crate::cluster::find_clusters::find_clusters;
use crate::output_manifest::{
  OutputState, OutputStatus, ResumeSummary,
  OUTPUT_LOG_FILE, RNG_SEED_FILE, SIGNAL_FILE,
  draw_main_seed, inspect_output,
};
use crate::signal::calculate_signal;
use crate::signal::{TauAxes,load_csv_to_vec_signals,write_vec_signals};

use num_complex::Complex;

/// This function follows the config instruction to central simulate spin 
/// decoherence.
/// It first checks the config with `Config::check_completeness()`, and the
/// output directory with `output_manifest::inspect_output()`: a run into a
/// directory with a complete output loads the saved signal, a run into a 
/// directory with a compatible incomplete output resumes it, and a run into
/// a directory with the output of a different simulation stops with an
/// error.
/// Without `run.rng_seed`, the main seed is drawn from the operating 
/// system's entropy; it is printed and saved in `rng_seed.txt`.
/// If the output is saved, everything printed during the run is also 
/// appended to `<output>/output_log.txt` (`output_log`); the caller ends 
/// the log with `output_log::finish()`.
/// It returns the time columns of the run (`TauAxes`, in units.time) and
/// the highest-order signal: `tau_axes.tau[i]` and `tau_axes.tau2[i]` are
/// the delays of `signal[i]`, as in the columns of the signal files.
pub fn run(mut config: Config) 
  -> Result<( TauAxes, Vec::<Complex::<f64>> ),CluEError>
{

  // Report missing settings and files before any work is done.
  config.check_completeness()?;

  // Read and hash the structure files once for the whole run.
  config.resolve_structures()?;

  config.set_tau_axis()?;

  // Decide what to do with the output directory before any work is done.
  let mut status = inspect_output(&config)?;

  // Draw the main seed, if it is not set, and name the default output 
  // directory after it.
  if config.rng_seed.is_none(){
    config.rng_seed = Some(draw_main_seed());
    status = inspect_output(&config)?;
  }
  let Some(main_seed) = config.rng_seed else{
    return Err(CluEError::NoRngSeed);
  };

  // The time columns of the signal files, in units.time.
  let tau_axes = TauAxes::from_config(&config)?;

  // The output log: what was printed so far is buffered, and is written 
  // to <output>/output_log.txt once the directory exists (output_log.rs).
  // A run that stops before this point writes nothing.
  if let OutputStatus::Complete(dir) = &status{
    output_log::open(&format!("{}/{}", dir, OUTPUT_LOG_FILE))?;
    clue_println!("Random number seed: {}.", main_seed);
    if let Some(signal) = load_complete_signal(dir, &config, &tau_axes){
      clue_println!("Output in {} is already complete; no calculations are \
needed.", dir);
      return Ok((tau_axes, signal));
    }
    clue_println!("The signal in {} cannot be read; it is averaged again from \
the saved instances.", dir);
  }

  let mut output = OutputState::start(&config, &status)?;

  match output.directory(){
    Some(dir) => {
      // Does nothing if the Complete branch above opened the same log.
      output_log::open(&format!("{}/{}", dir, OUTPUT_LOG_FILE))?;
      clue_println!("Random number seed: {} (saved in {}/{}).", 
          main_seed, dir, RNG_SEED_FILE);
    },
    None => {
      output_log::discard();
      clue_println!("Random number seed: {}.", main_seed);
    },
  }

  let order_n_signals 
    = calculate_signal::calculate_signals(&config, &mut output, &tau_axes)?;

  let max_size = order_n_signals.len();

  if let Some(save_dir) = output.directory(){
    let save_path = format!("{}/{}", save_dir, SIGNAL_FILE);
    let headers = (0..max_size)
      .map(|ii| format!("signal_{}",ii))
      .collect::<Vec::<String>>();
    write_vec_signals(&tau_axes, &order_n_signals, headers, &save_path)?;
  }

  let signal = order_n_signals[max_size-1].data.clone();
  Ok((tau_axes, signal))
}
//------------------------------------------------------------------------------
// This function loads the highest-order signal from the signal.csv of the
// complete output in `dir`, or returns `None` if the file does not have the
// signals of orders 0, ..., max_cluster_size with the run's number of time
// points, or does not have the time columns `tau_axes` of the run.
fn load_complete_signal(dir: &str, config: &Config, tau_axes: &TauAxes) 
  -> Option<Vec::<Complex::<f64>>>
{
  let max_cluster_size = config.max_cluster_size?;
  let n_tot = config.get_total_number_timesteps();
  let (read_axes, signals) 
    = load_csv_to_vec_signals(&format!("{}/{}", dir, SIGNAL_FILE)).ok()?;
  if signals.len() != max_cluster_size + 1 
    || signals.iter().any(|s| s.len() != n_tot)
    || read_axes.as_ref() != Some(tau_axes)
  {
    return None;
  }
  signals.last().map(|s| s.data.clone())
}
//------------------------------------------------------------------------------
/// `CheckReport` is the result of `check()`: the result of checking the 
/// input file, and, if the input is valid, what a run would do with its
/// output directory (`output_manifest::inspect_output()`).
#[derive(Debug,Clone,PartialEq)]
pub struct CheckReport{
  pub input: Result<(),CluEError>,
  pub output: Option<Result<OutputStatus,CluEError>>,
}

impl CheckReport{
  //----------------------------------------------------------------------------
  /// This function returns `true` if the input is valid and a run would not
  /// stop because of its output directory.
  pub fn is_ok(&self) -> bool{
    self.input.is_ok() && matches!(self.output, Some(Ok(_)))
  }
  //----------------------------------------------------------------------------
  /// This function returns the report as lines, each with `true` if it 
  /// reports an error: one line for the input, then one for the output.
  pub fn messages(&self) -> Vec::<(bool,String)>{
    let mut lines = Vec::<(bool,String)>::new();
    match &self.input{
      Ok(()) => lines.push((false, "input OK".to_string())),
      Err(err) => lines.push((true, format!("input: CluE Error: {}.", err))),
    }
    let line = match &self.output{
      None => (true, 
        "output: not checked, because the input has errors".to_string()),
      Some(Err(err)) => (true, format!("output: CluE Error: {}.", err)),
      Some(Ok(OutputStatus::NotSaved)) => (false, 
        "output not saved (output.directory = \"\"); nothing to check"
        .to_string()),
      Some(Ok(OutputStatus::NewWithoutSeed)) => (false, 
        "output: run.rng_seed is not set, so a run writes to a new directory \
CluE-<hash>; nothing to check".to_string()),
      Some(Ok(OutputStatus::New(dir))) => (false, format!("output {}: does \
not exist or is empty; a run starts a new simulation", dir)),
      Some(Ok(OutputStatus::Complete(dir))) => (false, format!("output {}: \
complete; a run will do no calculations", dir)),
      Some(Ok(OutputStatus::Resume(dir, summaries))) => (false, format!(
            "output {}: a run resumes: {}; the orientation seeds are checked \
by the run", dir, resume_text(summaries))),
    };
    lines.push(line);
    lines
  }
  //----------------------------------------------------------------------------
}
//------------------------------------------------------------------------------
// This function describes how a run resumes the structures of `summaries`.
fn resume_text(summaries: &[ResumeSummary]) -> String{
  let describe = |s: &ResumeSummary| if s.new_structure{
    format!("{} instances to compute", s.number_runs)
  }else{
    format!("{} of {} instances complete ({} orientations complete in \
unfinished instances), {} instances added", s.instances_complete, 
        s.number_runs, s.orientations_complete, s.instances_added)
  };
  if summaries.len() == 1 && !summaries[0].new_structure{
    describe(&summaries[0])
  }else{
    summaries.iter()
      .map(|s| format!("structure \"{}\" ({}{}): {}", s.file, s.id,
            if s.new_structure {", new"} else {""}, describe(s)))
      .collect::<Vec::<String>>().join("; ")
  }
}
//------------------------------------------------------------------------------
/// This function checks the input file named in `input` without running 
/// it, and then the output directory the run would use.
/// The input check parses the file (including the key check), sets the
/// defaults, and runs `Config::check_completeness()`. If the input is 
/// valid, `output_manifest::inspect_output()` reports whether a run would
/// start a new simulation, find a complete output, resume an output, or
/// stop because the directory holds the output of a different simulation.
/// Nothing is written.
pub fn check(input: CommandLineInput) -> CheckReport{
  let config = Config::read_input(input).and_then(|mut config|{
    config.set_defaults()?;
    config.check_completeness()?;
    config.resolve_structures()?;
    config.set_tau_axis()?;
    Ok(config)
  });
  match config{
    Ok(config) => CheckReport{
      input: Ok(()),
      output: Some(inspect_output(&config)),
    },
    Err(err) => CheckReport{input: Err(err), output: None},
  }
}
