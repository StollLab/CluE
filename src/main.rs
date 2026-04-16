///    CluE (Cluster Evolution) is a spin dynamics simulation program for
///    electron spin decoherence.
///    Copyright (C) 2022  Samuel M. Jahn
///
///    This program is free software: you can redistribute it and/or modify
///    it under the terms of the GNU General Public License as published by
///    the Free Software Foundation, version 3 of the License.
///
///    This program is distributed in the hope that it will be useful,
///    but WITHOUT ANY WARRANTY; without even the implied warranty of
///    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
///    GNU General Public License for more details.
///
///    You should have received a copy of the GNU General Public License
///    along with this program.  If not, see <https://www.gnu.org/licenses/>.
use clue_oxide::config::command_line_input::CommandLineInput;
use clue_oxide::config::Config;
use clue_oxide::info;

use std::time::Instant;
fn main() {

  // Record start time
  let now = Instant::now();

  // Collect input arguments.
  let args: Vec<String> = std::env::args().collect();

  // Parse comand line input.
  let input = CommandLineInput::new(args).unwrap_or_else(|err|{
    eprintln!("CluE Error: {}.",err);
    eprintln!("Try \"clue --help\" for more information.");
    std::process::exit(1);
  });

  // Decide what information to display.
  if input.show_help{
    info::help::print_help();
    return;
  }

  if input.show_license{
    info::license::print_license();
    return;
  }

  if input.show_warranty{
    info::warranty::print_warranty();
    return;
  }

  if input.show_version{
    info::version::print_version();
    return;
  }

  if input.show_title{
    info::title::print_title();
  }


  // Parse input files.
  let mut config = Config::read_input(input).unwrap_or_else(
      |err| {
        eprintln!("CluE Error: {}.",err);
        std::process::exit(1);
      });

  config.set_defaults().unwrap_or_else(
      |err| {
        eprintln!("CluE Error: {}.",err);
        std::process::exit(1);
      });

  /*
  config.set_time_axis().unwrap_or_else(
      |err| {
        eprintln!("CluE Error: {}.",err);
        std::process::exit(1);
      });
  */
  // Run simulations.
  clue_oxide::run(config).unwrap_or_else(
      |err| {
        eprintln!("CluE Error: {}.",err);
        std::process::exit(1);
      });

  // Exit.
  let elapsed_time = now.elapsed();
  eprintln!("CluE successfully completed in {} s.", elapsed_time.as_secs());
}




