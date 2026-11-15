use crate::clue_println;

/// This function prints the version.
pub fn print_version(){
  // Remember, when updating the version, also update it in Cargo.toml.
  clue_println!("CluE 1.0.0-dev");
}
