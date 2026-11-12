// Build script of CluE Oxide.
// It records the version of the Rust compiler in the environment variable
// CLUE_RUSTC_VERSION, which the code reads with env!("CLUE_RUSTC_VERSION")
// and writes to the record of an output directory (clue_manifest.toml).
// Different compilers can format a configuration differently, so the
// version helps to explain a configuration mismatch.
// If the compiler cannot be run, the version is "unknown"; the build never
// fails because of this script.
fn main(){
  // Cargo sets RUSTC for build scripts.
  let rustc = std::env::var("RUSTC").unwrap_or_else(|_| "rustc".to_string());
  let version = std::process::Command::new(rustc).arg("-V").output().ok()
    .and_then(|output| String::from_utf8(output.stdout).ok())
    .map(|s| s.trim().to_string())
    .filter(|s| !s.is_empty())
    .unwrap_or_else(|| "unknown".to_string());
  println!("cargo:rustc-env=CLUE_RUSTC_VERSION={}", version);
  println!("cargo:rerun-if-changed=build.rs");
}
