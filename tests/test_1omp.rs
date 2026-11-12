use clue_oxide::config::Config;
use clue_oxide::run;
use clue_oxide::signal::Signal;

const ERROR_THRESHOLD: f64 = 1e-12;
#[test]
fn int_test_1omp(){
 let mut config = Config::from_toml_file(
     "assets/1omp_K26R1_0.003228966616048703.toml").unwrap();

 // A fresh output directory, so that the run computes the signal.
 let out_dir = std::env::temp_dir().join(format!(
       "clue_oxide_int_test_1omp_{}", std::process::id()));
 let _ = std::fs::remove_dir_all(&out_dir);
 config.output_directory = Some(out_dir.to_str().unwrap().to_string());
 
 let (_time_axis,signal) = run(config).unwrap();
 let _ = std::fs::remove_dir_all(&out_dir);

 let signal = Signal{data: signal};

 let ref_signal = Signal::read_from_csv(
     "assets/1omp_K26R1_0.003228966616048703_signal.csv").unwrap();

 assert_eq!(signal.len(),ref_signal.len());

 let c = 1.0/(signal.len() as f64).sqrt(); 
 let rmsd = c*(&signal - &ref_signal).norm();

 println!("DB: RMSD = {}.",rmsd);
 assert!(rmsd < ERROR_THRESHOLD);
}
