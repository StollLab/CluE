use clue_oxide::config::Config;
use clue_oxide::run;
use clue_oxide::signal::Signal;

const ERROR_THRESHOLD: f64 = 1e-12;
#[test]
fn test_1omp(){
 let config = Config::from_toml_file(
     "assets/1omp_K26R1_0.003228966616048703.toml").unwrap();
 
 let (_time_axis,signal) = run(config).unwrap();

 let signal = Signal{data: signal};

 let ref_signal = Signal::read_from_csv(
     "assets/1omp_K26R1_0.003228966616048703_signal.csv").unwrap();

 assert_eq!(signal.len(),ref_signal.len());

 let c = 1.0/(signal.len() as f64).sqrt(); 
 let rmsd = c*(&signal - &ref_signal).norm();

 assert!(rmsd < ERROR_THRESHOLD);
}
