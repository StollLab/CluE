use num_complex::Complex;
use crate::CluEError;

// Numbers
pub const I: Complex<f64> = Complex::<f64>{re: 0.0, im: 1.0};
pub const ONE: Complex<f64> = Complex::<f64>{re: 1.0, im: 0.0};
pub const ZERO: Complex<f64> = Complex::<f64>{re: 0.0, im: 0.0};
pub const PI: f64 = std::f64::consts::PI;
pub const SQRT2: f64 = std::f64::consts::SQRT_2;
pub const SQRT2_INV: f64 = std::f64::consts::FRAC_1_SQRT_2;

// Units
pub const METER: f64 = 1.0;
pub const KILOMETER: f64 = METER*1e3;
pub const CENTIMETER: f64 = METER*1e-2;
pub const MILLIMETER: f64 = METER*1e-3;
pub const MICROMETER: f64 = METER*1e-6;
pub const MICRON: f64 = MICROMETER;
pub const NANOMETER: f64 = METER*1e-9;
pub const ANGSTROM: f64 = METER*1e-10;
pub const PICOMETER: f64 = METER*1e-12;
pub const FEMTOMETER: f64 = METER*1e-15;

pub const SECOND: f64 = 1.0;
pub const MILLISECOND: f64 = SECOND*1e-3;
pub const MICROSECOND: f64 = SECOND*1e-6;
pub const NANOSECOND: f64 = SECOND*1e-9;
pub const PICOSECOND: f64 = SECOND*1e-12;
pub const FEMTOSECOND: f64 = SECOND*1e-15;

pub const TESLA: f64 = 1.0;
pub const MILLITESLA: f64 = TESLA*1e-3;
pub const GAUSS: f64 = TESLA*1e-4;
pub const MILLIGAUSS: f64 = GAUSS*1e-3;
pub const KILOGAUSS: f64 = GAUSS*1e3;

pub const JOULE: f64 = 1.0;
pub const KILOJOULE: f64 = JOULE*1e3;
pub const CALORIE: f64 = 4.184*JOULE;
pub const KILOCALORIE: f64 = CALORIE*1e3;

pub const KELVIN: f64 = 1.0;

pub const HERTZ: f64 = 1.0/SECOND;
pub const KILOHERTZ: f64 = HERTZ*1e3;
pub const MEGAHERTZ: f64 = HERTZ*1e6;
pub const GIGAHERTZ: f64 = HERTZ*1e9;
pub const TERAHERTZ: f64 = HERTZ*1e12;


pub const KILOGRAM: f64 = JOULE/METER/METER*SECOND*SECOND;
pub const AMPERE: f64 = KILOGRAM/SECOND/SECOND/TESLA;
pub const NEWTON: f64 = KILOGRAM*METER/SECOND/SECOND;
pub const COULOMB: f64 = AMPERE*SECOND; 
pub const VOLT: f64 = JOULE/COULOMB;

// 2022 https://physics.nist.gov/cgi-bin/cuu/Value?h
pub const PLANCK: f64 = 6.62607015e-34*JOULE/HERTZ;

// 2021: https://physics.nist.gov/cgi-bin/cuu/Value?hbar 
//pub const HBAR: f64 = 1.054571817e-34 * JOULE * SECOND;
pub const HBAR: f64 = (PLANCK/2.0/PI) * JOULE * SECOND;

// 2022 https://www.physics.nist.gov/cgi-bin/cuu/Value?mub
pub const MUB:f64 = 9.2740100783e-24 * JOULE/TESLA;

// 2022 https://www.physics.nist.gov/cgi-bin/cuu/Value?mun
pub const MUN: f64 = 5.0507837461e-27 * JOULE/TESLA;

// 2022 https://www.physics.nist.gov/cgi-bin/cuu/Value?mu0
pub const MU0: f64 = 1.25663706212e-6 * NEWTON/AMPERE/AMPERE;

// 2022 https://physics.nist.gov/cgi-bin/cuu/Value?k
pub const BOLTZMANN: f64 = 1.380649e-23 * JOULE/KELVIN;

// 2022 https://physics.nist.gov/cgi-bin/cuu/Value?na
pub const AVOGADRO: f64 = 6.02214076e23;

pub const JOULE_PER_MOLE: f64 = JOULE/AVOGADRO;
pub const KILOJOULE_PER_MOLE: f64 = KILOJOULE/AVOGADRO;
pub const CALORIE_PER_MOLE: f64 = CALORIE/AVOGADRO;
pub const KILOCALORIE_PER_MOLE: f64 = KILOCALORIE/AVOGADRO;

// 2022 https://www.physics.nist.gov/cgi-bin/cuu/Value?gem 
pub const ELECTRON_G: f64 = 2.00231930436256;

// 2022 https://physics.nist.gov/cgi-bin/cuu/Value?e
pub const ELEMENTARY_CHARGE: f64 = 1.602176634e-19 * COULOMB;

pub const ELECTRON_VOLT: f64 = ELEMENTARY_CHARGE * VOLT;
pub const MILLIELECTRON_VOLT: f64 = ELECTRON_VOLT*1e-3;

// 2025 https://physics.nist.gov/cgi-bin/cuu/Value?hr
pub const HARTREE: f64 =  4.3597447222060e-18 * JOULE;

pub const BARN: f64 = 100.0 * FEMTOMETER * FEMTOMETER;

// 2025 https://physics.nist.gov/cgi-bin/cuu/Value?bohrrada0
pub const BOHR_RADIUS: f64 = 5.29177210544e-11 * METER;

// Conversions
pub const JOULES_TO_HERTZ: f64 = 1.0/PLANCK * HERTZ/JOULE;
pub const RAD_PER_S_TO_HZ: f64 = 0.5/PI;
pub const HZ_TO_RAD_PER_S: f64 = 2.0*PI;

pub const C3_TUNNEL_SPLITTING_TO_EXCHANGE_COUPLING: f64 = -2.0/3.0;


pub fn energy_unit_to_hertz(unit: &str) -> Result<f64,CluEError>
{
  match unit{
    "Hz" => Ok(HERTZ),
    "kHz" => Ok(KILOHERTZ),
    "MHz" => Ok(MEGAHERTZ),
    "GHz" => Ok(GIGAHERTZ),
    "THz" => Ok(TERAHERTZ),
    "meV" => Ok(ELECTRON_VOLT*JOULES_TO_HERTZ),
    "eV" => Ok(MILLIELECTRON_VOLT*JOULES_TO_HERTZ),
    "Eh" => Ok(HARTREE*JOULES_TO_HERTZ),
    "J" => Ok(JOULE*JOULES_TO_HERTZ),
    "kJ" => Ok(KILOJOULE*JOULES_TO_HERTZ),
    "kJ/mol" => Ok(KILOJOULE_PER_MOLE*JOULES_TO_HERTZ),
    "cal/mol" => Ok(CALORIE_PER_MOLE*JOULES_TO_HERTZ),
    "cal" => Ok(CALORIE*JOULES_TO_HERTZ),
    "kcal" => Ok(KILOCALORIE*JOULES_TO_HERTZ),
    _ => Err(CluEError::UnrecognizedUnit(unit.to_string())),
  }
}

pub fn distance_unit_to_meters(unit: &str) -> Result<f64,CluEError>
{
  match unit{
    "fm" => Ok(FEMTOMETER),
    "pm" => Ok(PICOMETER),
    "nm" => Ok(NANOMETER),
    "μm" => Ok(MICROMETER),
    "mm" => Ok(MILLIMETER),
    "cm" => Ok(CENTIMETER),
    "m" => Ok(METER),
    "km" => Ok(KILOMETER),
    "Å" => Ok(ANGSTROM),
    _ => Err(CluEError::UnrecognizedUnit(unit.to_string())),
  }
}

pub fn magnetic_field_unit_to_tesla(unit: &str) -> Result<f64,CluEError>
{
  match unit{
    "mT" => Ok(MILLITESLA),
    "T" => Ok(TESLA),
    "mG" => Ok(MILLIGAUSS),
    "G" => Ok(GAUSS),
    "kG" => Ok(KILOGAUSS),
    _ => Err(CluEError::UnrecognizedUnit(unit.to_string())),
  }
}

pub fn time_unit_to_seconds(unit: &str) -> Result<f64,CluEError>
{
  match unit{
    "ps" => Ok(PICOSECOND),
    "fs" => Ok(FEMTOSECOND),
    "ns" => Ok(NANOSECOND),
    "μs" => Ok(MICROSECOND),
    "ms" => Ok(MILLISECOND),
    "s" => Ok(SECOND),
    _ => Err(CluEError::UnrecognizedUnit(unit.to_string())),
  }
}

