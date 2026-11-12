pub const DEFAULT_UNIT_ENERGY: &str = "MHz";
pub const DEFAULT_UNIT_DISTANCE: &str = "Å";
pub const DEFAULT_UNIT_MAGNETIC_FIELD: &str = "T";
pub const DEFAULT_UNIT_TIME: &str = "μs";

// General Keys
pub const KEY_CUTOFF_COUPLING: &str = "coupling_xx_yy";
pub const KEY_CUTOFF_DELTA_HF: &str = "delta_hyperfine_zz";
pub const KEY_CUTOFF_DIPOLE_PERP: &str = "point_dipole_perpendicular";
pub const KEY_CUTOFF_DISTANCE: &str = "distance";
pub const KEY_CUTOFF_HAHN_MOD_DEPTH: &str = "hahn_mod_depth";
pub const KEY_CUTOFF_HAHN_TAYLOR_4: &str = "hahn_taylor_4";
pub const KEY_CUTOFF_DELTA_ZEEMAN: &str = "delta_zeeman";

pub const KEY_OUT_AUX_SIGS: &str = "auxiliary_signals";
pub const KEY_OUT_BATH: &str = "bath";
pub const KEY_OUT_DET_SPIN: &str = "detected_spin";
pub const KEY_OUT_CLUSTERS: &str = "clusters";
pub const KEY_OUT_INFO: &str = "info";
pub const KEY_OUT_EXCHANGE_GROUPS: &str = "exchange_groups";
pub const KEY_OUT_METHYL_PARTITIONS: &str = "methyl_partitions";
pub const KEY_OUT_ORI_GRID: &str = "orientation_grid";
pub const KEY_OUT_ORI_SIGS: &str = "orientation_signals";
pub const KEY_OUT_PART_TAB: &str = "partition_table";
pub const KEY_OUT_SANS_SPIN_SIGS: &str = "sans_spin_signals";
pub const KEY_OUT_STRUC_PDB: &str = "structure_pdb";
pub const KEY_OUT_TENSORS: &str = "tensors";

pub const KEY_DENSITY_MATRIX_ID: &str  = "uniform";
pub const KEY_DENSITY_MATRIX_RANDOM: &str  = "random";
pub const KEY_DENSITY_MATRIX_THERMAL: &str  = "thermal";
pub const KEY_DENSITY_MATRIX_ZEEMAN: &str  = "zeeman";

pub const KEY_DET_FRAME_EIGEN: &str  = "eigen";
pub const KEY_DET_FRAME_ZEEMAN: &str  = "zeeman";

pub const KEY_PARTITION_EX_GROUPS: &str = "exchange_groups";
pub const KEY_PARTITION_PARTICLE: &str = "singles";
pub const KEY_PARTITION_KMEANS: &str = "kmeans";
pub const KEY_PARTITION_RESTRICTED_KMEANS: &str = "restricted_kmeans";

pub const KEY_ORI_LEBEDEV: &str = "lebedev";
pub const KEY_ORI_RANDOM: &str = "random";
pub const KEY_ORI_FILE: &str = "file";
pub const KEY_ORI_VECTOR: &str = "vector";
pub const KEY_ORI_CIRCLE: &str = "circle";
pub const KEY_ORI_VECTORGRID: &str = "vector_grid";


pub const KEY_EIG_VALUES: &str = "values";
pub const KEY_EIG_AXES: &str = "axes";
pub const KEY_EIG_X_AXIS: &str = "x";
pub const KEY_EIG_Y_AXIS: &str = "y";
pub const KEY_EIG_Z_AXIS: &str = "z";

pub const KEY_IST_COEF: &str = "spherical_coefficients";
pub const KEY_STEVENS_COEF: &str = "stevens_coefficients";

// Group Keys
pub const KEY_NAME: &str = "name";

pub const KEY_DROP_PROB: &str = "drop_probability";
pub const KEY_ABUNDANCE_TOL: &str = "abundance_tolerance";

pub const KEY_ISO_COSUBSTITUTE: &str = "cosubstitute";


// Selection Keys
// The group keys that select particles, and the keys of their tables.
pub const KEY_CELLS: &str = "cells";
pub const KEY_INCLUDE: &str = "include";
pub const KEY_EXCLUDE: &str = "exclude";
pub const KEY_BONDED_TO: &str = "bonded_to";

pub const KEY_SELE_DISTANCE: &str = "distance";

// The atom attributes, used in include, exclude, and their bonded_to tables.
pub const KEY_SELE_CHAIN_IDS: &str = "chain_ids";
pub const KEY_SELE_ELEMENTS: &str = "elements";
pub const KEY_SELE_INDICES: &str = "indices";
pub const KEY_SELE_NAMES: &str = "names";
pub const KEY_SELE_RESIDUES: &str = "residues";
pub const KEY_SELE_RES_NUMS: &str = "residue_numbers";
pub const KEY_SELE_SERIALS: &str = "serials";

// Isotope Key
pub const KEY_ISO_ABUNDACE: &str = "abundance";
pub const KEY_ISO_ACTIVE: &str = "active";
pub const KEY_ISO_G_MATRIX: &str = "g_matrix";
pub const KEY_ISO_HYPERFINE: &str = "hyperfine";
pub const KEY_ISO_ELEC_QUADRUPOLE: &str = "electric_quadrupole";
pub const KEY_ISO_ZEROFIELD: &str = "zerofield";
pub const KEY_ISO_EXCHANGE_COUPLING: &str = "exchange_coupling";
pub const KEY_ISO_C3_TUNNEL_SPLITTING: &str = "c3_tunnel_splitting";
pub const KEY_ISO_: &str = "";

pub const KEY_VEC_SPECIFIER_FROM: &str = "from";
pub const KEY_VEC_SPECIFIER_FROM_BONDED_TO: &str = "from_bonded_to";
pub const KEY_VEC_SPECIFIER_FROM_SAME_MOLECULE_AS: &str = "from_same_molecule_as";
pub const KEY_VEC_SPECIFIER_TO: &str = "to";
pub const KEY_VEC_SPECIFIER_TO_BONDED_TO: &str = "to_bonded_to";
pub const KEY_VEC_SPECIFIER_TO_SAME_MOLECULE_AS: &str = "to_same_molecule_as";
pub const KEY_VEC_SPECIFIER_RANDOM: &str = "random";
