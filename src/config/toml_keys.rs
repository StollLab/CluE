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
pub const KEY_OUT_CLU_SIGS: &str = "cluster_signals";
pub const KEY_OUT_BATH: &str = "bath";
pub const KEY_OUT_DET_SPIN: &str = "detected_spin";
pub const KEY_OUT_CLUSTERS: &str = "clusters";
pub const KEY_OUT_CONFIG: &str = "config";
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

// Group Keys
pub const KEY_NAME: &str = "name";

pub const KEY_DROP_PROB: &str = "drop_probability";

pub const KEY_ISO_COSUBSTITUTE: &str = "cosubstitute";


// Filter Keys
pub const KEY_SELECTION: &str = "selection";
pub const KEY_CELL_TYPE: &str = "cell_type";


pub const KEY_SELE_INDICES: &str = "indices";
pub const KEY_SELE_NOT_INDICES: &str = "not_indices";

pub const KEY_SELE_CELL_IDS: &str = "cell_ids";
pub const KEY_SELE_NOT_CELL_IDS: &str = "not_cell_ids";

pub const KEY_SELE_CHAIN_IDS: &str = "chain_ids";
pub const KEY_SELE_NOT_CHAIN_IDS: &str = "not_chain_ids";

pub const KEY_SELE_PRIMARY_CELL: &str = "primary_cell";

pub const KEY_SELE_ELEMENTS: &str = "elements";
pub const KEY_SELE_NOT_ELEMENTS: &str = "not_elements";

pub const KEY_SELE_SERIALS: &str = "serials";
pub const KEY_SELE_NOT_SERIALS: &str = "not_serials";

pub const KEY_SELE_NAMES: &str = "names";
pub const KEY_SELE_NOT_NAMES: &str = "not_names";

pub const KEY_SELE_RESIDUES: &str = "residues";
pub const KEY_SELE_NOT_RESIDUES: &str = "not_residues";

pub const KEY_SELE_RES_SEQ_NUMS: &str = "residue_sequence_numbers";
pub const KEY_SELE_NOT_RES_SEQ_NUMS: &str = "not_residue_sequence_numbers";

pub const KEY_SELE_ISOTOPES: &str = "isotpes";
pub const KEY_SELE_NOT_ISOTOPES: &str = "not_isotopes";

pub const KEY_SELE_BONDED_INDICES: &str = "bonded_indices";
pub const KEY_SELE_NOT_BONDED_INDICES: &str = "not_bonded_indices";

pub const KEY_SELE_WITHIN_DISTANCE: &str = "within_distance";
pub const KEY_SELE_NOT_WITHIN_DISTANCE: &str = "not_within_distance";

pub const KEY_SELE_BONDED_CHAIN_IDS: &str = "bonded_chain_ids";
pub const KEY_SELE_NOT_BONDED_CHAIN_IDS: &str = "not_bonded_chain_ids";

pub const KEY_SELE_BONDED_ELEMENTS: &str = "bonded_elements";
pub const KEY_SELE_NOT_BONDED_ELEMENTS: &str = "not_bonded_elements";

pub const KEY_SELE_BONDED_SERIALS: &str = "bonded_serials";
pub const KEY_SELE_NOT_BONDED_SERIALS: &str = "not_bonded_serials";

pub const KEY_SELE_BONDED_NAMES: &str = "bonded_names";
pub const KEY_SELE_NOT_BONDED_NAMES: &str = "not_bonded_names";

pub const KEY_SELE_BONDED_RESIDUES: &str = "bonded_residues";
pub const KEY_SELE_NOT_BONDED_RESIDUES: &str = "not_bonded_residues";

pub const KEY_SELE_BONDED_RES_SEQ_NUMS: &str = "bonded_residue_sequence_numbers";
pub const KEY_SELE_NOT_BONDED_RES_SEQ_NUMS: &str 
    = "not_bonded_residue_sequence_numbers";

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

pub const ALLOWED_KEYS: [&str;364] = [
  "abundance",
  "active",
  "auxiliary_signals",
  "axes",
  "bath",
  "bonded_elements",
  "bonded_chain_ids",
  "bonded_indices",
  "bonded_names",
  "bonded_residues",
  "bonded_residue_sequence_numbers",
  "bonded_serials",
  "c3_tunnel_splitting",
  "cell_ids",
  "cell_type",
  "chain_ids",
  "cosubstitute",
  "coupling",
  "coupling_xx_yy",
  "connect_exchange_groups",
  "clash_distance",
  "clash_distance_pbc",
  "cluster_batch_size",
  "cluster_method",
  "cluster_source",
  "clusters",
  "cluster_signals",
  "config",
  "delta_hyperfine_zz",
  "delta_zeeman",
  "density_matrix",
  "detected_spin",
  "detection_operator",
  "distance",
  "drop_probability",
  "electric_quadrupole",
  "elements",
  "ensemble_cce",
  "exchange_coupling",
  "exchange_groups",
  "file",
  "from",
  "from_bonded_to",
  "from_same_molecule_as",
  "g_matrix",
  "grid",
  "groups",
  "hahn_mod_depth",
  "hahn_taylor_4",
  "hyperfine",
  "identity",
  "indices",
  "info",
  "input_structure_file",
  "isotpes",
  "kmeans_size",
  "lebedev",
  "magnetic_field",
  "matrix",
  "max_cell_size",
  "max_cluster_size",
  "max_spins",
  "mean_fields",
  "methyl_partitions",
  "min_cell_size",
  "multiplicity",
  "name",
  "names",
  "not_bonded_elements",
  "not_bonded_chain_ids",
  "not_bonded_indices",
  "not_bonded_names",
  "not_bonded_residues",
  "not_bonded_residue_sequence_numbers",
  "not_bonded_serials",
  "not_cell_ids",
  "not_chain_ids",
  "not_elements",
  "not_indices",
  "not_isotopes",
  "not_names",
  "not_residues",
  "not_residue_sequence_numbers",
  "not_serials",
  "not_within_distance",
  "number",
  "number_runs",
  "number_timepoints",
  "number_timepoints2",
  "orientation_grid",
  "orientations",
  "orientation_signals",
  "output",
  "output_directory",
  "pair_cutoffs",
  "partitioning",
  "partition_table",
  "pdb_model_index",
  "point_dipole_perpendicular",
  "populations",
  "position",
  "primary_cell",
  "pulse_sequence",
  "pulses",
  "radius",
  "random",
  "replicate_unit_cell",
  "residues",
  "residue_sequence_numbers",
  "rng_seed",
  "run_name",
  "run_in_parallel",
  "sans_spin_signals",
  "selection",
  "serials",
  "singles",
  "structure_pdb",
  "tau_increments",
  "tau2_increments",
  "temperature",
  "tensors",
  "thermal",
  "to",
  "to_bonded_to",
  "to_same_molecule_as",
  "transition",
  "unit_of_distance",
  "unit_of_energy",
  "unit_of_magnetic_field",
  "unit_of_time",
  "values",
  "vector",
  "vector_grid",
  "within_distance",
  "x",
  "y",
  "z",
  "zerofield",
  "e",
  "1H",
  "2H",
  "3H",
  "3He",
  "4He",
  "6Li",
  "7Li",
  "9Be",
  "10B",
  "11B",
  "12C",
  "13C",
  "14N",
  "15N",
  "16O",
  "17O",
  "19F",
  "20Ne",
  "21Ne",
  "22NA",
  "23Na",
  "24Mg",
  "25Mg",
  "27Al",
  "28Si",
  "29Si",
  "31P",
  "32S",
  "33S",
  "35Cl",
  "36Cl",
  "37Cl",
  "39Ar",
  "40Ar",
  "39K",
  "40K",
  "41K",
  "40Ca",
  "41Ca",
  "43Ca",
  "45Sc",
  "47Ti",
  "48Ti",
  "49Ti",
  "50V",
  "51V",
  "52Cr",
  "53Cr",
  "53Mn",
  "55Mn",
  "56Fe",
  "57Fe",
  "59Co",
  "60Co",
  "58Ni",
  "61Ni",
  "63Cu",
  "65Cu",
  "64Zn",
  "67Zn",
  "69Ga",
  "71Ga",
  "73Ge",
  "74Ge",
  "75As",
  "77Se",
  "78Se",
  "79Se",
  "79Br",
  "81Br",
  "83Kr",
  "84Kr",
  "85Kr",
  "85Rb",
  "87Rb",
  "87Sr",
  "88Sr",
  "89Y" ,
  "90Zr",
  "91Zr",
  "93Nb",
  "95Mo",
  "97Mo",
  "98Mo",
  "99Tc",
  "99Ru",
  "101Ru",
  "102Ru",
  "103Rh",
  "105Pd",
  "106Pd",
  "107Ag",
  "109Ag",
  "111Cd",
  "113Cd",
  "114Cd",
  "113In",
  "115In",
  "115Sn",
  "117Sn",
  "119Sn",
  "120Sn",
  "121Sb",
  "123Sb",
  "125Sb",
  "123Te",
  "125Te",
  "130Te",
  "127I",
  "129I",
  "129Xe",
  "131Xe",
  "132Xe",
  "133Cs",
  "134Cs",
  "135Cs",
  "137Cs",
  "133Ba",
  "134Ba",
  "137Ba",
  "138Ba",
  "137La",
  "138La",
  "139La",
  "140Ce",
  "141Pr",
  "143Nd",
  "144Nd",
  "145Nd",
  "141Pm",
  "147Sm",
  "149Sm",
  "151Sm",
  "152Sm",
  "151Eu",
  "152Eu",
  "153Eu",
  "154Eu",
  "155Eu",
  "155Gd",
  "1566Gd",
  "157Gd",
  "157Tb",
  "159Tb",
  "160Tb",
  "161Dy",
  "163Dy",
  "165Dy",
  "165Ho",
  "166Er",
  "167Er",
  "169Tm",
  "171Tm",
  "171Yb",
  "173Yb",
  "174Yb",
  "173Lu",
  "174Lu",
  "175Lu",
  "176Lu",
  "177Hf",
  "179Hf",
  "180Hf",
  "180Ta",
  "181Ta",
  "183W",
  "184W",
  "185Rh",
  "187Rh",
  "187Os",
  "189Os",
  "192Os",
  "191Ir",
  "193Ir",
  "194Pt",
  "195Pt",
  "197Au",
  "199Hg",
  "201Hg",
  "202Hg",
  "203Tl",
  "204Tl",
  "205Tl",
  "207Pb",
  "208Pb",
  "207Bi",
  "209Bi",
  "209Po",
  "210Po",
  "210AT",
  "222Rn",
  "223Fr",
  "226Ra",
  "227Ac",
  "229Th",
  "232Th",
  "231Pa",
  "235U",
  "238U",
  "237Np",
  "239Pu",
  "243Am",
  "247Cm",
  "247Bk",
  "251Cf",
  "252Es",
  "257Fm",
  "258Md",
  "259No",
  "266Lr",
  "267Rf",
  "268Db",
  "269Sg",
  "278Bh",
  "270Hs",
  "278Mt",
  "281Ds",
  "282Rg",
  "285Cn",
  "286Nh",
  "289Fl",
  "289Mc",
  "293Lv",
  "294Ts",
  "294Og",
];
