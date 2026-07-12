input_structure_file = "MD_300K_25ns-30ns_center_0001.pdb";
radius = 25; // angstroms.

detected_spin_position = centroid_over_serials([28,29]);
detected_spin_g_matrix = [2.0097, 2.0064, 2.0025];
detected_spin_g_x = diff(filter(tempo_n) , filter(tempo_o) );
detected_spin_g_y = diff(filter(tempo_c1) , filter(tempo_c19) );
number_timepoints = [101];

time_increments = [5e-8];
cluster_method = cce;
max_cluster_size = 3;
magnetic_field = 1.2; // T.
pulse_sequence = hahn;

orientation_grid = lebedev(14);

neighbor_cutoff_3_spin_hahn_mod_depth = 3.23e-5; 
neighbor_cutoff_3_spin_hahn_taylor_4 = 1e17; // (rad/s)^4.

#[group(tempo_h)]
  residues in [TEM];
  elements in [H];

#[spin_properties(tempo_h, 1H)]
  tunnel_splitting = 80e3; // Hz.

#[group(tempo_n)]
  residues in [TEM];
  elements in [N];

#[group(tempo_o)]
  residues in [TEM];
  elements in [O];

#[group(tempo_c1)]
  serials in [1];

#[group(tempo_c19)]
  serials in [19];
  
#[spin_properties(tempo_n, 14N)]
  hyperfine_coupling = [20,20,100]*1e6;
  hyperfine_x = diff(particle,bonded(tempo_o));
  hyperfine_y = diff(bonded(tempo_c1),bonded(tempo_c19));

  electric_quadrupole_coupling = [-0.28, -1.47, 1.75]*1e6;
  electric_quadrupole_x = diff(particle,bonded(tempo_o));
  electric_quadrupole_y = diff(bonded(tempo_c1),bonded(tempo_c19));

#[config]
write_orientation_signals = false;
write_clusters = false;


