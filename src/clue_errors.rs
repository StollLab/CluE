use crate::config::toml_keys::*;
use serde::Serialize;
use std::fmt;

/// `CluEError` contains the possible errors.
// After adding a new error, please update fmt() below as well.
#[derive(PartialEq,Debug,Clone,Serialize)]
pub enum CluEError{
  AllSignalsNotSameLength(String),
  AllVectorsNotSameLength(String),
  BondsAreNotDefined(String),
  BothStevensAndIST,
  CannotAddPointToGrid(usize,usize),
  CannotConvertSerialToIndex(u32),
  CannotCreateDir(String,String),
  CannotDiagonalizeHamiltonian(String),
  CannotDiagonalizeOperator(String),
  CannotExpandBlockClusters,
  CannotFindCellID(usize),
  CannotFindSpinOp(String),
  CannotFindParticleForRefIndex(usize),
  CannotFindBathIndexFromRefIndex(usize),
  CannotFindRefIndexFromBathIndex(usize),
  CannotFindRefIndexFromNthActive(usize),
  CannotNormalizeVector,
  CannotCreateFile(String,String),
  CannotOpenFile(String,String),
  CannotParseCellType(String),
  CannotParseClusterMethod(String),
  CannotParseElement(String),
  CannotParseGridEntry(String,usize,String),
  CannotParseLine(String),
  CannotParseLineInFile(String,usize,String),
  CannotParseSignalEntry(String,usize,String),
  CannotParseIsotope(String),
  CannotParseOrientations(String),
  CannotParsePartitioningMethod(String),
  CannotParsePulseSequence(String),
  CannotParseCosubstitute(String,String),
  CannotParseSpinOp(String),
  CannotPruneClustersMisMatchedSizes,
  CannotReadGrid(String,String),
  CannotReadTOML(String),
  CannotReadTOMLFile(String),
  CannotReadFile(String,String),
  CannotSampleBinomialDistribution(usize,f64),
  CannotTakeTrace(String),
  CannotWriteCSV(String,String),
  CannotWriteFile(String,String),
  CentroidGroupHasDistanceSelection(String),
  CentroidOfNoParticles,
  ClusterFileContainsNoHeader(String),
  ClusterHasNoSignal(String),
  ClusterLineFormatError(String),
  ClusterTOMLCannotRead(String),
  ClusterTOMLIncorrectNumberOfClusters,
  ClusterTOMLNoClusters,
  ClusterTOMLNoNumberClusters,
  PulseSequenceNotSupported(String,String),
  DensityMatrixTraceVanishes(String,f64),
  DetectedSpinDoesNotHaveAnActiveIndex,
  DetectedSpinPositionNotFinite(String),
  DetectedSpinSelectsNoParticles(String),
  DuplicateSerial(u32),
  ExpectedClusterSetWithNSizes(usize,usize),
  ExpectedTOMLArray(String),
  ExpectedTOMLFloat(String),
  ExpectedTOMLUInt(String),
  FailedKMeans,
  FiltersOverlap(String,String),
  GroupHasNoSelection(String),
  InvalidGroBox(String),
  IncompleteConfig(Vec::<String>),
  IncorrectNumberOfCellOffsets(usize,usize),
  InvalidAxes,
  InvalidClusterPartitionKey,
  InvalidDetectionFrame(String),
  InvalidDetectionOperator(String),
  InvalidIST(String),
  InvalidKMeansIterCap,
  InvalidKMeansSize,
  InvalidNumberOfKMeansClusters(usize,usize),
  InvalidPartitionTableEntry(usize,usize,String),
  InvalidPulseStep(usize,Box<CluEError>),
  InvalidAbundanceTolerance(String,f64),
  InvalidDropProbability(String,f64),
  InvalidIsotopeAbundance(String,String,f64),
  InvalidSelectionElement(String,usize,String,String),
  InvalidSpinGridWeights(f64),
  InvalidSpinGridWeightTolerance(f64),
  InvalidSpinMultiplicity(usize),
  InvalidOrientationGridWeightTolerance(f64),
  IsotopeAbundancesDoNotSumToOne(String,String,f64,f64),
  InvalidStevens(String),
  InvalidStevensCoefficientShape(usize,usize,usize),
  InvalidTOMLKeys(Vec::<String>),
  InvalidTOMLValue(String,String,String),
  LengthMismatchTimepointsIncrements(String,usize,usize),
  MeanFieldECCENotImplemented,
  MismatchedGroupNames(String,String),
  MissingDMatrix(usize),
  MissingFieldInCSVFile(String,String),
  MissingFilter(String),
  MissingGroupName,
  MissingHeader(usize,usize,String),
  MissingTOMLKey(String,String),
  MutuallyExclusiveKeys(String,String),
  NANTensorBathDipoleDipole(String,String,String,String),
  NANTensorBathZeeman(String,String),
  NANTensorDetectedZeeman,
  NANTensorExchangeCoupling(String,String,String,String),
  NANTensorHyperfine(String,String),
  NANTensorQuadrupole(String,String),
  NANTensorZerofield(String,String),
  NeighborListAndPartitionTableAreNotCompatible,
  NoApplyPBC,
  NoClusterBatchSize,
  NoClusterMethod,
  NoClusterSource,
  NoClusterDensityMatrixMethod,
  NoDetectionFrame,
  NoDetectedSpinDensityMatrix,
  NoDetectedSpinDetectionOperator,
  NoDetectedSpinIdentity,
  NoDetectedSpinMultiplicity,
  NoDetectedSpinPosition,
  NoDetectedSpinTransition,
  NoDetectedParticle,
  NoEnsembleCCE,
  NoGMatrixSpecifier,
  NoInputFile,
  NoKMeansIterCap,
  NoKMeansSize,
  NoMagneticField,
  NoMeanFields,
  NoModelIndex,
  NoMaxClusterSize,
  NoMaxISTRank,
  NoNumberSystemInstances,
  NoOrientationGrid,
  OrientationGridWeightsDoNotSumToOne(f64,f64),
  NoPartitioningMethod,
  NoPulseSequence,
  NoRadius,
  NoRunInParallel,
  NoDensityMatrixWithMultiplicity(usize),
  NoDetOpWithMultiplicity(usize),
  NoPulseOpWithMultiplicity(String,usize),
  NoSpinOpForClusterSize(usize,usize),
  NoSpinOpWithMultiplicity(usize),
  NoStevensOp(i32,i32),
  NoStructureFile,
  NotA3DRotationMatrix(String),
  NotA3DVector(usize),
  NotALebedevGrid(usize),
  NotAProperSubset(String,String),
  NoTemperature,
  NoTensorSpecifier,
  NoTensorValues,
  NoTimeAxis,
  NoTimeIncrements,
  NoTimeIncrements2,
  NoTimepoints,
  NoTimepoints2,
  NoUnitOfDistance,
  NoUnitOfEnergy,
  NoUnitOfMagneticField,
  NoUnitOfTime,
  ParticlesClash(String,String,String,String,f64,f64),
  Python(String),
  ParticleIsNotActive(usize),
  PropagatorAndDensityNotSameDimension(usize,usize),
  PartitionIsIncomplete,
  ResidueNumberAboveMaxInPDB(String,String,usize),
  SecondaryFilterRequiresAnIndex(String),
  SpaceGroupNotP1(String),
  SpinGridWeightsDoNotSumToOne(f64,f64),
  SpinAlreadyPartitioned(usize,i32),
  StructureLargerThanCell(String,f64),
  TensorNotSet(usize),
  ThermalizationRequiresEnsembleCCE,
  TooManyAtomsInPDB(String,usize),
  TooManyResiduesInPDB(String,String,usize),
  TOMLArrayContainsMultipleTypes,
  TOMLArrayDoesNotSpecifyATensor,
  TOMLArrayDoesNotSpecifyAVector,
  TOMLArrayIsEmpty,
  TOMLArrayIsNotAMatrix,
  TOMLArrayIsNotAComplexNumber,
  TomlPulseStepSpecifier(String),
  UnavailableSpinOp(usize,usize),
  UnassignedCosubstitutionGroup(usize),
  UnequalLengths(String,usize,String,usize),
  UnexpectedPulseStep(usize,String),
  UnknownPulse(String),
  UnrecognizedOption(String),
  UnrecognizedVectorSpecifier(String),
  UnrecognizedUnit(String,String,String),
  VectorSpecifierDoesNotSpecifyUniqueVector(String),
  WrongBatchSignalLength(String,usize,usize),
  WrongClusterSizeForAnalyticCCE(usize),
  WrongDelayDimension(usize),
  WrongNumberOfBatchSignals(String,usize,usize),
  WrongNumberOfAxes(usize,usize),
  WrongNumberOfCellEdges(usize),
  WrongTauDimension(usize,usize),
  ZeroClusterBatchSize,
  ZeroMaxClusterSize,
}

impl fmt::Display for CluEError{
  // This function defines the error messages.
  fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
    match self{

      CluEError::AllSignalsNotSameLength(filename) => write!(f,
          "for \"{}\", signals must all have the same length", filename),

      CluEError::AllVectorsNotSameLength(filename) => write!(f,
          "cannot write \"{}\": all data columns must have the same length",
          filename),

      CluEError:: BondsAreNotDefined(label) => write!(f,
          "{} needs bonds, but the particle has none; bonds are read \
from the CONECT records of the structure file", label),

      CluEError:: BothStevensAndIST => write!(f,
          "{} and {} are mutually exclusive",KEY_IST_COEF,KEY_STEVENS_COEF),

      CluEError::CannotAddPointToGrid(point_dim, grid_dim) => write!(f,
          "cannot add {}D point to {}D grid",point_dim, grid_dim),

      CluEError::CannotConvertSerialToIndex(serial) => write!(f,
          "cannot convert serial id, {}, to an index",serial),

      CluEError::CannotDiagonalizeHamiltonian(matrix) => write!(f,
          "cannot diagonalize the Hamiltonian, {}",matrix),

      CluEError::CannotDiagonalizeOperator(matrix) => write!(f,
          "cannot diagonalize the operator, {}",matrix),

      CluEError::CannotCreateDir(path,reason) => write!(f,
          "cannot create directory \"{}\": {}",path,reason),

      CluEError::CannotExpandBlockClusters => write!(f,
          "cannot expand block clusters"),

      CluEError::CannotFindCellID(idx) => write!(f,
          "cannot determine cell id for particle {}",idx),

      CluEError::CannotFindSpinOp(sop) => write!(f,
          "cannot find spin operator \"{}\"",sop),

      CluEError::CannotFindParticleForRefIndex(ref_index) => write!(f,
          "no active bath particle matches reference index \"{}\"",
          ref_index),

      CluEError::CannotFindBathIndexFromRefIndex(ref_index) => write!(f,
          "cannot find bath index for reference index \"{}\"", ref_index),

      CluEError::CannotFindRefIndexFromBathIndex(bath_index) => write!(f,
          "cannot find reference index for bath index \"{}\"", bath_index),

      CluEError::CannotFindRefIndexFromNthActive(n) => write!(f,
          "cannot find reference index for the {}th active particle", n),

      CluEError::CannotNormalizeVector => write!(f,
          "cannot normalize vector"),

      CluEError::CannotCreateFile(file,reason) => write!(f,
          "cannot create \"{}\": {}", file, reason),

      CluEError::CannotOpenFile(file,reason) => write!(f,
          "cannot open \"{}\": {}", file, reason),

      CluEError::CannotParseCellType(cell_type) => write!(f,
          "cannot parse \"{}\" as cells; the options are \"all\", \
\"primary\", and \"extra\"", cell_type),

      CluEError::CannotParseClusterMethod(method) => write!(f,
          "cannot parse \"{}\" as a cluster method; the options are \
\"cce\", \"gcce\", \"appa\", and \"lce\"", method),

      CluEError::CannotParseElement(element) => write!(f,
          "cannot parse \"{}\" as an element", element),

      CluEError::CannotParseGridEntry(filename,line,entry) => write!(f,
          "cannot read grid from \"{}\": \"{}\" on line {} is not a number",
          filename, entry, line),

      CluEError::CannotParseOrientations(ori) => write!(f,
          "cannot parse orientations.grid = \"{}\"; the options are \
\"lebedev\", \"random\", \"file\", \"vector\", \"circle\", and \
\"vector_grid\"", ori),

      CluEError::CannotParseLine(line) => write!(f,
          "cannot parse line \"{}\"", line),

      CluEError::CannotParseLineInFile(file,line_number,line) => write!(f,
          "cannot parse line {} of \"{}\": \"{}\"", line_number, file, line),

      CluEError::CannotParseSignalEntry(file,row,entry) => write!(f,
          "cannot read signals from \"{}\": \"{}\" in data row {} is not \
a complex number", file, entry, row),

      CluEError::CannotParseIsotope(isotope) => write!(f,
          "cannot parse \"{}\" as an isotope", isotope),

      CluEError::CannotParsePartitioningMethod(part_method) => write!(f,
          "cannot parse partitioning method \"{}\"; the options are \
\"{}\", \"{}\", \"{}\", and \"{}\"", part_method, 
          KEY_PARTITION_PARTICLE, KEY_PARTITION_EX_GROUPS, 
          KEY_PARTITION_KMEANS, KEY_PARTITION_RESTRICTED_KMEANS),

      CluEError::CannotParsePulseSequence(seq) => write!(f,
          "cannot parse \"{}\" as a pulse sequence; the options are \
\"hahn\", \"fid\", \"free_evolution\", \"refocused_hahn_echo\", \
\"cp-N\" for an integer N >= 0, or an array of pulse steps", seq),

      CluEError::CannotParseCosubstitute(group,cosubstitute) => write!(f,
          "group \"{}\": cannot parse cosubstitute = \"{}\"; the options \
are \"same_molecule\", \"group\", and \"particle\" (\"bonded\" is not \
supported)", group, cosubstitute),

      CluEError::CannotParseSpinOp(s) => write!(f,
          "cannot parse \"{}\" as a spin operator; the options are \
\"Sx\", \"Sy\", \"Sz\", \"S+\", \"S-\", and \"S^2\" \
(or the same in lowercase)", s),

      CluEError::CannotSampleBinomialDistribution(n,p) => write!(f,
          "cannot sample from the binomial distribution B(n={},p={})",
          n,p),

      CluEError::CannotTakeTrace(matrix) => write!(f,
          "cannot take the trace of {}",matrix),

      CluEError::CannotPruneClustersMisMatchedSizes => write!(f,
          "cannot prune clusters because clusters_to_keep has the wrong size"),

      CluEError::CannotReadGrid(filename,reason) => write!(f,
          "cannot read grid from \"{}\": {}; \
the grid should be specified as a csv file with one column per dimension, \
followed by a column for the weights", filename,reason),

      CluEError::CannotReadTOML(msg) => write!(f,
          "cannot parse TOML: {}",msg),

      CluEError::CannotReadTOMLFile(msg) => write!(f,
          "cannot parse the config as TOML: {}",msg),

      CluEError::CannotReadFile(file,reason) => write!(f,
          "cannot read \"{}\": {}", file, reason),

      CluEError::ClusterFileContainsNoHeader(file) => write!(f,
          "cluster file \"{}\" does not contain the correct header: \
\"#[clusters, number_clusters = [N1,N2,...Nn] ]\", where Ni is the number \
of clusters of size i in the file, and the list runs from clusters of \
size 1 to clusters of size n", file),

      CluEError::ClusterHasNoSignal(cluster) => write!(f,
          "expected cluster {} to have a signal, but found none", cluster),

      CluEError::ClusterLineFormatError(line) => write!(f,
          "cluster line \"{}\" is not formatted correctly", line),

      CluEError::ClusterTOMLCannotRead(msg) => write!(f,
          "cannot read cluster toml: \"{}\"", msg),

      CluEError::ClusterTOMLNoClusters => write!(f,
          "no clusters in clusters.toml"),

      CluEError::ClusterTOMLNoNumberClusters => write!(f,
          "no number_clusters in clusters.toml"),

      CluEError::ClusterTOMLIncorrectNumberOfClusters=> write!(f,
          "in clusters.toml, number_clusters is incorrect"),

      CluEError::CannotWriteCSV(file,reason) => write!(f,
          "cannot write CSV file \"{}\": {}", file, reason),

      CluEError::CannotWriteFile(file,reason) => write!(f,
          "cannot write to \"{}\": {}", file,reason),

      CluEError::CentroidOfNoParticles => write!(f,
          "cannot find the centroid of an empty set of particles"),

      CluEError::CentroidGroupHasDistanceSelection(label) => write!(f,
          "group \"{}\" is used to place the detected spin and cannot use \
include.distance, since distances are measured from the detected spin",
          label),

      CluEError::PulseSequenceNotSupported(method,seq) => write!(f,
        "cluster method \"{}\" only supports pulse_sequence.sequence = \
\"{}\"",
        method,seq),

      CluEError::DuplicateSerial(srl) => write!(f,
          "serial \"{}\" has multiple entries", srl),

      CluEError::DensityMatrixTraceVanishes(trace,value) => write!(f,
          "cannot normalize the thermal density matrix: |{}| = {:e} is \
below 1e-12; check the temperature and the energy scale of the Hamiltonian",
          trace, value),

      CluEError::DetectedSpinDoesNotHaveAnActiveIndex => write!(f,
          "the detected spin (reference index 0) is always active and has no \
position among the active bath particles"),

      CluEError::DetectedSpinPositionNotFinite(position) => write!(f,
          "the detected spin position, {} angstroms, is not finite; check \
detected_spin.position and the coordinates in the structure file",
          position),

      CluEError::DetectedSpinSelectsNoParticles(selection) => write!(f,
          "the detected spin position is the centroid of {}", selection),

      CluEError::ExpectedClusterSetWithNSizes(n_exp,n_act) => write!(f,
          "expected a cluster set with {} sizes, but got {} sizes",
          n_exp, n_act),

      CluEError::ExpectedTOMLArray(type_str) => write!(f,
          "expected array, but got {}",type_str),

      CluEError::ExpectedTOMLFloat(type_str) => write!(f,
          "expected float, but got {}",type_str),

      CluEError::ExpectedTOMLUInt(type_str) => write!(f,
          "expected int > 0, but got {}",type_str),

      CluEError::FailedKMeans => write!(f,
          "k-means failed"),

      CluEError::FiltersOverlap(label0,label1) => write!(f,
          "groups \"{}\" and \"{}\" overlap: \
particles must not reside in more than one group",label0,label1),

      CluEError::GroupHasNoSelection(label) => write!(f,
          "group \"{}\" has no include or exclude table, so it selects no \
atoms; include = {{}} selects all atoms (of the cells set by cells)",
          label),

      CluEError::InvalidGroBox(edges) => write!(f,
          "the unit cell edges {} (nm) cannot be written in GRO format, \
which needs the first edge along x and the second edge in the xy plane",
          edges),

      CluEError::IncompleteConfig(problems) => write!(f,
          "the input is incomplete:\n  {}", problems.join("\n  ")),

      CluEError::IncorrectNumberOfCellOffsets(n,n_ref)=> write!(f,
          "expected {} cell offsets, but {} were provided",n_ref, n),

      CluEError::InvalidAxes => write!(f,
          "invalid axes; the axes must be among x, y, and z"),

      CluEError::InvalidClusterPartitionKey => write!(f,
          "invalid cluster partition key"),

      CluEError::InvalidDetectionFrame(frame) => write!(f,
          "cannot parse detected_spin.frame = \"{}\"; the options are \
\"{}\" and \"{}\"", frame, KEY_DET_FRAME_EIGEN, KEY_DET_FRAME_ZEEMAN),

      CluEError::InvalidDetectionOperator(reason) => write!(f,
          "invalid detected-spin operator (detection operator, density matrix, \
or pulse): {}",reason),

      CluEError::InvalidIST(ist) => write!(f,
          "invalid {} \"{}\"",KEY_IST_COEF,ist),

      CluEError::InvalidPartitionTableEntry(serial,block,reason) => write!(f,
          "cluster_expansion.clustering.partition_table entry (serial {}, \
block {}): {}",
          serial,block,reason),

      CluEError::InvalidPulseStep(index,error) => write!(f,
          "pulse_sequence.sequence[{}]: {}",index,error),

      CluEError::InvalidKMeansIterCap => write!(f,
          "cluster_expansion.clustering.kmeans_iter_cap must be at least 1"),

      CluEError::InvalidKMeansSize => write!(f,
          "cluster_expansion.clustering.kmeans_size must be at least 1"),

      CluEError::InvalidNumberOfKMeansClusters(k,n) => write!(f,
          "k-means cannot split {} points into {} clusters; the number of \
clusters must be between 1 and the number of points",n,k),

      CluEError::InvalidAbundanceTolerance(group,tolerance) => write!(f,
          "group \"{}\": abundance_tolerance = {} must be at least 0 and \
less than 1", group, tolerance),

      CluEError::InvalidDropProbability(group,p) => write!(f,
          "group \"{}\": drop_probability = {} is outside the interval \
[0, 1]", group, p),

      CluEError::InvalidIsotopeAbundance(group,isotope,abundance) => write!(f,
          "group \"{}\": {}.abundance = {} is not a non-negative number",
          group, isotope, abundance),

      CluEError::IsotopeAbundancesDoNotSumToOne(group,abundances,sum,tol) 
        => write!(f,
          "group \"{}\": the isotope abundances ({}) sum to {}, which \
differs from 1 by more than abundance_tolerance = {}; list every isotope, \
for example 1H.abundance = 0.7 and 2H.abundance = 0.3", 
          group, abundances, sum, tol),

      CluEError::InvalidSelectionElement(key,index,expected,found)
        => write!(f,
          "{}[{}] must be {}, but got {}",
          key, index, expected, found),

      CluEError::InvalidSpinGridWeights(sum) => write!(f,
          "the weights of the detected spin grid sum to {}, but they must \
sum to a positive number", sum),

      CluEError::InvalidSpinGridWeightTolerance(tolerance) => write!(f,
          "detected_spin.weight_tolerance = {} must be at least 0 and less \
than 1", tolerance),

      CluEError::InvalidOrientationGridWeightTolerance(tolerance) => write!(f,
          "orientations.weight_tolerance = {} must be at least 0 and less \
than 1", tolerance),

      CluEError::InvalidSpinMultiplicity(n) => write!(f,
          "invalid spin multiplicity \"{}\"",n),

      CluEError::InvalidStevens(ist) => write!(f,
          "invalid {} \"{}\"",KEY_STEVENS_COEF,ist),

      CluEError::InvalidStevensCoefficientShape(index,n_rows,n_cols) 
        => write!(f,
          "Stevens coefficient entry {} is a {}×{} matrix, but each entry \
must have one column", index, n_rows, n_cols),

      CluEError::InvalidTOMLKeys(problems) => write!(f,
          "the input contains keys that are not allowed where they \
appear:\n  {}", problems.join("\n  ")),

      CluEError::InvalidTOMLValue(key,expected,found) => write!(f,
          "{} must be {}, but got {}",key,expected,found),

      CluEError::LengthMismatchTimepointsIncrements(suffix,n_dts,dts) 
        => write!(f,
          "pulse_sequence.number_timepoints{} has {} entries, but \
pulse_sequence.tau{}_increments has {}; they must have the same length", 
          suffix, n_dts, suffix, dts),

      CluEError::MeanFieldECCENotImplemented => write!(f,
          "mean field averaging is not implemented for ensemble CCE \
(cluster_expansion.mean_fields = true with cluster_expansion.ensemble = true)"),

      CluEError::MismatchedGroupNames(name0,name1) => write!(f,
          "group names \"{}\" and \"{}\" do not match",name0,name1),

      CluEError::MissingDMatrix(l) => write!(f,
          "no D matrix for l = {}",l),

      CluEError::MissingFilter(label) => write!(f,
          "no group with label \"{}\"",label),

      CluEError::MissingFieldInCSVFile(field,file_name) => write!(f,
          "CSV file \"{}\" has no column \"{}\"",file_name,field),

      CluEError::MissingTOMLKey(key,required_by) => write!(f,
          "{} is required for {}",key,required_by),

      CluEError::MutuallyExclusiveKeys(key0,key1) => write!(f,
          "{} and {} are mutually exclusive",key0,key1),

      CluEError::MissingGroupName => write!(f,
          "no group name specified"),

      CluEError::MissingHeader(n_headers, n_cols,filename) => write!(f,
          "in \"{}\", every entry must have a header, \
but there are {} headers and {} columns of data",filename,n_headers,n_cols),

      CluEError::NANTensorBathDipoleDipole(idx0,isotope0,idx1,isotope1) 
        => write!(f,
          "the dipole-dipole tensor for {} ({}) and {} ({}) contains NaNs",
          idx0,isotope0,idx1,isotope1),

      CluEError::NANTensorBathZeeman(idx,isotope) => write!(f,
          "the Zeeman tensor for {} ({}) contains NaNs",idx,isotope),

      CluEError::NANTensorDetectedZeeman => write!(f,
          "the Zeeman tensor for the detected particle contains NaNs"),

      CluEError::NANTensorExchangeCoupling(idx0,isotope0,idx1,isotope1) 
        => write!(f,
          "the exchange tensor for {} ({}) and {} ({}) contains NaNs",
          idx0,isotope0,idx1,isotope1),

      CluEError::NANTensorHyperfine(idx,isotope) => write!(f,
          "the hyperfine tensor for {} ({}) contains NaNs",idx,isotope),

      CluEError::NANTensorQuadrupole(idx,isotope) => write!(f,
          "the electric quadrupole tensor for {} ({}) contains NaNs",
          idx,isotope),

      CluEError::NANTensorZerofield(idx,isotope) => write!(f,
          "the zero-field tensor for {} ({}) contains NaNs",
          idx,isotope),

      CluEError::NeighborListAndPartitionTableAreNotCompatible  => write!(f,
          "neighbor list and partition table are incompatible"),

      CluEError::NoApplyPBC => write!(f,
          "structure.replicate_unit_cell is not set (true, false, or \
[n1, n2, n3])"),

      CluEError::NoClusterBatchSize => write!(f,
          "batch size for clusters not specified"),

      CluEError::NoClusterMethod => write!(f,
          "no cluster method specified"),

      CluEError::NoClusterSource => write!(f,
          "no cluster source specified"),

      CluEError::NoClusterDensityMatrixMethod=> write!(f,
          "no density matrix method specified"),

      CluEError::NoDetectionFrame=> write!(f,
          "no frame for the detected spin"),

      CluEError::NoDetectedSpinDensityMatrix=> write!(f,
          "no density matrix for the detected spin"),

      CluEError::NoDetectedSpinDetectionOperator=> write!(f,
          "no detection operator"),

      CluEError::NoDetectedSpinIdentity => write!(f,
          "detected_spin.identity is not set"),

      CluEError::NoDetectedSpinMultiplicity => write!(f,
          "the detected spin multiplicity is not set"),

      CluEError::NoDetectedSpinPosition => write!(f,
          "detected_spin.position is not set"),

      CluEError::NoDetectedSpinTransition => write!(f,
          "detected_spin.transition is not set"),

      CluEError::NoDetectedParticle => write!(f,
          "the structure has no detected spin; the detected spin is added \
when the structure is built from the config"),

      CluEError::NoEnsembleCCE => write!(f,
          "cluster_expansion.ensemble is not set"),

      CluEError::NoGMatrixSpecifier => write!(f,
          "no g-matrix specifier"),

      CluEError::NoInputFile => write!(f,
          "no input file"),

      CluEError::NoKMeansIterCap => write!(f,
          "k-means requires cluster_expansion.clustering.kmeans_iter_cap to \
be set"),

      CluEError::NoKMeansSize => write!(f,
          "k-means requires cluster_expansion.clustering.kmeans_size to be \
set"),

      CluEError::NoMagneticField => write!(f,
          "please specify the applied magnetic field"),

      CluEError::NoMeanFields => write!(f,
          "cluster_expansion.mean_fields is not set"),

      CluEError::NoModelIndex => write!(f,
          "PDB model not selected"),

      CluEError::NoMaxClusterSize => write!(f,
          "maximum cluster size not set"),

      CluEError::NoMaxISTRank => write!(f,
          "maximum irreducible spherical tensor rank not set"),

      CluEError::NoNumberSystemInstances => write!(f,
          "run.number_runs is not set"),

      CluEError::NoOrientationGrid => write!(f,
          "orientations.grid is not set"),

      CluEError::OrientationGridWeightsDoNotSumToOne(sum,tolerance) => write!(f,
          "the weights of the orientation grid sum to {}, which differs \
from 1 by more than orientations.weight_tolerance = {}", sum, tolerance),

      CluEError::NoPartitioningMethod => write!(f,
          "no partitioning method is defined"),

      CluEError::NoPulseSequence => write!(f,
          "no pulse sequence is defined"),

      CluEError::NoRadius => write!(f,
          "structure.radius is not set"),

      CluEError::NotA3DRotationMatrix(matrix) => write!(f,
          "the matrix\n{}\ndoes not correspond to a valid rotation matrix", 
          matrix),

      CluEError::NotA3DVector(dim) => write!(f,
          "vector is {}-dimensional, not 3-dimensional", 
          dim),

      CluEError::NotALebedevGrid(n_ori) => write!(f,
          "orientations.number = {} is not a valid Lebedev grid size; \
please choose one of {{6, 14, 26, 38, 50, 74, 86, 110, 146, 170, \
194, 230, 266, 302, 350, 434, 590, 770, 974, 1202, \
1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810}}",
          n_ori),

      CluEError::NotAProperSubset(subcluster, cluster) => write!(f,
          "{} is not a proper subset of {}", subcluster, cluster),

      CluEError::NoTemperature => write!(f,
          "conditions.temperature is not set"),

      CluEError::NoTensorSpecifier => write!(f,
          "no specifier for tensor"),

      CluEError::NoTensorValues => write!(f,
          "no values specified for tensor"),

      CluEError::NoTimeAxis => write!(f,
          "the time-axis has not been built"),

      CluEError::NoTimeIncrements => write!(f,
          "pulse_sequence.tau_increments is not set"),

      CluEError::NoTimeIncrements2 => write!(f,
          "pulse_sequence.tau2_increments is not set"),

      CluEError::NoTimepoints => write!(f,
          "pulse_sequence.number_timepoints is not set; please specify how \
many timepoints there are for each tau increment"),

      CluEError::NoTimepoints2 => write!(f,
          "pulse_sequence.number_timepoints2 is not set; please specify how \
many timepoints there are for each tau2 increment"),

      CluEError::NoUnitOfDistance => write!(f,
          "please specify a unit of distance"),

      CluEError::NoUnitOfEnergy => write!(f,
          "please specify a unit of energy"),

      CluEError::NoUnitOfMagneticField => write!(f,
          "please specify a unit of magnetic field"),

      CluEError::NoUnitOfTime => write!(f,
          "please specify a unit of time"),

      CluEError::NoRunInParallel => write!(f,
          "run.parallel is not set"),

      CluEError::NoDensityMatrixWithMultiplicity(spin_multiplicity) => write!(f,
          "no bath-spin-{} density matrices are built for the detected spin", 
          (*spin_multiplicity as f64 - 1.0)/2.0),

      CluEError::NoDetOpWithMultiplicity(spin_multiplicity) => write!(f,
          "no bath-spin-{} operators are built for the detected spin", 
          (*spin_multiplicity as f64 - 1.0)/2.0),

      CluEError::NoPulseOpWithMultiplicity(p,spin_multiplicity) => write!(f,
          "pulse {} has no bath-spin-{} operator", 
          p,(*spin_multiplicity as f64 - 1.0)/2.0),

      CluEError::NoSpinOpForClusterSize(cluster_size,max_size) => write!(f,
          "no spin operators for clusters of size {} are built, \
only sizes up to {}",cluster_size,max_size),

      CluEError::NoSpinOpWithMultiplicity(spin_multiplicity) => write!(f,
          "no spin-{} operators are built", 
          (*spin_multiplicity as f64 - 1.0)/2.0),

      CluEError::NoStevensOp(k,q) => write!(f,
          "cannot get Stevens operator O_{}^{}",k,q),

      CluEError::NoStructureFile => write!(f,
          "no structure file defined"),

      CluEError::ParticlesClash(idx0,elmt0,idx1,elmt1,r,r_clash) => write!(f,
          "{} ({}) and {} ({}) are {} Å apart, closer than the \
clash distance of {} Å",idx0,elmt0,idx1,elmt1,r,r_clash),

      CluEError::Python(msg) => write!(f,
          "{}",msg),

      CluEError::ParticleIsNotActive(ref_index) => write!(f,
          "particle \"{}\" is not active", ref_index),

      CluEError::PropagatorAndDensityNotSameDimension(h,rho) => write!(f,
          "the propagator and density matrix are \"{}\" and \"{}\" dimensional \
respectively: they must have the same dimension",
          h,rho),


      CluEError::PartitionIsIncomplete => write!(f,
          "the partition must assign every spin to a block, with blocks \
numbered 0, 1, ..., n-1 without gaps; block numbers in \
cluster_expansion.clustering.partition_table must be 1, 2, ..., k"),

      CluEError::ResidueNumberAboveMaxInPDB(filename,line,max) => write!(f,
          "\"{}\" has a residue sequence number that is not a decimal \
number, as written for numbers above {}, in line \"{}\"", 
          filename, max, line),

      CluEError::SecondaryFilterRequiresAnIndex(group) => write!(f,
          "secondary particle group, \"{}\", requires a particle index",
          group),

      CluEError::SpaceGroupNotP1(space_group) => write!(f,
          "the structure has space group \"{}\", but \
structure.replicate_unit_cell only applies lattice translations; expand the \
structure to space group P 1 or set structure.replicate_unit_cell = false",
          space_group),

      CluEError::SpinGridWeightsDoNotSumToOne(sum,tolerance) => write!(f,
          "the weights of the detected spin grid sum to {}, which differs \
from 1 by more than detected_spin.weight_tolerance = {}", sum, tolerance),

      CluEError::SpinAlreadyPartitioned(idx,blk) => write!(f,
          "spin {} already partitioned into block {}", idx,blk),

      CluEError::StructureLargerThanCell(axis,span) => write!(f,
          "the structure spans {:.2} unit cells along cell edge {}, \
so the unit cell in the structure file does not describe its periodicity; \
correct or remove the unit cell (CRYST1) record",span,axis),

      CluEError::TensorNotSet(index) => write!(f,
          "no tensor set for particle {}",index),

      CluEError::ThermalizationRequiresEnsembleCCE => write!(f,
          "cluster_expansion.thermalize_cluster_populations = true requires \
cluster_expansion.ensemble = true"),

      CluEError::TooManyAtomsInPDB(filename,max) => write!(f,
          "\"{}\" has more than {} atoms in the selected model, \
which the PDB format cannot number with unique serial numbers",
          filename, max),

      CluEError::TooManyResiduesInPDB(filename,chain_id,max) => write!(f,
          "chain \"{}\" of \"{}\" has more than {} residues in the \
selected model, which the PDB format cannot number with unique residue \
sequence numbers", chain_id, filename, max),

      CluEError::TOMLArrayContainsMultipleTypes => write!(f,
          "TOML array contains data of different types"),

      CluEError::TOMLArrayDoesNotSpecifyATensor => write!(f,
          "cannot interpret TOML array as a tensor"),

      CluEError::TOMLArrayDoesNotSpecifyAVector => write!(f,
          "cannot interpret TOML array as a vector"),

      CluEError::TOMLArrayIsEmpty => write!(f,
          "TOML array is empty"),

      CluEError::TOMLArrayIsNotAMatrix => write!(f,
          "TOML array is not a matrix"),

      CluEError::TOMLArrayIsNotAComplexNumber => write!(f,
          "TOML array is not a complex number"),

      CluEError::TomlPulseStepSpecifier(msg) => write!(f,
          "{}",msg),

      CluEError::UnassignedCosubstitutionGroup(index)=> write!(f,
          "particle {} cannot be assigned to a cosubstitution group",
          index),

      CluEError::UnavailableSpinOp(op_pos,n_ops) => write!(f,
          "spin operator index, {}, exceeds vector length of {}",
          op_pos,n_ops),

      CluEError::UnequalLengths(vec0,len0,vec1,len1) => write!(f,
          "unequal lengths, {} has length {}, but {} has length {}",
          vec0,len0,vec1,len1),

      CluEError::UnexpectedPulseStep(position,expected) => write!(f,
          "pulse step {} should be {}",position,expected),

      CluEError::UnknownPulse(pulse) => write!(f,
          "no pulse named \"{}\"",pulse),

      CluEError::UnrecognizedOption(option) => write!(f,
          "unrecognized option \"{}\"",option),

      CluEError::UnrecognizedVectorSpecifier(specifier) => write!(f,
          "unrecognized vector specifier \"{}\"",specifier),

      CluEError::UnrecognizedUnit(quantity,unit,options) => write!(f,
          "unrecognized unit of {} \"{}\"; the options are {}",
          quantity,unit,options),

      CluEError::VectorSpecifierDoesNotSpecifyUniqueVector(specifier) 
        => write!(f,
          "vector specifier does not specify a unique vector: {}",specifier),

      CluEError::WrongClusterSizeForAnalyticCCE(given_size) => write!(f,
          "analytic 2-CCE cannot work with clusters of size {}",given_size),

      CluEError::WrongDelayDimension(d_given) => write!(f,
          "delay{} was provided, but only 1 and 2 are available",d_given),

      CluEError::WrongNumberOfAxes(num_axes, expected_num) => write!(f,
          "{} axes were provided, but {} are expected",num_axes, expected_num),

      CluEError::WrongNumberOfCellEdges(n) => write!(f,
          "expected 3 unit cell edges, but got {}",n),

      CluEError::WrongTauDimension(delay_dim, tau_dim) => write!(f,
          "\"{}\" cannot be used with \"{}\"; use \"delay\" with \"tau\" \
or \"-tau\", and \"delay2\" with \"tau2\" or \"-tau2\"",
          if *delay_dim == 2 {"delay2"} else {"delay"},
          if *tau_dim == 2 {"tau2"} else {"tau"}),

      CluEError::WrongBatchSignalLength(filename,found,expected) => write!(f,
          "the signals in \"{}\" have {} time points, but the run has {}; \
the file was written by a run with a different time axis",
          filename, found, expected),

      CluEError::WrongNumberOfBatchSignals(filename,found,expected) => write!(f,
          "\"{}\" has {} signals, but the batch has {} clusters; the file was \
written for a different cluster set or run.cluster_batch_size",
          filename, found, expected),

      CluEError::ZeroClusterBatchSize => write!(f,
          "run.cluster_batch_size must be at least 1"),

      CluEError::ZeroMaxClusterSize => write!(f,
          "the maximum cluster size must be at least 1"),
    }
  }
}

/// This implementation lets a `CluEError` be used wherever a
/// `std::error::Error` is expected, for example as `Box<dyn Error>`.
/// `source()` is left as `None`: `InvalidPulseStep` already includes the
/// message of the error it wraps in its own message.
impl std::error::Error for CluEError{}

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


#[cfg(test)]
mod tests{
  use super::*;
  //----------------------------------------------------------------------------
  // This test checks the text of messages whose wording, spacing, or 
  // key names were corrected.
  #[test]
  fn test_error_messages(){
    let s = |x: &str| x.to_string();
    let cases = vec![
      (CluEError::NoSpinOpForClusterSize(3,2),
       "no spin operators for clusters of size 3 are built, only sizes up to 2"),
      (CluEError::PropagatorAndDensityNotSameDimension(4,2),
       "the propagator and density matrix are \"4\" and \"2\" dimensional \
respectively: they must have the same dimension"),
      (CluEError::ParticlesClash(s("bath particle 0"),s("H"),
          s("bath particle 1"),s("H"),0.5,1.0),
       // The message uses the Angstrom sign, U+212B.
       "bath particle 0 (H) and bath particle 1 (H) are 0.5 \u{212B} apart, \
closer than the clash distance of 1 \u{212B}"),
      (CluEError::DetectedSpinDoesNotHaveAnActiveIndex,
       "the detected spin (reference index 0) is always active and has no \
position among the active bath particles"),
      (CluEError::NoDetectedSpinMultiplicity,
       "the detected spin multiplicity is not set"),
      (CluEError::AllVectorsNotSameLength(s("a.csv")),
       "cannot write \"a.csv\": all data columns must have the same length"),
      (CluEError::AllSignalsNotSameLength(s("a.csv")),
       "for \"a.csv\", signals must all have the same length"),
      (CluEError::CannotFindParticleForRefIndex(7),
       "no active bath particle matches reference index \"7\""),
      (CluEError::PulseSequenceNotSupported(s("lce"),s("hahn")),
       "cluster method \"lce\" only supports pulse_sequence.sequence = \
\"hahn\""),
      (CluEError::NoApplyPBC,
       "structure.replicate_unit_cell is not set (true, false, or \
[n1, n2, n3])"),
      (CluEError::CentroidGroupHasDistanceSelection(s("ring")),
       "group \"ring\" is used to place the detected spin and cannot use \
include.distance, since distances are measured from the detected spin"),
      (CluEError::PartitionIsIncomplete,
       "the partition must assign every spin to a block, with blocks \
numbered 0, 1, ..., n-1 without gaps; block numbers in \
cluster_expansion.clustering.partition_table must be 1, 2, ..., k"),
      (CluEError::NoDetectedSpinIdentity, "detected_spin.identity is not set"),
      (CluEError::NoDetectedSpinPosition, "detected_spin.position is not set"),
      (CluEError::NoDetectedParticle,
       "the structure has no detected spin; the detected spin is added \
when the structure is built from the config"),
      (CluEError::ThermalizationRequiresEnsembleCCE,
       "cluster_expansion.thermalize_cluster_populations = true requires \
cluster_expansion.ensemble = true"),
      (CluEError::InvalidStevensCoefficientShape(1,5,2),
       "Stevens coefficient entry 1 is a 5×2 matrix, but each entry must \
have one column"),
      (CluEError::CannotReadTOML(s("expected `=`")),
       "cannot parse TOML: expected `=`"),
      (CluEError::Python(s("invalid python code")), "invalid python code"),
      (CluEError::CannotWriteCSV(s("a.csv"),s("permission denied")),
       "cannot write CSV file \"a.csv\": permission denied"),
      (CluEError::CannotWriteFile(s("tensors.txt"),
          s("No space left on device (os error 28)")),
       "cannot write to \"tensors.txt\": No space left on device \
(os error 28)"),
      (CluEError::CannotCreateDir(s("out/info"),s("File exists (os error 17)")),
       "cannot create directory \"out/info\": File exists (os error 17)"),
      (CluEError::CannotReadGrid(s("grid.csv"),
          s("the file contains no grid points")),
       "cannot read grid from \"grid.csv\": the file contains no grid points; \
the grid should be specified as a csv file with one column per dimension, \
followed by a column for the weights"),
      (CluEError::NoDetectedSpinTransition,
       "detected_spin.transition is not set"),
      (CluEError::NoEnsembleCCE, "cluster_expansion.ensemble is not set"),
      (CluEError::NoMeanFields, "cluster_expansion.mean_fields is not set"),
      (CluEError::ZeroClusterBatchSize, 
       "run.cluster_batch_size must be at least 1"),
      (CluEError::WrongNumberOfBatchSignals(s("b.csv"),3,2),
       "\"b.csv\" has 3 signals, but the batch has 2 clusters; the file was \
written for a different cluster set or run.cluster_batch_size"),
      (CluEError::WrongBatchSignalLength(s("b.csv"),5,8),
       "the signals in \"b.csv\" have 5 time points, but the run has 8; \
the file was written by a run with a different time axis"),
      (CluEError::NoOrientationGrid, "orientations.grid is not set"),
      (CluEError::NoTimeIncrements2, 
       "pulse_sequence.tau2_increments is not set"),
      (CluEError::MissingHeader(1,2,s("a.csv")),
       "in \"a.csv\", every entry must have a header, but there are 1 \
headers and 2 columns of data"),
      (CluEError::CannotParseCellType(s("all_cells")),
       "cannot parse \"all_cells\" as cells; the options are \"all\", \
\"primary\", and \"extra\""),
      (CluEError::DensityMatrixTraceVanishes(s("Tr[exp(-βH)]"),1e-15),
       "cannot normalize the thermal density matrix: |Tr[exp(-βH)]| = 1e-15 \
is below 1e-12; check the temperature and the energy scale of the \
Hamiltonian"),
      (CluEError::CannotCreateFile(s("out/a.csv"),s("permission denied")),
       "cannot create \"out/a.csv\": permission denied"),
      (CluEError::GroupHasNoSelection(s("methyl")),
       "group \"methyl\" has no include or exclude table, so it selects no \
atoms; include = {} selects all atoms (of the cells set by cells)"),
      (CluEError::MutuallyExclusiveKeys(s("exchange_coupling"),
          s("c3_tunnel_splitting")),
       "exchange_coupling and c3_tunnel_splitting are mutually exclusive"),
      (CluEError::UnexpectedPulseStep(1,s("a pi pulse")),
       "pulse step 1 should be a pi pulse"),
      (CluEError::UnknownPulse(s("pi/4")), "no pulse named \"pi/4\""),
      (CluEError::WrongNumberOfCellEdges(2),
       "expected 3 unit cell edges, but got 2"),
      (CluEError::CannotParseLineInFile(s("a.pdb"),5,s("CONECT x")),
       "cannot parse line 5 of \"a.pdb\": \"CONECT x\""),
    ];
    for (error, message) in cases{
      assert_eq!(error.to_string(), message, "{:?}", error);
    }

    // CluEError implements std::error::Error, with the same text as Display
    // and no source, also for the variant that wraps another error.
    let boxed: Box<dyn std::error::Error> = Box::new(CluEError::NoApplyPBC);
    assert_eq!(boxed.to_string(), CluEError::NoApplyPBC.to_string());
    assert!(boxed.source().is_none());
    let step = CluEError::InvalidPulseStep(2,Box::new(CluEError::NoApplyPBC));
    assert!(std::error::Error::source(&step).is_none());
    let as_box = || -> Result<(),Box<dyn std::error::Error>>{
      Err(CluEError::NoApplyPBC)?
    };
    assert_eq!(as_box().unwrap_err().to_string(),
        CluEError::NoApplyPBC.to_string());

    // Messages do not end with a period, since main.rs adds one.
    let header = CluEError::ClusterFileContainsNoHeader(s("c.txt")).to_string();
    assert!(!header.ends_with('.'));

    // Units list the options for the quantity.
    let err = crate::physical_constants::distance_unit_to_meters("in")
      .unwrap_err();
    assert_eq!(err.to_string(), 
        "unrecognized unit of distance \"in\"; the options are \
fm, pm, nm, μm, mm, cm, m, km, \u{212B}");

    // Lebedev grid sizes name the key.
    let err = crate::integration_grid::IntegrationGrid::lebedev(7)
      .unwrap_err();
    assert!(err.to_string().starts_with(
          "orientations.number = 7 is not a valid Lebedev grid size"));
    assert!(err.to_string().contains("194, 230,"));
    assert!(err.to_string().contains("4802, 5294,"));
  }
  //----------------------------------------------------------------------------
}
