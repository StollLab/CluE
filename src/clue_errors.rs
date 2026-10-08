use crate::config::toml_keys::*;
use serde::Serialize;
use std::fmt;

/// `CluEError` contains the possible errors.
// After adding a new error, please update fmt() below as well.
#[derive(PartialEq,Debug,Clone,Serialize)]
pub enum CluEError{
  AllSignalsNotSameLength(String),
  AllVectorsNotSameLength(String),
  BondsAreNotDefined,
  BothStevensAndIST,
  CannotAddPointToGrid(usize,usize),
  CannotConvertSerialToIndex(u32),
  CannotCreateDir(String),
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
  CannotOpenFile(String),
  CannotParseCellType(String),
  CannotParseClusterMethod(String),
  CannotParseElement(String),
  CannotParseLine(String),
  CannotParseIsotope(String),
  CannotParseOrientations(String),
  CannotParsePartitioningMethod(usize,String),
  CannotParsePulseSequence(String),
  CannotParseSecondaryParticleFilter(String),
  CannotParseSpinOp(String),
  CannotPruneClustersMisMatchedSizes,
  CannotReadGrid(String),
  CannotReadTOML(String),
  CannotReadTOMLFile(String),
  CannotSampleBinomialDistribution(usize,f64),
  CannotTakeTrace(String),
  CannotWriteFile(String),
  ClusterFileContainsNoHeader(String),
  ClusterHasNoSignal(String),
  ClusterLineFormatError(String),
  ClusterTOMLCannotRead(String),
  ClusterTOMLIncorrectNumberOfClusters,
  ClusterTOMLNoClusters,
  ClusterTOMLNoNumberClusters,
  PulseSequenceNotSupported(String,String),
  DetectedSpinDoesNotHaveAnActiveIndex,
  Error(String),
  ExpectedClusterSetWithNSizes(usize,usize),
  ExpectedTOMLArray(String),
  ExpectedTOMLBool(String),
  ExpectedTOMLFloat(String),
  ExpectedTOMLUInt(String),
  ExpectedTOMLString(String),
  ExpectedTOMLTable(String),
  FailedKMeans,
  FilterAlreadySet(String),
  FilterNeedsALabel,
  FiltersOverlap(String,String),
  Generic(String),
  IncorrectNumberOfAxes(usize,usize),
  InorrectNumberOfCellOffsets(usize,usize),
  InvalidAxes,
  InvalidClusterPartitionKey,
  InvalidDensityMatrix,
  InvalidDetectionFrame,
  InvalidDetectionOperator,
  InvalidIST(String),
  InvalidKMeansSize,
  InvalidPulse(String),
  InvalidSpinMultiplicity(usize),
  InvalidStevens(String),
  LenghMismatchTimepointsIncrements(usize,usize),
  MeanFieldECCENotImplemented,
  MismatchedGroupNames(String,String),
  MissingDMatrix(usize),
  MissingFieldInCSVFile(String,String),
  MissingFilter(String),
  MissingGroupName,
  MissingHeader(usize,usize,String),
  NANTensorBathDipoleDipole(usize,String,usize,String),
  NANTensorBathZeeman(usize,String),
  NANTensorDetectedZeeman,
  NANTensorExchangeCoupling(usize,String,usize,String),
  NANTensorHyperfine(usize,String),
  NANTensorQuadrupole(usize,String),
  NANTensorZerofield(usize,String),
  NeighborListAndPartitionTableAreNotCompatible,
  NoApplyPBC,
  NoCentralSpin,
  NoCentralSpinCoor,
  NoCentralSpinIdentity,
  NoCentralSpinTransition,
  NoClusterBatchSize,
  NoClusterMethod,
  NoClustersOfSize(usize),
  NoClusterSource,
  NoClusterDensityMatrixMethod,
  NoDetectionFrame,
  NoDetectedSpinDensityMatrix,
  NoDetectedSpinDetectionOperator,
  NoDetectedSpinIdentity,
  NoDetectedSpinMultiplicity,
  NoDetectedSpinNotSet,
  NoDetectedSpinTransition,
  NoEnsembleCCE,
  NoGMatrixSpecifier,
  NoInputFile,
  NoKMeansSize,
  NoMagneticField,
  NoMeanFields,
  NoModelIndex,
  NoMaxClusterSize,
  NoMaxISTRank,
  NoNumberSystemInstances,
  NoOrientationGrid,
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
  NoUnitOfClustering,
  NoUnitOfDistance,
  NoUnitOfEnergy,
  NoUnitOfMagneticField,
  NoUnitOfTime,
  ParticlesClash(usize,String,usize,String,f64,f64),
  ParticleIsNotActive(usize),
  PropagatorAndDensityNotSameDimension(usize,usize),
  PartitionIsIncomplette,
  SecondaryFilterRequiresAnIndex(String),
  SpinAlreadyPartitioned(usize,i32),
  SpinPropertiesNeedsALabel,
  SpinPropertiesNeedsAnIsotope(String),
  StructurePropertiesNeedsALabel,
  StructurePropertiesDoesNotNeedAnIsotope(String),
  TensorNotSet(usize),
  TooManyExchangeCouplingsSpecified(usize),
  TOMLArrayContainsMultipleTypes,
  TOMLArrayDoesNotSpecifyATensor,
  TOMLArrayDoesNotSpecifyAVector,
  TOMLArrayIsEmpty,
  TOMLArrayIsNotAMatrix,
  TOMLArrayIsNotAComplexNumber,
  TOMLValueDoesNotSpecifyAPartitionTable,
  ErrorPulseSequence(String),
  TomlPulseStepSpecifier(String),
  UnavailableSpinOp(usize,usize),
  UnassignedCosubstitutionGroup(usize),
  UnequalLengths(String,usize,String,usize),
  UnrecognizedOption(String),
  UnrecognizedVectorSpecifier(String),
  UnrecognizedUnit(String),
  VectorSpecifierDoesNotSpecifyUniqueVector(String),
  WrongClusterSizeForAnalyticCCE(usize),
  WrongDelayDimension(usize),
  WrongNumberOfAxes(usize,usize),
  WrongTauDimension(usize,usize),
}

impl fmt::Display for CluEError{
  // This function defines the error messages.
  fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
    match self{

      CluEError::AllSignalsNotSameLength(filename) => write!(f,
          "for \"{}\",signals must all have the same length", filename),

      CluEError::AllVectorsNotSameLength(filename) => write!(f,
          "for \"{}\",signals must all have the same length", filename),

      CluEError:: BondsAreNotDefined => write!(f,
          "no chemical bonds are established"),

      CluEError:: BothStevensAndIST => write!(f,
          "{} and {} are mutually exclusive",KEY_IST_COEF,KEY_STEVENS_COEF),

      CluEError::CannotAddPointToGrid(point_dim, grid_dim) => write!(f,
          "cannot add {}D point t0 {}D grid",point_dim, grid_dim),

      CluEError::CannotConvertSerialToIndex(serial) => write!(f,
          "cannot convert serial id, {}, to an index",serial),

      CluEError::CannotDiagonalizeHamiltonian(matrix) => write!(f,
          "cannot diagonalize Hamiltonian,\n{}",matrix),

      CluEError::CannotDiagonalizeOperator(matrix) => write!(f,
          "cannot diagonalize \n{}",matrix),

      CluEError::CannotCreateDir(path) => write!(f,
          "cannot create directory \"{}\"",path),

      CluEError::CannotExpandBlockClusters => write!(f,
          "cannot expand block clusters"),

      CluEError::CannotFindCellID(idx) => write!(f,
          "cannot determine cell id for particle {}",idx),

      CluEError::CannotFindSpinOp(sop) => write!(f,
          "cannot find spin operator \"{}\"",sop),

      CluEError::CannotFindParticleForRefIndex(ref_index) => write!(f,
          "cannot find bath index for reference index \"{}\"", ref_index),

      CluEError::CannotFindBathIndexFromRefIndex(ref_index) => write!(f,
          "cannot find bath index for reference index \"{}\"", ref_index),

      CluEError::CannotFindRefIndexFromBathIndex(bath_index) => write!(f,
          "cannot find reference index for bath index \"{}\"", bath_index),

      CluEError::CannotFindRefIndexFromNthActive(n) => write!(f,
          "cannot find reference index for the {}th active particle", n),

      CluEError::CannotNormalizeVector => write!(f,
          "cannot normalize vector"),

      CluEError::CannotOpenFile(file) => write!(f,
          "cannot open \"{}\"", file),

      CluEError::CannotParseCellType(cell_type) => write!(f,
          "cannot parse \"{}\" a a cell type", cell_type),

      CluEError::CannotParseClusterMethod(method) => write!(f,
          "cannot parse \"{}\" as a cluster method", method),

      CluEError::CannotParseElement(element) => write!(f,
          "cannot parse \"{}\" as an element", element),

      CluEError::CannotParseOrientations(ori) => write!(f,
          "cannot parse \"{}\" as an orientation averaging scheme", ori),

      CluEError::CannotParseLine(line) => write!(f,
          "cannot parse line \"{}\"", line),

      CluEError::CannotParseIsotope(isotope) => write!(f,
          "cannot parse \"{}\" as an isotope", isotope),

      CluEError::CannotParsePartitioningMethod(line_num, part_method) 
          => if *line_num > 0 {
            write!(f, "line {}: cannot parse partitioning method \"{}\"", 
            line_num, part_method)
          }else{
            write!(f, "cannot parse partitioning method \"{}\"", 
            part_method)
          },

      CluEError::CannotParsePulseSequence(seq) => write!(f,
          "cannot parse \"{}\" as a pulse sequence", seq),

      CluEError::CannotParseSecondaryParticleFilter(group) => write!(f,
          "cannot parse secondary particle group \"{}\"", group),

      CluEError::CannotParseSpinOp(s) => write!(f,
          "cannot parse \"{}\" as a spin operator", s),

      CluEError::CannotSampleBinomialDistribution(n,p) => write!(f,
          "cannot sample from the binomial distribution B(n={},p={})",
          n,p),

      CluEError::CannotTakeTrace(matrix) => write!(f,
          "cannot take trace of \"{}\"",matrix),

      CluEError::CannotPruneClustersMisMatchedSizes => write!(f,
          "cannot prune clusters because clusters_to_keep has the wrong size"),

      CluEError::CannotReadGrid(filename) => write!(f,
          "cannot read grid from \"{}\": \
the grid should be specified as a csv file with one column per dimension, \
followed by a column for the weights", filename),

      CluEError::CannotReadTOML(string) => write!(f,
          "cannot parse \"{}\" as TOML",string),

      CluEError::CannotReadTOMLFile(msg) => write!(f,
          "cannot parse TOML file \"{}\"",msg),

      CluEError::ClusterFileContainsNoHeader(file) => write!(f,
          "Cluster file \"{}\" does not contain the correct header: \
\"#[clusters, number_clusters = [N1,N2,...Nn] ]\", where Ni is the number \
of clusters of size i in the file, and the list runs from to clusters of \
size n.", file),

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

      CluEError::CannotWriteFile(file) => write!(f,
          "cannot write to \"{}\"", file),

      CluEError::PulseSequenceNotSupported(fun,seq) => write!(f,
        "\"{}\" does not support {}",fun,seq),

      CluEError::Error(error_message) => write!(f,
          "{}", error_message),

      CluEError::DetectedSpinDoesNotHaveAnActiveIndex => write!(f,
          "the detected spin is always active and so does have an index \
fo nth active"),

      CluEError::ExpectedClusterSetWithNSizes(n_exp,n_act) => write!(f,
          "expected a cluster set with {} sizes, but got {} sizes",
          n_exp, n_act),

      CluEError::ExpectedTOMLArray(type_str) => write!(f,
          "expected array, but got {}",type_str),

      CluEError::ExpectedTOMLBool(type_str) => write!(f,
          "expected bool, but got {}",type_str),

      CluEError::ExpectedTOMLFloat(type_str) => write!(f,
          "expected float, but got {}",type_str),

      CluEError::ExpectedTOMLUInt(type_str) => write!(f,
          "expected int > 0, but got {}",type_str),

      CluEError::ExpectedTOMLString(type_str) => write!(f,
          "expected string, but got {}",type_str),

      CluEError::ExpectedTOMLTable(type_str) => write!(f,
          "expected table, but got {}",type_str),

      CluEError::FailedKMeans => write!(f,
          "k-means failed"),

      CluEError::FilterAlreadySet(filter) => write!(f,
          "{} has already been set",filter),

      CluEError::FilterNeedsALabel => write!(f,
          "group requires a label to be set: #[group(label = LABEL)]"),

      CluEError::FiltersOverlap(label0,label1) => write!(f,
          "groups \"{}\" and \"{}\" overlap: \
particles must not reside in more than one group",label0,label1),

      CluEError::Generic(s) => write!(f,"{}",s),

      CluEError::IncorrectNumberOfAxes(n,n_ref)=> write!(f,
          "expected {} axes, but {} were provided",n_ref, n),

      CluEError::InorrectNumberOfCellOffsets(n,n_ref)=> write!(f,
          "expected {} cell offsets, but {} were provided",n_ref, n),

      CluEError::InvalidAxes => write!(f,
          "invalid axes"),

      CluEError::InvalidClusterPartitionKey => write!(f,
          "invalid cluster partition key"),

      CluEError::InvalidDensityMatrix => write!(f,
          "invalid density matrix"),

      CluEError::InvalidDetectionFrame => write!(f,
          "invalid detection frame"),

      CluEError::InvalidDetectionOperator => write!(f,
          "invalid detection operator"),

      CluEError::InvalidIST(ist) => write!(f,
          "invalid {} \"{}\"",KEY_IST_COEF,ist),

      CluEError::InvalidKMeansSize => write!(f,
          "kmeans_size must be at least 1"),

      CluEError::InvalidPulse(pulse) => write!(f,
          "invalid {}",pulse),

      CluEError::InvalidSpinMultiplicity(n) => write!(f,
          "invalid spin multiplicity \"{}\"",n),

      CluEError::InvalidStevens(ist) => write!(f,
          "invalid {} \"{}\"",KEY_STEVENS_COEF,ist),

      CluEError::LenghMismatchTimepointsIncrements(n_dts,dts) => write!(f,
          "there are {} timepoint specifications, but {} time increments",
          n_dts,dts),

      CluEError::MeanFieldECCENotImplemented => write!(f,
          "mean field averaging is not implemented for ensemble CCE"),

      CluEError::MismatchedGroupNames(name0,name1) => write!(f,
          "group names \"{}\" and \"{}\" do not match",name0,name1),

      CluEError::MissingDMatrix(l) => write!(f,
          "no D matrix for l = {}",l),

      CluEError::MissingFilter(label) => write!(f,
          "no group with label \"{}\"",label),

      CluEError::MissingFieldInCSVFile(field,file_name) => write!(f,
          "CSV file \"{}\" has no column \"{}\"",file_name,field),

      CluEError::MissingGroupName => write!(f,
          "no group name specified"),

      CluEError::MissingHeader(n_headers, n_cols,filename) => write!(f,
          "in \"{}\", every entry must have a header, \
but there are {} headers and {} columns of data.",filename,n_headers,n_cols),

      CluEError::NANTensorBathDipoleDipole(idx0,isotope0,idx1,isotope1) 
        => write!(f,
          "dipole-dipole tensor for particles {} {} and {} {} contains NANs",
          idx0,isotope0,idx1,isotope1),

      CluEError::NANTensorBathZeeman(idx,isotope) => write!(f,
          "zeeman tensor for particle {} {} contains NANs",idx,isotope),

      CluEError::NANTensorDetectedZeeman => write!(f,
          "Zeeman tensor for the detected particle contains NANs"),

      CluEError::NANTensorExchangeCoupling(idx0,isotope0,idx1,isotope1) 
        => write!(f,
          "exchange tensor for particles {} {} and {} {} contains NANs",
          idx0,isotope0,idx1,isotope1),

      CluEError::NANTensorHyperfine(idx,isotope) => write!(f,
          "hyperfine tensor for particle {} {} contains NANs",idx,isotope),

      CluEError::NANTensorQuadrupole(idx,isotope) => write!(f,
          "electric quadrupole tensor for particle {} {} contains NANs",
          idx,isotope),

      CluEError::NANTensorZerofield(idx,isotope) => write!(f,
          "zerofield tensor for particle {} {} contains NANs",
          idx,isotope),

      CluEError::NeighborListAndPartitionTableAreNotCompatible  => write!(f,
          "neighbor list and partition table are incompatible"),

      CluEError::NoApplyPBC => write!(f,
          "please set \"do_replicate_unit_cell: bool\" to select if \
periodic boundary conditions should be applied"),

      CluEError::NoCentralSpin => write!(f,
          "no detected spin defined"),

      CluEError::NoCentralSpinCoor => write!(f,
          "coordinates for the detected spin were not defined"),

      CluEError::NoCentralSpinIdentity => write!(f,
          "the detected spin's identity was not defined"),

      CluEError::NoCentralSpinTransition => write!(f,
          "the detected spin transition was not defined"),

      CluEError::NoClusterBatchSize => write!(f,
          "batch size for clusters not specified"),

      CluEError::NoClusterMethod => write!(f,
          "no cluster method specified"),

      CluEError::NoClustersOfSize(size) => write!(f,
          "cannot find any clusters of size {}", size),

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
          "detected_spin_identity is not set"),

      CluEError::NoDetectedSpinMultiplicity => write!(f,
          "detected_spin_identity is not set"),

      CluEError::NoDetectedSpinNotSet => write!(f,
          "detected_spin is not set"),

      CluEError::NoDetectedSpinTransition => write!(f,
          "detected_transition is not set"),

      CluEError::NoEnsembleCCE => write!(f,
          "ensemble_cce: bool is not set"),

      CluEError::NoGMatrixSpecifier => write!(f,
          "no g-matrix specifier"),

      CluEError::NoInputFile => write!(f,
          "no input file"),

      CluEError::NoKMeansSize => write!(f,
          "k-means requires kmeans_size to be set"),

      CluEError::NoMagneticField => write!(f,
          "please specify the applied magnetic field"),

      CluEError::NoMeanFields => write!(f,
          "please whether or not mean fields should be used"),

      CluEError::NoModelIndex => write!(f,
          "PDB model not selected"),

      CluEError::NoMaxClusterSize => write!(f,
          "maximum cluster size not set"),

      CluEError::NoMaxISTRank => write!(f,
          "maximum irreducible spherical tensor rank not set"),

      CluEError::NoNumberSystemInstances => write!(f,
          "number_runs is not set"),

      CluEError::NoOrientationGrid => write!(f,
          "orientation_grid is not defined"),

      CluEError::NoPartitioningMethod => write!(f,
          "no partitioning method is defined"),

      CluEError::NoPulseSequence => write!(f,
          "no pulse sequence is defined"),

      CluEError::NoRadius => write!(f,
          "system radius not set"),

      CluEError::NotA3DRotationMatrix(matrix) => write!(f,
          "matrix \"\n{}\n\"does not correspond to a valid rotation matrix", 
          matrix),

      CluEError::NotA3DVector(dim) => write!(f,
          "vector is {}-dimensional, not 3-dimensional", 
          dim),

      CluEError::NotALebedevGrid(n_ori) => write!(f,
          "\"{}\" is not a valid Lebedev grid; please choose n_ori in \
{{6, 14, 26, 38, 50, 74, 86, 110, 146, 170, \
194,  230, 266, 302, 350, 434, 590, 770, 974, 1202, \
1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802,5294, 5810}}",n_ori),

      CluEError::NotAProperSubset(subcluster, cluster) => write!(f,
          "{} is not a proper subset of {}", subcluster, cluster),

      CluEError::NoTemperature => write!(f,
          "no temperature specified"),

      CluEError::NoTensorSpecifier => write!(f,
          "no specifier for tensor"),

      CluEError::NoTensorValues => write!(f,
          "no values specified for tensor"),

      CluEError::NoTimeAxis => write!(f,
          "the time-axis has not been built"),

      CluEError::NoTimeIncrements => write!(f,
          "no tau increments defined"),

      CluEError::NoTimeIncrements2 => write!(f,
          "no taa2 increments defined"),

      CluEError::NoTimepoints => write!(f,
          "please specify how many timepoints there are for each tau increment"),

      CluEError::NoTimepoints2 => write!(f,
          "please specify how many timepoints there are for each tau2 increment"),

      CluEError::NoUnitOfClustering => write!(f,
          "please specify a unit of clustering"),

      CluEError::NoUnitOfDistance => write!(f,
          "please specify a unit of distance"),

      CluEError::NoUnitOfEnergy => write!(f,
          "please specify a unit of energy"),

      CluEError::NoUnitOfMagneticField => write!(f,
          "please specify a unit of magnetic field"),

      CluEError::NoUnitOfTime => write!(f,
          "please specify a unit of time"),

      CluEError::NoRunInParallel => write!(f,
          "run_in_parallel is not set"),

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
          "no spin operators for clusters of size {} are built,\
          only sizes up to {}",cluster_size,max_size),

      CluEError::NoSpinOpWithMultiplicity(spin_multiplicity) => write!(f,
          "no spin-{} operators are built", 
          (*spin_multiplicity as f64 - 1.0)/2.0),

      CluEError::NoStevensOp(k,q) => write!(f,
          "Cannot get Stevens operator O_{}^{}",k,q),

      CluEError::NoStructureFile => write!(f,
          "no structure file defined"),

      CluEError::ParticlesClash(idx0,elmt0,idx1,elmt1,r,r_clash) => write!(f,
          "particle {} {} and particle {} {} are {} Å apart, closer than the
clash distance of {} Å",idx0,elmt0,idx1,elmt1,r,r_clash),

      CluEError::ParticleIsNotActive(ref_index) => write!(f,
          "particle \"{}\" is not active", ref_index),

      CluEError::PropagatorAndDensityNotSameDimension(h,rho) => write!(f,
          "the propagator and density matrix are \"{}\" and \"{}\" dimensional\
respectively: they must have the same dimension",
          h,rho),


      CluEError::PartitionIsIncomplette => write!(f,
          "partition table should contain n cells indexed by [0,n-1]"),

      CluEError::SecondaryFilterRequiresAnIndex(group) => write!(f,
          "secondary particle group, \"{}\", requires a particle index",
          group),

      CluEError::SpinAlreadyPartitioned(idx,blk) => write!(f,
          "spin {} already partitioned into block {}", idx,blk),

      CluEError::SpinPropertiesNeedsALabel => write!(f,
          "spin_properties requires a label to be set: 
#[spin_properties(label = LABEL, isotope = ISOTOPE)]"),

      CluEError::SpinPropertiesNeedsAnIsotope(label) => write!(f,
          "spin_properties requires an isotope to be set: 
#[spin_properties(label = {}, isotope = ISOTOPE)]",label),

      CluEError::StructurePropertiesNeedsALabel => write!(f,
          "structure_properties requires a label to be set: 
#[structure_properties(label = LABEL)]"),

      CluEError::StructurePropertiesDoesNotNeedAnIsotope(label) => write!(f,
          "structure_properties does not require an isotope to be set: 
#[spin_properties(label = {})]",label),

      CluEError::TensorNotSet(index) => write!(f,
          "no tensor set for particle {}",index),

      CluEError::TooManyExchangeCouplingsSpecified(n) => write!(f,
          "expected 1 exchange coupling, but got {}", 
          n),

      CluEError::TOMLArrayContainsMultipleTypes => write!(f,
          "TOML array contains data of different types", 
          ),

      CluEError::TOMLArrayDoesNotSpecifyATensor => write!(f,
          "cannot interpret TOML array as a tensor", 
          ),

      CluEError::TOMLArrayDoesNotSpecifyAVector => write!(f,
          "cannot interpret TOML array as a vector", 
          ),

      CluEError::TOMLArrayIsEmpty => write!(f,
          "TOML array is empty", 
          ),

      CluEError::TOMLArrayIsNotAMatrix => write!(f,
          "TOML array is not a matrix", 
          ),

      CluEError::TOMLArrayIsNotAComplexNumber => write!(f,
          "TOML array is not a complex number", 
          ),

      CluEError::TOMLValueDoesNotSpecifyAPartitionTable => write!(f,
          "TOML value is not a partition table", 
          ),

      CluEError::ErrorPulseSequence(msg) => write!(f,
          "{}",msg),

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

      CluEError::UnrecognizedOption(option) => write!(f,
          "unrecognized option \"{}\"",option),

      CluEError::UnrecognizedVectorSpecifier(specifier) => write!(f,
          "unrecognized vector specifier \"{}\"",specifier),

      CluEError::UnrecognizedUnit(unit) => write!(f,
          "unrecognized unit \"{}\"",unit),

      CluEError::VectorSpecifierDoesNotSpecifyUniqueVector(specifier) 
        => write!(f,
          "vector specifier does not specify unique vector \"{}\"",specifier),

      CluEError::WrongClusterSizeForAnalyticCCE(given_size) => write!(f,
          "analytic 2-CCE cannot work with clusters of size {}",given_size),

      CluEError::WrongDelayDimension(d_given) => write!(f,
          "delay{} was provided, but only 1 and 2 are available",d_given),

      CluEError::WrongNumberOfAxes(num_axes, expected_num) => write!(f,
          "{} axes were provided, but {} are expected",num_axes, expected_num),

      CluEError::WrongTauDimension(d_given, d_expected) => write!(f,
          "tau{} was provided, but tau{} is expected",d_given, d_expected),
    }
  }
}

