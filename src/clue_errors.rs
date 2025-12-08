use serde::Serialize;
use std::fmt;

/// `CluEError` contains the possible errors.
// After adding a new error, please update fmt() below as well.
#[derive(PartialEq,Debug,Clone,Serialize)]
pub enum CluEError{
  AllSignalsNotSameLength(String),
  AtomDoesNotSpecifyElement(usize),
  AllVectorsNotSameLength(String),
  BondsAreNotDefined,
  CannotAddPointToGrid(usize,usize),
  CannotAddTokens,
  CannontAugmentFilter(usize,String),
  CannotCombineTokens(usize),
  CannotConvertToFloat(usize,String),
  CannotConvertSerialToIndex(u32),
  CannotConvertToVector(usize),
  CannotCreateDir(String),
  CannotDiagonalizeHamiltonian(String),
  CannotDivTokens,
  CannotExpandBlockClusters,
  CannotFindCellID(usize),
  CannotFindSpinOp(String),
  CannotFindParticleForRefIndex(usize),
  CannotFindBathIndexFromRefIndex(usize),
  CannotFindRefIndexFromBathIndex(usize),
  CannotFindRefIndexFromNthActive(usize),
  CannotInferEigenvalues(usize),
  CannotMatchVertexToIndex(usize),
  CannotMulTokens,
  CannotOpenFile(String),
  CannotParseCellType(String),
  CannotParseClusterMethod(String),
  CannotParseDensityMatrix(String),
  CannotParseElement(String),
  CannotParseLine(String),
  CannotParseIsotope(String),
  CannotParseOrientations(String),
  CannotParsePartitioningMethod(usize,String),
  CannotParsePulseSequence(String),
  CannotParseRHS(usize),
  CannotParseSecondaryParticleFilter(String),
  CannotParseUnitOfClustering(usize,String),
  CannotPowTokens,
  CannotPruneClustersMisMatchedSizes,
  CannotReadGrid(String),
  CannotReadTOML(String),
  CannotReadTOMLFile(String),
  CannotSampleBinomialDistribution(usize,f64),
  CannotSetExchangeCoupling(usize),
  CannotSubTokens,
  CannotTakeTrace(String),
  CannotWriteFile(String),
  CIFCannotDeclareDataItem(usize,String),
  CIFCannotFindSymmetryEquivalentSite(String),
  CIFCannotParseFloat(String),
  CIFCannotStartLoop(usize),
  CIFDataItemAlreadyDeclared(usize,String),
  CIFInconsistantCoorDataSizes(usize),
  CIFNoDeclaredDataItem(usize,String),
  CIFNoField(String),
  CIFNoFractX,
  CIFNoFractY,
  CIFNoFractZ,
  CIFNoTypeSymbol,
  ClusterFileContainsNoHeader(String),
  ClusterHasNoSignal(String),
  ClusterLineFormatError(String),
  ClusterTOMLCannotRead(String),
  ClusterTOMLIncorrectNumberOfClusters,
  ClusterTOMLNoClusters,
  ClusterTOMLNoNumberClusters,
  ClusterTOMLWrongClusterSize,  
  ConfigModeNotRecognized(String),
  DeprecatedKeywordReplaced(usize,String,String),
  DetectedSpinDoesNotHaveAnActiveIndex,
  EmptyVector(usize),
  Error(String),
  ExpectedClusterSetWithNSizes(usize,usize),
  ExpectedEquality(usize),
  ExpectedBoolRHS(usize),
  ExpectedFloatRHS(usize),
  ExpectedIntRHS(usize),
  ExpectedNonNegativeIntRHS(usize),
  ExpectedVecOfNFloatsRHS(usize,usize),
  ExpectedNumber(usize),
  ExpectedTOMLArray(String),
  ExpectedTOMLBool(String),
  ExpectedTOMLFloat(String),
  ExpectedTOMLString(String),
  ExpectedTOMLTable(String),
  FilterNoMaxDistance(String),
  FilterNeedsALabel,
  FilterNoMinDistance(String),
  FiltersOverlap(String,String),
  InconsistentExhangeGroupActiveStatus(usize),
  IncorrectFormattingIsotopeAbundances(usize),
  IncorrectNumberOfAxes(usize,usize),
  InorrectNumberOfCellOffsets(usize,usize),
  IndexOutOfBounds(usize,usize,usize),
  InvalidArgument(usize,String),
  InvalidAxes,
  InvalidClusterPartitionKey,
  InvalidConfigFile(String),
  InvalidPulseSequence(usize),
  InvalidSecondaryFilter(usize,String),
  InvalidToken(usize,String),
  IsotopeAbundancesCannotBeNormalized(usize),
  IsotopeAbundancesMustBeNonnegative(usize),
  LenghMismatchTimepointsIncrements(usize,usize),
  MismatchedGroupNames(String,String),
  MissingFieldInCSVFile(String,String),
  MissingFilter(String),
  MissingFilterArgument(usize,String),
  MissingFilterLabel(usize),
  MissingGroupName,
  MissingHeader(usize,usize,String),
  MissingProperties(String),
  MissingPropertiesLabel(usize),
  ModeAttributeWrongBrackets,
  ModeAttributeWrongOption(String),
  ModeAttributeWrongSharp,
  MultipleCosubstitutionGroups(usize),
  NANTensorBathDipoleDipole(usize,String,usize,String),
  NANTensorBathZeeman(usize,String),
  NANTensorDetectedZeeman,
  NANTensorExchangeCoupling(usize,String,usize,String),
  NANTensorHyperfine(usize,String),
  NANTensorQuadrupole(usize,String),
  NeighborListAndPartitionTableAreNotCompatible,
  NoApplyPBC,
  NoArgument(usize),
  NoBathGMatrixSpecifier(String,String),
  NoCentralSpin,
  NoCentralSpinCoor,
  NoCentralSpinIdentity,
  NoCentralSpinTransition,
  NoClashDistancePBC,
  NoClusterBatchSize,
  NoClusterMethod,
  NoClustersOfSize(usize),
  NoClusterSource,
  NoDensityMatrixMethod,
  NoDetectedSpinIdentity,
  NoDetectedSpinMultiplicity,
  NoDetectedSpinNotSet,
  NoDetectedSpinTransition,
  NoExtracellIsotopicDistribution(String),
  NoGMatrixSpecifier,
  NoGMatrixValues,
  NoHyperfineSpecifier(String,String),
  NoInputFile,
  NoMagneticField,
  NoModelIndex,
  NoMaxClusterSize,
  NoNeighborCutoffDistance,
  NoNumberSystemInstances,
  NoOrientationGrid,
  NoPartitioningMethod,
  NoPulseSequence,
  NoQuadrupoleSpecifier(String,String),
  NoRadius,
  NoRelationalOperators(usize),
  NoRHS(usize),
  NoSpinOpForClusterSize(usize,usize),
  NoSpinOpWithMultiplicity(usize),
  NoStructureFile,
  NotA3DVector(usize),
  NotALebedevGrid(usize),
  NotAnOperator(usize,String),
  NotAProperSubset(String,String),
  NoTemperature,
  NoTensorSpecifier,
  NoTensorValues,
  NoTimeAxis,
  NoTimeIncrements,
  NoTimepoints,
  NoUnitOfClustering,
  NoUnitOfDistance,
  NoUnitOfEnergy,
  NoUnitOfMagneticField,
  NoUnitOfTime,
  OptionAlreadySet(usize,String),
  ParticlesClash(usize,String,usize,String,f64,f64),
  ParticleIsNotActive(usize),
  PartitionIsIncomplette,
  SaveNameEmpty,
  SaveNameNotSet,
  SecondaryFilterRequiresAnIndex(String),
  SpinAlreadyPartitioned(usize,i32),
  SpinPropertiesNeedsALabel,
  SpinPropertiesNeedsAnIsotope(String),
  StructurePropertiesNeedsALabel,
  StructurePropertiesDoesNotNeedAnIsotope(String),
  TensorNotSet(usize),
  TooFewRHSArguments(usize),
  TooManyRelationalOperators(usize),
  TooManyExchangeCouplingsSpecified(usize),
  TooManyRHSArguments(usize),
  TOMLArrayContainsMultipleTypes,
  TOMLArrayDoesNotSpecifyATensor,
  TOMLArrayDoesNotSpecifyAVector,
  TOMLArrayIsEmpty,
  TOMLValueDoesNotSpecifyAPartitionTable,
  TOMLValueIsNotABool,
  UnavailableSpinOp(usize,usize),
  UnassignedCosubstitutionGroup(usize),
  UnequalLengths(String,usize,String,usize),
  UnmatchedBlockComment(usize),
  UnmatchedDelimiter(usize),
  UnrecognizedOption(String),
  UnrecognizedVectorSpecifier(String),
  UnrecognizedUnit(String),
  VectorSpecifierDoesNotSpecifyUniqueVector(String),
  WrongClusterSizeForAnalyticCCE(usize),
  WrongNumberOfAxes(usize,usize),
  WrongOrientationGridDim(usize,usize,usize),
  WrongProbabilityDistributionDim(usize,usize,usize),
  WrongVectorLength(usize,usize,usize)
}

impl fmt::Display for CluEError{
  // This function defines the error messages.
  fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
    match self{

      CluEError::AllSignalsNotSameLength(filename) => write!(f,
          "for \"{}\",signals must all have the same length", filename),

      CluEError::AtomDoesNotSpecifyElement(serial) => write!(f,
          "atom {} does not specify an element",serial),

      CluEError::AllVectorsNotSameLength(filename) => write!(f,
          "for \"{}\",signals must all have the same length", filename),

      CluEError:: BondsAreNotDefined => write!(f,
          "no chemical bonds are established"),

      CluEError::CannotAddPointToGrid(point_dim, grid_dim) => write!(f,
          "cannot add {}D point t0 {}D grid",point_dim, grid_dim),

      CluEError::CannotAddTokens => write!(f,
          "cannot add tokens meaningfully"),

      CluEError::CannontAugmentFilter(index,secondary_filter) => write!(f,
          "cannot use secondary group \"{}\" with particle {}",
         secondary_filter,index),

      CluEError::CannotCombineTokens(line_number) => write!(f,
          "line {}, cannot combine tokens meaningfully", line_number),

      CluEError::CannotConvertSerialToIndex(serial) => write!(f,
          "cannot convert serial id, {}, to an index",serial),

      CluEError::CannotDiagonalizeHamiltonian(matrix) => write!(f,
          "cannot diagonalize Hamiltonian,\n {}",matrix),

      CluEError::CannotCreateDir(path) => write!(f,
          "cannot create directory \"{}\"",path),

      CluEError::CannotDivTokens => write!(f,
          "cannot divide tokens meaningfully"),

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

      CluEError::CannotInferEigenvalues(line_number) => write!(f,
          "line {}, cannot infer eigenvalues from input",line_number),
      
      CluEError::CannotMatchVertexToIndex(vertex) => write!(f,
          "cannot match vertex {} to an index",vertex),
      
      CluEError::CannotMulTokens => write!(f,
          "cannot multiply tokens meaningfully"),

      CluEError::CannotOpenFile(file) => write!(f,
          "cannot open \"{}\"", file),

      CluEError::CannotParseCellType(cell_type) => write!(f,
          "cannot parse \"{}\" a a cell type", cell_type),

      CluEError::CannotParseClusterMethod(method) => write!(f,
          "cannot parse \"{}\" as a cluster method", method),

      CluEError::CannotParseDensityMatrix(d) => write!(f,
          "cannot parse \"{}\" as a density matrix method", d),

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

      CluEError::CannotParseRHS(line_number) => write!(f,
          "line {}, cannot parse line right hand side", line_number),

      CluEError::CannotParseSecondaryParticleFilter(group) => write!(f,
          "cannot parse secondary particle group \"{}\"", group),

      CluEError::CannotParseUnitOfClustering(line_num, unit) 
          => write!(f, "line {}: cannot parse unit of clustering \"{}\"", 
          line_num, unit),

      CluEError::CannotSampleBinomialDistribution(n,p) => write!(f,
          "cannot sample from the binomial distribution B(n={},p={})",
          n,p),

      CluEError::CannotSetExchangeCoupling(line_number) => write!(f,
          "line {}, cannot set exchange coupling",line_number),

      CluEError::CannotSubTokens => write!(f,
          "cannot subtract tokens meaningfully"),

      CluEError::CannotTakeTrace(matrix) => write!(f,
          "cannot take trace of \"{}\"",matrix),

      CluEError::CannotPowTokens => write!(f,
          "cannot do token^token meaningfully"),

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

      CluEError::CannotConvertToFloat(line_number, token) => write!(f,
          "line {}, cannot convert \"{}\" to type float", line_number,token),

      CluEError::CannotConvertToVector(line_number) => write!(f,
          "line {}, cannot find vector", line_number),

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

      CluEError:: ClusterTOMLWrongClusterSize=> write!(f,
          "in clusters.toml, there is at least one cluster of the wrong size"),

      CluEError::CannotWriteFile(file) => write!(f,
          "cannot write to \"{}\"", file),

      CluEError::CIFCannotDeclareDataItem(line_number,item) => write!(f,
          "line {}, cannot declare data item \"{}\"", line_number,item),

      CluEError::CIFCannotFindSymmetryEquivalentSite(item) => write!(f,
          "cannot find symmetry equivalent site from \"{}\"", item),

      CluEError::CIFCannotParseFloat(item) => write!(f,
          "cannot parse \"{}\" as a float", item),

      CluEError::CIFCannotStartLoop(line_number) => write!(f,
          "line {}, cannot start loop_", line_number),

      CluEError::CIFDataItemAlreadyDeclared(line_number,item) => write!(f,
          "line {}, cannot declare data item again\"{}\"", line_number,item),

      CluEError::CIFNoDeclaredDataItem(line_number,item) => write!(f,
          "line {}, no declared data item for \"{}\"", line_number,item),

      CluEError::CIFInconsistantCoorDataSizes(n_atoms)=> write!(f,
          "x, y and z coordinates must have {} elements in cif file",
          n_atoms),

      CluEError::CIFNoField(field) => write!(f,
          "no \"{}\" in cif file",field),

      CluEError::CIFNoFractX => write!(f,
          "no \"_atom_site_fract_x\" in cif file"),

      CluEError::CIFNoFractY => write!(f,
          "no \"_atom_site_fract_y\" in cif file"),

      CluEError::CIFNoFractZ => write!(f,
          "no \"_atom_site_fract_z\" in cif file"),

      CluEError::CIFNoTypeSymbol => write!(f,
          "no \"_atom_site_type_symbol\" in cif file"),

      CluEError::ConfigModeNotRecognized(mode) => write!(f,
          "#[{}] is not recognized",mode),

      CluEError::DeprecatedKeywordReplaced(line_number, 
          deprecated, replaced) => write!(f,
          "line {}, \"{}\" is deprecated and is replaced by \"{}\"", 
          line_number,deprecated, replaced),

      CluEError::EmptyVector(line_number) => write!(f,
          "line {}, supplied vector is empty", line_number),

      CluEError::Error(error_message) => write!(f,
          "{}", error_message),

      CluEError::DetectedSpinDoesNotHaveAnActiveIndex => write!(f,
          "the detected spin is always active and so does have an index \
fo nth active"),

      CluEError::ExpectedClusterSetWithNSizes(n_exp,n_act) => write!(f,
          "expected a cluster set with {} sizes, but got {} sizes",
          n_exp, n_act),

      CluEError::ExpectedEquality(line_number) => write!(f,
          "line {}, expected an equaliy",line_number),

      CluEError::ExpectedBoolRHS(line_number) => write!(f,
          "line {}, expected a bool on the right hand side",line_number),

      CluEError::ExpectedFloatRHS(line_number) => write!(f,
          "line {}, expected a float on the right hand side",line_number),

      CluEError::ExpectedIntRHS(line_number) => write!(f,
          "line {}, expected an integer on the right hand side",line_number),

      CluEError::ExpectedNonNegativeIntRHS(line_number) => write!(f,
          "line {}, expected a non-negative integer on the right hand side",
          line_number),

      CluEError::ExpectedNumber(line_number) => write!(f,
          "line {}, expected a number",line_number),

      CluEError::ExpectedTOMLArray(type_str) => write!(f,
          "expected array, but got {}",type_str),

      CluEError::ExpectedTOMLBool(type_str) => write!(f,
          "expected bool, but got {}",type_str),

      CluEError::ExpectedTOMLFloat(type_str) => write!(f,
          "expected float, but got {}",type_str),

      CluEError::ExpectedTOMLString(type_str) => write!(f,
          "expected string, but got {}",type_str),

      CluEError::ExpectedTOMLTable(type_str) => write!(f,
          "expected table, but got {}",type_str),

      CluEError::ExpectedVecOfNFloatsRHS(line_number,n) => write!(f,
          "line {}, expected a vector of {} floats on the right hand side",
          line_number,n),

      CluEError::FilterNeedsALabel => write!(f,
          "group requires a label to be set: #[group(label = LABEL)]"),

      CluEError::FilterNoMaxDistance(label) => write!(f,
          "\"#[group(label = {})]\", has no max distance.",
          label),

      CluEError::FilterNoMinDistance(label) => write!(f,
          "\"#[group(label = {})]\", has no min distance.",
          label),

      CluEError::FiltersOverlap(label0,label1) => write!(f,
          "groups \"{}\" and \"{}\" overlap: \
particles must not reside in more than one group",label0,label1),

      CluEError::InconsistentExhangeGroupActiveStatus(ex_grp_id) => write!(f,
          "in exchange group {}, spins should all be either active or inactive",
          ex_grp_id),

      CluEError::IncorrectFormattingIsotopeAbundances(line_number) 
        => write!(f,"line {}, isotope distributions are expected as
\"isotope_distribution = {{X0: p0, X1: p1 }}\", where X0 and X1 are isotopes 
and p0,p1 > 0 are abundances",line_number),

      CluEError::IncorrectNumberOfAxes(n,n_ref)=> write!(f,
          "expected {} axes, but {} were provided",n_ref, n),

      CluEError::InorrectNumberOfCellOffsets(n,n_ref)=> write!(f,
          "expected {} cell offsets, but {} were provided",n_ref, n),

      CluEError::IndexOutOfBounds(line_number,idx, len) => write!(f,
          "line {}, cannot access element {} from array of length {}", 
          line_number, idx, len),

      CluEError::InvalidAxes => write!(f,
          "invalid axes"),

      CluEError::InvalidArgument(line_number,expected_arg) => write!(f,
          "line {}, argument should be a(n) {}",line_number, expected_arg),

      CluEError::InvalidClusterPartitionKey => write!(f,
          "invalid cluster partition key"),

      CluEError::InvalidConfigFile(filename) => write!(f,
          "cannot not read config file \"{}\"", filename),

      CluEError::InvalidPulseSequence(line_number) => write!(f,
          "line {}, invalid pulse sequence",line_number),

      CluEError::InvalidSecondaryFilter(line_number, arg) => write!(f,
          "line {}, invalid secondary filter \"{}\"",line_number, arg),

      CluEError::InvalidToken(line_number,err_token) => write!(f,
          "line {}, invalid token \"{}\"",line_number, err_token),

      CluEError::IsotopeAbundancesCannotBeNormalized(line_number) => write!(f,
          "line {}, isotope abundances cannot be normalized",line_number),

      CluEError::IsotopeAbundancesMustBeNonnegative(line_number) => write!(f,
          "line {}, isotope abundances must be non-negative",line_number),

      CluEError::LenghMismatchTimepointsIncrements(n_dts,dts) => write!(f,
          "there are {} timepoint specifications, but {} time increments",
          n_dts,dts),

      CluEError::MismatchedGroupNames(name0,name1) => write!(f,
          "group names \"{}\" and \"{}\" do not match",name0,name1),

      CluEError::MissingFilter(label) => write!(f,
          "no group with label \"{}\"",label),

      CluEError::MissingFieldInCSVFile(field,file_name) => write!(f,
          "CSV file \"{}\" has no column \"{}\"",file_name,field),

      CluEError::MissingFilterArgument(line_number,fn_name) => write!(f,
          "line {}, missing group argument in \"{}()\"",line_number,fn_name),

      CluEError::MissingFilterLabel(line_number) => write!(f,
          "line {}, missing label in group",line_number),

      CluEError::MissingGroupName => write!(f,
          "no group name specified"),

      CluEError::MissingHeader(n_headers, n_cols,filename) => write!(f,
          "in \"{}\", every entry must have a header, \
but there are {} headers and {} columns of data.",filename,n_headers,n_cols),

      CluEError::MissingProperties(label) => write!(f,
          "no properties with label \"{}\"",label),

      CluEError::MissingPropertiesLabel(line_number) => write!(f,
          "line {}, missing label in properties",line_number),

      CluEError::ModeAttributeWrongBrackets => write!(f,
          "modes details should with square brackets"),

      CluEError::ModeAttributeWrongOption(mode)=> write!(f,
          "#[{}...] contains an invalid option",mode),

      CluEError::MultipleCosubstitutionGroups(index)=> write!(f,
          "particle {} is assigned to more than one cosubstitution group",
          index),

      CluEError::ModeAttributeWrongSharp => write!(f,
          "modes are specified with a single '#' at the start"),

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

      CluEError::NeighborListAndPartitionTableAreNotCompatible  => write!(f,
          "neighbor list and partition table are incompatible"),

      CluEError::NoApplyPBC => write!(f,
          "please set \"do_replicate_unit_cell: bool\" to select if \
periodic boundary conditions should be applied"),

      CluEError::NoArgument(line_number) => write!(f,
          "line {}, expected a function argument",line_number),

      CluEError::NoBathGMatrixSpecifier(label,isotope) => write!(f,
          "no g-matrix specifier found for {} {}",label,isotope),
      
      CluEError::NoCentralSpin => write!(f,
          "no detected spin defined"),
      
      CluEError::NoCentralSpinCoor => write!(f,
          "coordinates for the detected spin were not defined"),
      
      CluEError::NoCentralSpinIdentity => write!(f,
          "the detected spin's identity was not defined"),
      
      CluEError::NoCentralSpinTransition => write!(f,
          "the detected spin transition was not defined"),

      CluEError::NoClashDistancePBC => write!(f,
          "clash_distance_pbc is not defined"),

      CluEError::NoClusterBatchSize => write!(f,
          "batch size for clusters not specified"),

      CluEError::NoClusterMethod => write!(f,
          "no cluster method specified"),
      
      CluEError::NoClustersOfSize(size) => write!(f,
          "cannot find any clusters of size {}", size),
      
      CluEError::NoClusterSource => write!(f,
          "no cluster source specified"),
      
      CluEError::NoDensityMatrixMethod=> write!(f,
          "no density matrix method specified"),
      
      CluEError::NoDetectedSpinIdentity => write!(f,
          "detected_spin_identity is not set"),

      CluEError::NoDetectedSpinMultiplicity => write!(f,
          "detected_spin_identity is not set"),

      CluEError::NoDetectedSpinNotSet => write!(f,
          "detected_spin is not set"),

      CluEError::NoDetectedSpinTransition => write!(f,
          "detected_transition is not set"),

      CluEError::NoExtracellIsotopicDistribution(label) => write!(f,
          "extracel_isotopic_distribution is not set for {}",label),

      CluEError::NoGMatrixSpecifier => write!(f,
          "no g-matrix specifier"),
      
      CluEError::NoGMatrixValues => write!(f,
          "g-matrix values are not set"),
      
      CluEError::NoInputFile => write!(f,
          "no input file"),
      
      CluEError::NoHyperfineSpecifier(label,isotope) => write!(f,
          "no hyperfine specifier found for {} {}",label,isotope),
      
      CluEError::NoMagneticField => write!(f,
          "please specify the applied magnetic field"),

      CluEError::NoModelIndex => write!(f,
          "PDB model not selected"),

      CluEError::NoMaxClusterSize => write!(f,
          "maximum cluster size not set"),

      CluEError::NoNeighborCutoffDistance => write!(f,
          "neighbor_cutoff_distance is not defined"),
      
      CluEError::NoNumberSystemInstances => write!(f,
          "number_runs is not set"),

      CluEError::NoOrientationGrid => write!(f,
          "orientation_grid is not defined"),
      
      CluEError::NoPartitioningMethod => write!(f,
          "no partitioning method is defined"),
      
      CluEError::NoPulseSequence => write!(f,
          "no pulse sequence is defined"),
      
      CluEError::NoQuadrupoleSpecifier(label,isotope) => write!(f,
          "no quadrupole specifier found for {} {}",label,isotope),
      
      CluEError::NoRadius => write!(f,
          "system radius not set"),

      CluEError::NoRelationalOperators(line_number) => write!(f,
          "line {}, no relational operators (=, <, >, in, ...), are present", 
          line_number),

      CluEError::NotA3DVector(dim) => write!(f,
          "vector is {}-dimensional, not 3-dimensional", 
          dim),

      CluEError::NotALebedevGrid(n_ori) => write!(f,
          "\"{}\" is not a valid Lebedev grid; please choose n_ori in \
{{6, 14, 26, 38, 50, 74, 86, 110, 146, 170, \
194,  230, 266, 302, 350, 434, 590, 770, 974, 1202, \
1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802,5294, 5810}}",n_ori),

      CluEError::NotAnOperator(line_number, token) => write!(f,
          "line {}, cannot interpret \"{}\" as an operator", line_number,token),

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
          "no time increments defined"),
      
      CluEError::NoTimepoints => write!(f,
          "please specify how many timepoints there are for each increment"),
      
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
      
      CluEError::NoRHS(line_number) => write!(f,
          "line {}, cannot read right hand side",line_number),

      CluEError::NoSpinOpForClusterSize(cluster_size,max_size) => write!(f,
          "no spin operators for clusters of size {} are built,\
          only sizes up to {}",cluster_size,max_size),

      CluEError::NoSpinOpWithMultiplicity(spin_multiplicity) => write!(f,
          "no spin-{} operators are built", 
          (*spin_multiplicity as f64 - 1.0)/2.0),

      CluEError::NoStructureFile => write!(f,
          "no structure file defined"),

      CluEError::ParticlesClash(idx0,elmt0,idx1,elmt1,r,r_clash) => write!(f,
          "particle {} {} and particle {} {} are {} Å apart, closer than the
clash distance of {} Å",idx0,elmt0,idx1,elmt1,r,r_clash),

      CluEError::ParticleIsNotActive(ref_index) => write!(f,
          "particle \"{}\" is not active", ref_index),

      CluEError::PartitionIsIncomplette => write!(f,
          "partition table should contain n cells indexed by [0,n-1]"),

      CluEError::OptionAlreadySet(line_number,err_token) => write!(f,
          "line {}, \"{}\" has already been set",line_number, err_token),

      CluEError::UnmatchedBlockComment(line_number) => write!(f,
          "line {}, unmatched \"*/\"", line_number),

      CluEError::UnmatchedDelimiter(line_number) => write!(f,
          "line {}, unmatched delimiter", line_number),

      CluEError::SaveNameEmpty => write!(f,
          "save name is empty"),

      CluEError::SaveNameNotSet => write!(f,
          "save name not set"),

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

      CluEError::TooFewRHSArguments(line_number) => write!(f,
          "line {}, too few arguments on the right hand side", 
          line_number),

      CluEError::TooManyExchangeCouplingsSpecified(n) => write!(f,
          "expected 1 exchange coupling, but got {}", 
          n),

      CluEError::TooManyRelationalOperators(line_number) => write!(f,
          "line {}, too many relational operators (=, <, >, in, ...), are present", 
          line_number),

      CluEError::TooManyRHSArguments(line_number) => write!(f,
          "line {}, too many arguments on the right hand side", 
          line_number),

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

      CluEError::TOMLValueDoesNotSpecifyAPartitionTable => write!(f,
          "TOML value is not a partition table", 
          ),

      CluEError::TOMLValueIsNotABool => write!(f,
          "TOML value is not a bool", 
          ),

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

      CluEError::WrongNumberOfAxes(num_axes, expected_num) => write!(f,
          "{} axes were provided, but {} are expected",num_axes, expected_num),

      CluEError::WrongOrientationGridDim(line_number, expected, actual) 
        => write!(f,
          "line {}, cannot add a {}D, point to a {}D grid", 
          line_number, expected,actual),

      CluEError::WrongProbabilityDistributionDim(line_number, 
          expected, actual) => write!(f,
          "line {}, expected probability distribution with {} dimensions, \
but found {} dimensions", line_number, expected,actual),

      CluEError::WrongVectorLength(line_number, expected, actual) => write!(f,
          "line {}, expected vector of length {}, but recieved a length of {}", 
          line_number, expected,actual),
    }
  }
}
