//! The structure files of a run: the entries of `structure.file` (paths,
//! globs, and weights csv files, each with a weight) and their expansion
//! into a list of files with weights (`Config::expand_structure_files()`).
//! The expansion does not read the structure files; their contents are
//! hashed and grouped by `Config::structure_inputs()`.
use crate::clue_errors::CluEError;
use crate::config::Config;
use crate::config::validate::check_path;

use std::collections::HashMap;
use std::path::{Path,PathBuf};

/// The column of a weights csv file that holds the structure files.
pub const WEIGHTS_CSV_FILE_COLUMN: &str = "file";
/// The column of a weights csv file that holds the weights.
pub const WEIGHTS_CSV_WEIGHT_COLUMN: &str = "weight";

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// `StructureFileSpec` is one entry of `structure.file`: a path, a glob,
/// or a weights csv file (`pattern`), and its weight, if one is given.
/// The weight of a glob applies to each file it matches.
#[derive(Debug,Clone,PartialEq)]
pub struct StructureFileSpec{
  pub pattern: String,
  pub weight: Option<f64>,
}

impl StructureFileSpec{
  //----------------------------------------------------------------------------
  /// This function returns the entries of the parsed `structure.file`:
  /// a string, or an array whose elements are strings or tables with a
  /// string `file` and an optional number `weight`.
  /// An entry without a weight has `weight: None` (weight 1).
  /// Values of other types are errors that name the element, such as
  /// `structure.file[1]`. Other keys of the tables are reported by the 
  /// key check (`check_toml_keys()`) before this function is called.
  pub fn from_toml(file: toml::Value) 
    -> Result<Vec::<StructureFileSpec>,CluEError>
  {
    const KEY: &str = "structure.file";
    let entries = match file{
      toml::Value::String(pattern) 
        => return Ok(vec![StructureFileSpec{pattern, weight: None}]),
      toml::Value::Array(entries) => entries,
      other => return Err(CluEError::InvalidTOMLValue(KEY.to_string(),
            "a path, a glob, a weights csv file, or an array of them (each \
a string or a table {file = ..., weight = ...})".to_string(),
            other.type_str().to_string())),
    };

    let mut specs = Vec::<StructureFileSpec>::with_capacity(entries.len());
    for (ii, entry) in entries.into_iter().enumerate(){
      let path = format!("{}[{}]", KEY, ii);
      let mut table = match entry{
        toml::Value::String(pattern) => {
          specs.push(StructureFileSpec{pattern, weight: None});
          continue;
        },
        toml::Value::Table(table) => table,
        other => return Err(CluEError::InvalidTOMLValue(path,
              "a path or a table {file = ..., weight = ...}".to_string(),
              other.type_str().to_string())),
      };
      let pattern = match table.remove("file"){
        Some(toml::Value::String(pattern)) => pattern,
        Some(other) => return Err(CluEError::InvalidTOMLValue(
              format!("{}.file", path), "a string".to_string(),
              other.type_str().to_string())),
        None => return Err(CluEError::MissingTOMLKey(
              format!("{}.file", path), format!("the table {}", path))),
      };
      let weight = match table.remove("weight"){
        None => None,
        Some(toml::Value::Integer(w)) => Some(w as f64),
        Some(toml::Value::Float(w)) => Some(w),
        Some(value) => return Err(CluEError::InvalidStructureWeight(
              pattern, value.to_string())),
      };
      specs.push(StructureFileSpec{pattern, weight});
    }
    Ok(specs)
  }
  //----------------------------------------------------------------------------
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// `StructureSource` is one structure file of a run and its relative
/// weight, as given in the input (not normalized).
#[derive(Debug,Clone,PartialEq)]
pub struct StructureSource{
  pub file: String,
  pub weight: f64,
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

impl Config{
  //----------------------------------------------------------------------------
  /// This function sets `structure.file` to the single entry `pattern`
  /// (a path, glob, or weights csv file) without a weight, and clears the
  /// structures read by `resolve_structures()`.
  pub fn set_structure_file(&mut self, pattern: &str){
    self.structure_files = Some(vec![StructureFileSpec{
      pattern: pattern.to_string(), weight: None}]);
    self.resolved_structures = None;
  }
  //----------------------------------------------------------------------------
  /// This function expands the entries of `structure.file` into the list
  /// of structure files with their weights, in input order:
  /// a weights csv file (a path ending in `.csv`) is replaced by its rows,
  /// and a glob (a pattern with `*`, `?`, or `[` that is not the name of 
  /// an existing file) by the files it matches, sorted lexicographically;
  /// any other path must be an existing file.
  /// Weights must be finite and greater than 0; the weight of a glob
  /// applies to each file it matches.
  /// A file that is given more than once, also by different paths to the
  /// same file (`fs::canonicalize()`), is an error.
  /// Relative paths are relative to the working directory.
  /// No file at all (`file = []`) is the error `NoStructureFile`.
  /// The structure files are not read; a weights csv file is.
  pub fn expand_structure_files(&self)
    -> Result<Vec::<StructureSource>,CluEError>
  {
    let Some(specs) = &self.structure_files else{
      return Err(CluEError::NoStructureFile);
    };

    // Each source with the entry that gave it, for duplicate messages.
    let mut sources = Vec::<StructureSource>::new();
    let mut entries = Vec::<String>::new();

    for spec in specs.iter(){
      if is_weights_csv(&spec.pattern){
        if spec.weight.is_some(){
          return Err(CluEError::StructureWeightsFileWithWeight(
                spec.pattern.clone()));
        }
        for (pattern, weight) in read_structure_weights(&spec.pattern)?{
          if is_weights_csv(&pattern){
            return Err(CluEError::NestedStructureWeightsFile(
                  spec.pattern.clone(), pattern));
          }
          let entry = format!("{} (in {})", pattern, spec.pattern);
          check_weight(&entry, weight)?;
          for file in expand_pattern(&pattern)?{
            sources.push(StructureSource{file, weight});
            entries.push(entry.clone());
          }
        }
      }else{
        let weight = spec.weight.unwrap_or(1.0);
        check_weight(&spec.pattern, weight)?;
        for file in expand_pattern(&spec.pattern)?{
          sources.push(StructureSource{file, weight});
          entries.push(spec.pattern.clone());
        }
      }
    }

    // An empty array, or weights files without rows.
    if sources.is_empty(){
      return Err(CluEError::NoStructureFile);
    }

    check_duplicates(&sources, &entries)?;

    Ok(sources)
  }
  //----------------------------------------------------------------------------
}

//------------------------------------------------------------------------------
// This function returns `true` if `pattern` names a weights csv file:
// its extension is `csv`, in any case.
fn is_weights_csv(pattern: &str) -> bool{
  pattern.to_lowercase().ends_with(".csv")
}
//------------------------------------------------------------------------------
// This function returns `true` if `pattern` is a glob.
fn is_glob(pattern: &str) -> bool{
  pattern.contains(['*', '?', '['])
}
//------------------------------------------------------------------------------
// This function returns an error if `weight` of the entry `entry` is not
// finite or not greater than 0.
fn check_weight(entry: &str, weight: f64) -> Result<(),CluEError>{
  if weight.is_finite() && weight > 0.0{
    Ok(())
  }else{
    Err(CluEError::InvalidStructureWeight(entry.to_string(),
          weight.to_string()))
  }
}
//------------------------------------------------------------------------------
// This function returns the files of `pattern`: the files that the glob
// `pattern` matches, sorted lexicographically, or `pattern` itself, which
// must be an existing file.
// A pattern that names an existing file is that file, also if it contains
// glob characters, such as `frame[1].pdb`.
// Directories that a glob matches are skipped.
fn expand_pattern(pattern: &str) -> Result<Vec::<String>,CluEError>{
  if !is_glob(pattern) || Path::new(pattern).is_file(){
    if let Some(problem) = check_path(pattern, false){
      return Err(CluEError::InvalidStructureFile(problem));
    }
    return Ok(vec![pattern.to_string()]);
  }

  let paths = match glob::glob(pattern){
    Ok(paths) => paths,
    Err(err) => return Err(CluEError::InvalidStructureGlob(
          pattern.to_string(), err.to_string())),
  };

  let mut files = Vec::<String>::new();
  for path in paths{
    let path = match path{
      Ok(path) => path,
      Err(err) => return Err(CluEError::InvalidStructureGlob(
            pattern.to_string(), err.to_string())),
    };
    if !path.is_file(){
      continue;
    }
    match path.to_str(){
      Some(file) => files.push(file.to_string()),
      None => return Err(CluEError::InvalidStructureGlob(pattern.to_string(),
            format!("the path {} is not valid UTF-8", path.display()))),
    }
  }

  if files.is_empty(){
    return Err(CluEError::StructureGlobNoMatch(pattern.to_string()));
  }

  files.sort();
  Ok(files)
}
//------------------------------------------------------------------------------
// This function reads the weights csv file `filename`: a header with the
// columns `file` and `weight`, then one structure file (a path or a glob)
// per row. A blank weight is 1. Lines starting with `#` are comments, 
// blank lines are skipped, and spaces around entries are ignored. 
// A quoted entry cannot span lines.
// Errors name the line of the file, counting comment and blank lines.
fn read_structure_weights(filename: &str) 
  -> Result<Vec::<(String,f64)>,CluEError>
{
  let text = match std::fs::read_to_string(filename){
    Ok(text) => text,
    Err(err) => return Err(CluEError::CannotOpenFile(
        filename.to_string(), err.to_string())),
  };
  // Comment and blank lines are removed here, not by the csv reader, 
  // whose line numbers do not count the lines it skips. `file_lines[i]` is
  // the line of the file of line i + 1 of `data`. A byte-order mark is
  // removed.
  let mut file_lines = Vec::<usize>::new();
  let mut kept = Vec::<&str>::new();
  for (ii, line) in text.strip_prefix('\u{feff}').unwrap_or(&text).lines()
    .enumerate()
  {
    if line.starts_with('#') || line.trim().is_empty(){
      continue;
    }
    file_lines.push(ii + 1);
    kept.push(line);
  }
  if kept.is_empty(){
    return Err(CluEError::CannotReadStructureWeights(filename.to_string(), 1,
          format!("the file has no header \"{},{}\"", 
            WEIGHTS_CSV_FILE_COLUMN, WEIGHTS_CSV_WEIGHT_COLUMN)));
  }
  let data = kept.join("\n");
  let mut rdr = csv::ReaderBuilder::new()
    .has_headers(true)
    .trim(csv::Trim::All)
    .from_reader(data.as_bytes());

  let line_at = |pos: &csv::Position| file_lines.get((pos.line() as usize)
      .saturating_sub(1)).copied().unwrap_or(0);
  let error = |line: usize, msg: String| CluEError::CannotReadStructureWeights(
      filename.to_string(), line, msg);
  let read_error = |err: csv::Error| {
    let line = err.position().map_or(0, line_at);
    match err.kind(){
      csv::ErrorKind::UnequalLengths{expected_len, len, ..} => error(line,
          format!("the row has {} entries, but the header has {}", len,
            expected_len)),
      _ => error(line, err.to_string()),
    }
  };

  let header = match rdr.headers(){
    Ok(header) => header.clone(),
    Err(err) => return Err(read_error(err)),
  };
  let header_line = header.position().map_or(1, line_at);
  let headers = header.iter().collect::<Vec::<&str>>();
  let expected = [WEIGHTS_CSV_FILE_COLUMN, WEIGHTS_CSV_WEIGHT_COLUMN];
  if headers != expected{
    return Err(error(header_line, format!("the header must be \"{}\", not \
\"{}\"", expected.join(","), headers.join(","))));
  }

  let mut rows = Vec::<(String,f64)>::new();
  for result in rdr.records(){
    // Rows with a different number of entries are reading errors.
    let record = match result{
      Ok(record) => record,
      Err(err) => return Err(read_error(err)),
    };
    let line = record.position().map_or(0, line_at);

    let file = record[0].to_string();
    if file.is_empty(){
      return Err(error(line, "the file is blank".to_string()));
    }
    let weight = if record[1].is_empty(){
      1.0
    }else{
      match record[1].parse::<f64>(){
        Ok(weight) => weight,
        Err(_) => return Err(error(line, format!(
                "cannot read the weight \"{}\" as a number", &record[1]))),
      }
    };
    rows.push((file, weight));
  }
  Ok(rows)
}
//------------------------------------------------------------------------------
// This function returns an error if two of `sources` are the same file
// (`fs::canonicalize()`); `entries[i]` is the entry of `structure.file`
// that gave `sources[i]`.
fn check_duplicates(sources: &[StructureSource], entries: &[String])
  -> Result<(),CluEError>
{
  let mut seen = HashMap::<PathBuf,usize>::with_capacity(sources.len());
  for (ii, source) in sources.iter().enumerate(){
    let path = match std::fs::canonicalize(&source.file){
      Ok(path) => path,
      Err(err) => return Err(CluEError::CannotOpenFile(source.file.clone(),
            err.to_string())),
    };
    if let Some(&first) = seen.get(&path){
      return Err(CluEError::DuplicateStructureFile(
            path.display().to_string(), entries[first].clone(),
            entries[ii].clone()));
    }
    seen.insert(path, ii);
  }
  Ok(())
}
//------------------------------------------------------------------------------


//==============================================================================
#[cfg(test)]
mod tests{
  use super::*;
  use crate::io::FromTOMLString;

  //----------------------------------------------------------------------------
  // This function returns an empty directory for test `name` under the 
  // system temp directory.
  fn temp_dir(name: &str) -> String{
    let dir = std::env::temp_dir().join(format!(
          "clue_oxide_test_structure_files_{}_{}", name, std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir.to_str().unwrap().to_string()
  }
  //----------------------------------------------------------------------------
  // This function returns the config of the [structure] table `structure`.
  fn config_of(structure: &str) -> Result<Config,CluEError>{
    Config::from_toml_string(&format!("[structure]\n{}\n", structure))
  }
  //----------------------------------------------------------------------------
  fn spec(pattern: &str, weight: Option<f64>) -> StructureFileSpec{
    StructureFileSpec{pattern: pattern.to_string(), weight}
  }
  //----------------------------------------------------------------------------
  fn source(file: &str, weight: f64) -> StructureSource{
    StructureSource{file: file.to_string(), weight}
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_structure_file_syntax(){
    let files = |structure: &str| config_of(structure).unwrap()
      .structure_files.unwrap();

    assert_eq!(files("file = \"a.pdb\""), vec![spec("a.pdb", None)]);
    assert_eq!(files("file = \"frames/*.pdb\""), 
        vec![spec("frames/*.pdb", None)]);
    assert_eq!(files("file = [\"a.pdb\", \"b.pdb\"]"), 
        vec![spec("a.pdb", None), spec("b.pdb", None)]);
    // The mixed array of decision 1; integer weights are accepted.
    assert_eq!(files("file = [{file = \"file1.pdb\", weight = 0.1}, \
{file = \"path/*.pdb\", weight = 0.001}, \"file3.pdb\", \
{file = \"file4.pdb\", weight = 2}, {file = \"file5.pdb\"}]"), 
        vec![spec("file1.pdb", Some(0.1)), spec("path/*.pdb", Some(0.001)),
          spec("file3.pdb", None), spec("file4.pdb", Some(2.0)),
          spec("file5.pdb", None)]);
    assert_eq!(files("file = \"weights.csv\""), 
        vec![spec("weights.csv", None)]);

    // A weight that is not a number.
    assert_eq!(config_of("file = [{file = \"a.pdb\", weight = \"2\"}]")
        .unwrap_err(), 
        CluEError::InvalidStructureWeight("a.pdb".to_string(), 
          "\"2\"".to_string()));

    // Values of other types name the element (audit F2).
    let table_or_path = "a path or a table {file = ..., weight = ...}";
    let invalid = |key: &str, expected: &str, found: &str| 
      CluEError::InvalidTOMLValue(key.to_string(), expected.to_string(), 
          found.to_string());
    assert_eq!(config_of("file = [{file = \"a.pdb\"}, 3]").unwrap_err(),
        invalid("structure.file[1]", table_or_path, "integer"));
    assert_eq!(config_of("file = [[\"a.pdb\"]]").unwrap_err(),
        invalid("structure.file[0]", table_or_path, "array"));
    assert_eq!(config_of("file = [{file = 3}]").unwrap_err(),
        invalid("structure.file[0].file", "a string", "integer"));
    assert_eq!(config_of("file = [{weight = 2}]").unwrap_err(),
        CluEError::MissingTOMLKey("structure.file[0].file".to_string(),
          "the table structure.file[0]".to_string()));
    let whole = "a path, a glob, a weights csv file, or an array of them \
(each a string or a table {file = ..., weight = ...})";
    assert_eq!(config_of("file = {file = \"a.pdb\"}").unwrap_err(),
        invalid("structure.file", whole, "table"));
    assert_eq!(config_of("file = 3").unwrap_err(),
        invalid("structure.file", whole, "integer"));
    assert_eq!(config_of("file = [{file = \"a.pdb\", weight = true}]")
        .unwrap_err(), 
        CluEError::InvalidStructureWeight("a.pdb".to_string(), 
          "true".to_string()));
    // The message names the element and the allowed forms.
    assert_eq!(config_of("file = [\"a.pdb\", 3]").unwrap_err().to_string(),
        "structure.file[1] must be a path or a table {file = ..., \
weight = ...}, but got integer");

    // Unknown keys of an entry table are reported by the key check.
    match config_of("file = [\"b.pdb\", {file = \"a.pdb\", wieght = 2}]"){
      Err(CluEError::InvalidTOMLKeys(errors)) => {
        assert_eq!(errors.len(), 1);
        assert!(errors[0].starts_with("structure.file[1].wieght"), 
            "{}", errors[0]);
      },
      other => panic!("{:?}", other),
    }
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_expand_structure_files(){
    let dir = temp_dir("expand");
    let path = |name: &str| format!("{}/{}", dir, name);
    for (name, contents) in [("f_1.pdb", "1"), ("f_2.pdb", "2"), 
      ("f_10.pdb", "10"), ("g.pdb", "1")]
    {
      std::fs::write(path(name), contents).unwrap();
    }
    // Directories that a glob matches are skipped.
    std::fs::create_dir(path("f_dir.pdb")).unwrap();
    // A file whose name contains glob characters.
    std::fs::write(path("h[1].pdb"), "h").unwrap();
    std::fs::write(path("h1.pdb"), "h1").unwrap();

    let mut config = Config::new();
    let expand = |config: &Config| config.expand_structure_files();

    // A glob is sorted lexicographically, and its weight applies to each
    // file.
    config.structure_files = Some(vec![spec(&path("f_*.pdb"), Some(2.0)),
      spec(&path("g.pdb"), None)]);
    assert_eq!(expand(&config).unwrap(), vec![source(&path("f_1.pdb"), 2.0),
      source(&path("f_10.pdb"), 2.0), source(&path("f_2.pdb"), 2.0),
      source(&path("g.pdb"), 1.0)]);

    // An existing file name with glob characters is that file (audit F4);
    // otherwise the pattern is a glob.
    config.set_structure_file(&path("h[1].pdb"));
    assert_eq!(expand(&config).unwrap(), vec![source(&path("h[1].pdb"), 1.0)]);
    config.set_structure_file(&path("h[12].pdb"));
    assert_eq!(expand(&config).unwrap(), vec![source(&path("h1.pdb"), 1.0)]);

    // A glob without a match, an invalid glob, and a missing file.
    config.set_structure_file(&path("h_*.pdb"));
    assert_eq!(expand(&config).unwrap_err(), 
        CluEError::StructureGlobNoMatch(path("h_*.pdb")));
    config.set_structure_file(&path("f_[.pdb"));
    assert!(matches!(expand(&config).unwrap_err(), 
          CluEError::InvalidStructureGlob(..)));
    config.set_structure_file(&path("h.pdb"));
    assert_eq!(expand(&config).unwrap_err(), 
        CluEError::InvalidStructureFile(format!("\"{}\" does not exist", 
            path("h.pdb"))));
    config.set_structure_file(&path("f_dir.pdb"));
    assert_eq!(expand(&config).unwrap_err(), 
        CluEError::InvalidStructureFile(format!("\"{}\" is not a file", 
            path("f_dir.pdb"))));

    // Duplicates: an explicit file that a glob also matches, and two
    // paths to the same file.
    config.structure_files = Some(vec![spec(&path("f_2.pdb"), None),
      spec(&path("f_*.pdb"), None)]);
    let canonical = std::fs::canonicalize(path("f_2.pdb")).unwrap();
    assert_eq!(expand(&config).unwrap_err(), 
        CluEError::DuplicateStructureFile(canonical.display().to_string(),
          path("f_2.pdb"), path("f_*.pdb")));
    config.structure_files = Some(vec![spec("assets/TEMPO.pdb", None),
      spec("./assets/TEMPO.pdb", None)]);
    assert!(matches!(expand(&config).unwrap_err(), 
          CluEError::DuplicateStructureFile(..)));

    // Weights must be finite and greater than 0.
    for weight in [0.0, -1.0, f64::NAN, f64::INFINITY]{
      config.structure_files = Some(vec![spec(&path("g.pdb"), Some(weight))]);
      assert_eq!(expand(&config).unwrap_err(), 
          CluEError::InvalidStructureWeight(path("g.pdb"), 
            weight.to_string()));
    }

    // No file at all.
    config.structure_files = Some(Vec::<StructureFileSpec>::new());
    assert_eq!(expand(&config).unwrap_err(), CluEError::NoStructureFile);
    config.structure_files = None;
    assert_eq!(expand(&config).unwrap_err(), CluEError::NoStructureFile);

    let _ = std::fs::remove_dir_all(&dir);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_structure_weights_csv(){
    let dir = temp_dir("weights_csv");
    let path = |name: &str| format!("{}/{}", dir, name);
    for name in ["a.pdb", "b_1.pdb", "b_2.pdb"]{
      std::fs::write(path(name), name).unwrap();
    }
    let write_csv = |name: &str, text: &str| {
      std::fs::write(path(name), text).unwrap();
      path(name)
    };

    // Comments, blank weights, globs, and spaces around entries.
    let csv = write_csv("w.CSV", &format!("# frames\nfile,weight\n\
{}, 0.5\n# b frames\n{},\n", path("a.pdb"), path("b_*.pdb")));
    let mut config = Config::new();
    config.structure_files = Some(vec![spec(&csv, None)]);
    assert_eq!(config.expand_structure_files().unwrap(), vec![
      source(&path("a.pdb"), 0.5), source(&path("b_1.pdb"), 1.0),
      source(&path("b_2.pdb"), 1.0)]);

    // A weight on a weights file.
    config.structure_files = Some(vec![spec(&csv, Some(1.0))]);
    assert_eq!(config.expand_structure_files().unwrap_err(),
        CluEError::StructureWeightsFileWithWeight(csv.clone()));

    // A weights file in a weights file.
    let nested = write_csv("nested.csv", &format!("file,weight\n{},2\n", csv));
    config.set_structure_file(&nested);
    assert_eq!(config.expand_structure_files().unwrap_err(),
        CluEError::NestedStructureWeightsFile(nested.clone(), csv.clone()));

    // Reading errors name the line.
    let read_error = |text: &str| {
      let mut config = Config::new();
      config.set_structure_file(&write_csv("bad.csv", text));
      config.expand_structure_files().unwrap_err()
    };
    assert!(matches!(read_error("path,weight\na.pdb,1\n"), 
          CluEError::CannotReadStructureWeights(_, 1, _)));
    assert!(matches!(read_error(&format!("file,weight\n{},x\n", 
              path("a.pdb"))), 
          CluEError::CannotReadStructureWeights(_, 2, _)));
    assert!(matches!(read_error("file,weight\n,1\n"), 
          CluEError::CannotReadStructureWeights(_, 2, _)));
    assert!(matches!(read_error(&format!("file,weight\n{},1,2\n", 
              path("a.pdb"))), 
          CluEError::CannotReadStructureWeights(_, 2, _)));
    assert_eq!(read_error(&format!("file,weight\n{},0\n", path("a.pdb"))),
        CluEError::InvalidStructureWeight(format!("{} (in {})", 
            path("a.pdb"), path("bad.csv")), "0".to_string()));

    // The lines are those of the file, counting comment and blank lines
    // (audit F9), also with CRLF line ends and a byte-order mark.
    let weights_error = |line: usize, msg: &str| 
      CluEError::CannotReadStructureWeights(path("bad.csv"), line, 
          msg.to_string());
    let not_a_number = "cannot read the weight \"x\" as a number";
    assert_eq!(read_error(&format!("# a\n# b\nfile,weight\n# c\n\n{},x\n", 
            path("a.pdb"))), weights_error(6, not_a_number));
    assert_eq!(read_error("# a\n\n  \nfile,wait\n"), 
        weights_error(4, "the header must be \"file,weight\", not \
\"file,wait\""));
    assert_eq!(read_error(&format!("# a\nfile,weight\n\n{},1,2\n", 
            path("a.pdb"))), 
        weights_error(4, "the row has 3 entries, but the header has 2"));
    assert_eq!(read_error(&format!("file,weight\r\n# c\r\n{},x\r\n", 
            path("a.pdb"))), weights_error(3, not_a_number));
    assert_eq!(read_error(&format!("\u{feff}# c\nfile,weight\n{},x\n", 
            path("a.pdb"))), weights_error(3, not_a_number));
    assert_eq!(read_error("# only a comment\n"), 
        weights_error(1, "the file has no header \"file,weight\""));
    // A byte-order mark before the header, and a line of spaces.
    let csv = write_csv("bom.csv", &format!("\u{feff}file,weight\n  \n{},\n",
          path("a.pdb")));
    config.set_structure_file(&csv);
    assert_eq!(config.expand_structure_files().unwrap(), 
        vec![source(&path("a.pdb"), 1.0)]);

    let _ = std::fs::remove_dir_all(&dir);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_structure_inputs_group_identical_contents(){
    let dir = temp_dir("group");
    let path = |name: &str| format!("{}/{}", dir, name);
    for (name, contents) in [("a.pdb", "A"), ("b.pdb", "B"), ("c.pdb", "A")]{
      std::fs::write(path(name), contents).unwrap();
    }
    let mut config = Config::new();
    config.structure_files = Some(vec![spec(&path("a.pdb"), Some(1.0)),
      spec(&path("b.pdb"), Some(3.0)), spec(&path("c.pdb"), Some(0.5))]);

    let structures = config.structure_inputs().unwrap();
    assert_eq!(structures.len(), 2);
    assert_eq!(structures[0].file, path("a.pdb"));
    assert_eq!(structures[0].weight, 1.5);
    assert_eq!(structures[0].sources, vec![source(&path("a.pdb"), 1.0),
      source(&path("c.pdb"), 0.5)]);
    assert_eq!(structures[1].file, path("b.pdb"));
    assert_eq!(structures[1].weight, 3.0);
    assert_eq!(structures[1].sources, vec![source(&path("b.pdb"), 3.0)]);
    assert_ne!(structures[0].content_hash, structures[1].content_hash);

    // The structures are read once by resolve_structures().
    config.resolve_structures().unwrap();
    std::fs::remove_file(path("b.pdb")).unwrap();
    assert_eq!(config.structure_inputs().unwrap(), structures);
    // Without the cache, the missing file is an error.
    let mut unresolved = config.clone();
    unresolved.resolved_structures = None;
    assert!(unresolved.structure_inputs().is_err());
    // A new structure.file clears the cache.
    config.set_structure_file(&path("a.pdb"));
    assert_eq!(config.structure_inputs().unwrap().len(), 1);
    assert_eq!(config.structure_inputs().unwrap()[0].weight, 1.0);

    let _ = std::fs::remove_dir_all(&dir);
  }
  //----------------------------------------------------------------------------
}
