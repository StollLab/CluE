use crate::clue_errors::CluEError;
use substring::Substring;

/// The struct contains the processed command line input.
pub struct CommandLineInput{
  pub config_file: Option<String>,
  //pub output_file: Option<String>,
  pub config_options: Option<String>,
  pub show_help: bool,
  pub show_license: bool,
  pub show_title: bool,
  pub show_version: bool,
  pub show_warranty: bool,
  next_arg: NextArg
}

enum NextArg{
  ConfigOption,
  InputConfig,
  //Output,
  Unknown,
}

impl CommandLineInput{

  /// This function creates a new instance of `CommandLineInput`, with 
  /// field values from the input `args`.
  pub fn new(args: Vec::<String>) -> Result<Self,CluEError> {
    let mut cli = CommandLineInput{
      config_file: None,
      //output_file: None,
      config_options: None,
      show_help: false,
      show_license: false,
      show_title: true,
      show_version: false,
      show_warranty: false,
      next_arg: NextArg::InputConfig,
    };

    let mut it = args.iter();
    it.next();
    for arg in it{

      if (*arg).substring(0,2)== "--"{
        cli.parse_long_option(arg)?;
        continue;
      }
      
      if(*arg).substring(0,1)== "-"{
        cli.parse_short_options((*arg).substring(1,(*arg).len()))?;
        continue;
      }

      cli.parse_next_argument(arg)?;
    }
    Ok(cli)
  }
  //----------------------------------------------------------------------------
  // This function sets the `next_arg` field.
  fn parse_next_argument(&mut self, arg: &str) -> Result<(),CluEError>{
    match self.next_arg{
      NextArg::ConfigOption => {
        let opt: String = match &self.config_options{
          Some(opt_str) => opt_str.clone(),
          None => String::new(),
        };
        self.config_options = Some(format!("{}\n{}",opt,arg));
      }
      NextArg::InputConfig => self.config_file = Some(arg.to_string()),
      //NextArg::Output => self.output_file = Some(arg.to_string()),
      NextArg::Unknown =>
        return Err(CluEError::UnrecognizedOption(arg.to_string() )),
    }

    if self.config_file.is_none() {
      self.next_arg = NextArg::InputConfig;
    //}else if self.output_file.is_none(){
      //self.next_arg = NextArg::Output;
    }else{
      self.next_arg = NextArg::Unknown;
    }
    Ok(())
  }
  //----------------------------------------------------------------------------
  // This function updates `self` from a single word option.
  fn parse_long_option(&mut self, option: &str) -> Result<(),CluEError>{
  
    match option {
      "--help" => self.show_help = true,
      "--license" => self.show_license = true,
      "--hide-title" => self.show_title = false,
      //"--output" => self.next_arg = NextArg::Output,
      "--option" => self.next_arg = NextArg::ConfigOption,
      "--version" => self.show_version = true,
      "--warranty" => self.show_warranty = true,
      _ => return Err(CluEError::UnrecognizedOption(option.to_string() )),  
    }
    Ok(())
  }
  //----------------------------------------------------------------------------
  // This function updates `self` from a string of single character options.
  fn parse_short_options(&mut self, options: &str) ->  Result<(),CluEError>{
    
    let n_opts = options.len();

    for ii in 0..n_opts{
      let opt = options.substring(ii,ii+1);

      self.parse_short_option(opt)?;
    }
    Ok(())
  
  }
  //----------------------------------------------------------------------------
  // This function updates `self` from single character option.
  fn parse_short_option(&mut self, option: &str) ->  Result<(),CluEError>{

    match option{
      "h" => self.show_help = true,
      "l" => self.show_license = true,
      "H" => self.show_title = false,
      //"o" => self.next_arg = NextArg::Output,
      "O" => self.next_arg = NextArg::ConfigOption,
      "V" => self.show_version = true,
      "W" => self.show_warranty = true,
      _ => return Err(CluEError::UnrecognizedOption(option.to_string() )),  
    }
    Ok(())
  }  
  //----------------------------------------------------------------------------
}
