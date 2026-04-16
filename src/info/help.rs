/// This function prints the help message.
pub fn print_help(){
  println!("\
Usage:\n    \
    clue_oxide input [-o \"options\"]\n\
\n\
Options:\n    \
    -h, --help       Prints help information.\n    \
    -H, --hide-title Do not display title information.\n    \
    -l, --license    Prints license.\n    \
    -O, --option     input additional options.\n    \
    -V, --version    Prints version information.\n    \
    -W, --warranty   Prints warranty information.\n    \
");
}
