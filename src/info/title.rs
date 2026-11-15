use crate::clue_println;
use crate::info::version;

/// This function prints the title info. 
pub fn print_title() {  
clue_println!("\
##############################################################################\n\
#                                                                            #\n\
#                      ______    _             ________                      #\n\
#                     /  ____|  | |           | _______|                     #\n\
#                    |  |       | |   _   _   | |_____                       #\n\
#                    |  |       | |  | | | |  |  _____|                      #\n\
#                    |  |____   | |  | |_| |  | |______                      #\n\
#                     \\______|  |_|  \\_____/  |________|                     #\n\
#                                                                            #\n\
#                 ___               ______               ___                 #\n\
#                    \\             /      \\             /                    #\n\
#                     \\           /        \\           /                     #\n\
#                      \\    O    /          \\    O    /                      #\n\
#                       \\       /            \\       /                       #\n\
#                        \\     /              \\     /                        #\n\
#                ()       \\   /       ()       \\   /       ()                #\n\
#                          \\ /                  \\ /                          #\n\
#                           V                    V                           #\n\
#                                                                            #\n\
#                                                                            #\n\
##############################################################################\n");
  clue_println!("CLuE Copyright (C) 2022--2024 Samuel M. Jahn");
  version::print_version();
  clue_println!("\nThis program comes with ABSOLUTELY NO WARRANTY; for details run \n\
\"clue_oxide --warranty\".");
  clue_println!("\nThis is free software distributed under the GPL v3.0 license; for details run ");
  clue_println!("\"clue_oxide --license\".");
  clue_println!("\nFor usage information, run \n\"clue_oxide --help\".\n");


}
