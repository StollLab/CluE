use crate::info::version;

/// This function prints the title info. 
pub fn print_title() {  
println!("\
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
  println!("CLuE Copyright (C) 2022--2024 Samuel M. Jahn");
  version::print_version();
  println!("\nThis program comes with ABSOLUTELY NO WARRANTY; for details run \n\
\"clue_oxide --warrenty\".");
  println!("\nThis is free software distributed under the GPL v3.0 license; for details run ");
  println!("\"clue_oxide --license\".");
  println!("\nFor usage information, run \n\"clue_oxide --help\".\n");


}
