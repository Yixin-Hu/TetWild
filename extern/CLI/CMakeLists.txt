add_library(cli11_tiny INTERFACE)
target_include_directories(cli11_tiny INTERFACE ..)
add_library(CLI11::CLI11 ALIAS cli11_tiny)
