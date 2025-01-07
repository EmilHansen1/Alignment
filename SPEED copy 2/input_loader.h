#ifndef HAVE_INPUT_LOADER_H
#define HAVE_INPUT_LOADER_H   

#include "utils.h"

solver_params load_input(int argc, char *argv[]);

solver_params load_command_line_arguments(int argc, char *argv[]);

solver_params load_from_file(char *filename);

solver_params load_directly();

void get_custom_pulse_spline(solver_params *params);

#endif // HAVE_INPUT_LOADER_H