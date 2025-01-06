#ifndef HAVE_INPUT_LOADER_H
#define HAVE_INPUT_LOADER_H   

#include "utils.h"

void load_input(int argc, char *argv[], molecule_params *molecule, field_params *field, solver_params *params, avg_params *avg);

void load_command_line_arguments(int argc, char *argv[], molecule_params *molecule, field_params *field, solver_params *params, avg_params *avg);

void load_from_file(char *filename, molecule_params *molecule, field_params *field, solver_params *params, avg_params *avg);

void load_directly(molecule_params *molecule, field_params *field, solver_params *params, avg_params *avg);

void get_custom_pulse_spline(char *filename, solver_params *params);

#endif // HAVE_INPUT_LOADER_H