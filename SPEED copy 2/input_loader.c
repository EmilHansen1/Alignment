#include <assert.h>
#include <stdio.h>
#include <math.h>
#include "input_loader.h"
#include "planar_rotor.h"
#include "units.h"
#include "utils.h"

#define NUM_COMMAND_LINE_ARGS 69
#define MAX_CC_LENGTH 1024

/// @brief Load in all parameters for the simulation into the parameter structs
/// @param argc Command line argument length
/// @param argv Command line arguments
/// @param molecule The struct for containing molecular parameters
/// @param field The struct for contatinig field parameters
/// @param params The struct for containing all simulation parameters
void load_input(int argc, char *argv[], molecule_params molecule, field_params field, solver_params params, avg_params avg)
{
    if (argc == 1)
    {
        // Debug mode - hardcoded parameters from function
        load_directly(molecule, field, params, avg);
    }
    else if (argc == 2)
    {
        // Load parameters from file
        load_from_file(argv[1], molecule, field, params, avg);
    } 
    else 
    {
        // Load parameters from command line arguments
        assert(argc == NUM_COMMAND_LINE_ARGS);
        load_command_line_arguments(argc, argv, molecule, field, params, avg);
    }
}



void load_command_line_arguments(int argc, char *argv[], molecule_params molecule, field_params field, solver_params params, avg_params avg)
{
    // Load data into the structs based on the command line arguments
    //molecule->j_max = (size_t) atoi(argv[1]);
}



void load_from_file(char *filename, molecule_params molecule, field_params field, solver_params params, avg_params avg)
{
    // ... To be implemented
}


/// @brief Super cursed function for loading in parameters directly w/o CLA or input file
/// @brief Intended for debugging purposes
/// @param molecule 
/// @param field 
/// @param params 
void load_directly(molecule_params molecule, field_params field, solver_params params, avg_params avg)
{

    size_t j_max = (size_t) 30;  
    double B = 2.0 * M_PI * 0.69 * 1e9 * au_frequency_per_SI; // (3Li 5.25) (3Na 1.33) (3K 0.69) (3Rb 0.28) (1Na 4.60) (1.78 1K) (0.69 1Rb)
    double delta_alpha = 68.54 * pow(au_per_angstrom, 3.0); // (3Li 65.73) (3Na 31.06) (3K 68.54) (3Rb 67.82) (1Na 25.8) (1K 49.45) (1Rb 58.31)
    double intensity_SI = 11.0e9 * 1.0e4; // (Li 22.217) (Na 13.435) (K 11.0) (Rb 15.454)
    double e_field_amplitude = sqrt(2.0 * intensity_SI / (speed_of_light_SI * vaccuum_permitivity_SI)) * au_e_field_strength_per_au;
    double fwhm = 169 * au_per_ps; // (Li 191.0) (Na 175.0) (K 169.0) (Rb 179.0)
    double dt = 0.25 * au_per_ps;
    double even_abundance = 5.0; // Triplet = highest first, singlet = lowest first. (5, 3) for all except Rb (7, 5)
    double odd_abundance = 3.0;
    double temperature = 0.37;
    double t_start = -2 * fwhm;
    double t_end = 1250 * au_per_ps;
    size_t n_lw = 15;
    double Delta_B = 200.0 * 1.e6 * 2.0 * M_PI * au_frequency_per_SI; 
    //size_t n_focal = 20;
    //double w_pump = 90.0;
    double w_probe = 25.0;
    size_t n_r = 10, n_theta = 10;
    double w_x = 42.5, w_y = 113.0; // (Li 35 58) (Na 40 90) (K 42.5 113) (Rb 32 90)
    bool custom_pulse_flag = true;
    char *custom_pulse_fname = "Pulses/K_pulse_smooth_taper.dat";
    double homonuclear_abundance = 1.0; //0.933; // The heteronuclear is just 1 - homonuclear_abundance

    // Find the number of steps for during and after the pulse
    double field_duration = 4.0 * fwhm; //40.0 * fwhm;
    size_t n_field = (size_t) ceil(field_duration / dt);
    size_t n_field_free = (size_t) ceil((t_end + t_start) / dt);
    size_t n_steps = n_field + n_field_free;

    rotor_type run_type = PLANAR_ROTOR;
    size_t dim = (run_type == LINEAR_ROTOR) ? (j_max * j_max) : (2 * j_max + 1);
    double E_rot[dim];
    get_rot_energies_2D(j_max, B, E_rot);

    // Set up the molecule parameter struct
    molecule_params m = {
        .B = B, 
        .delta_alpha = delta_alpha,
        .even_abundance = even_abundance,
        .odd_abundance = odd_abundance,
        .homonuclear_abundance = homonuclear_abundance,
        .delta_B = 0.0,
        .E_rot = E_rot,
        .j_max = j_max
    };

    // Set up laser parameters/interaction parameters
    field_params f = {
        .n_pulses = 1,
        .custom_pulse_flag = custom_pulse_flag,
        .e_field_squared = pow(e_field_amplitude, 2),
        .fwhm = fwhm
    };

    // Set up all the parameters in a nice struct
    solver_params p = {
        .dim = dim,
        .field = &f,
        .molecule = &m,
        .temperature = temperature,
        .beta = 1.0/(boltzmann_constant_SI * temperature * au_energy_per_SI),
        .type = run_type,
        .n_steps_field = n_field,
        .n_steps_field_free = n_field_free,
        .dt = dt,
        .t_start = t_start,
        .t_end = t_end
        //.diagonals=&diags
    };

    // Parameters for the various averaging techniques
    avg_params a = {
        .T = temperature,
        .beta = 1.0/(boltzmann_constant_SI * temperature * au_energy_per_SI),
        .n_lw = n_lw,
        .delta_B = Delta_B,
        .n_r = n_r,
        .n_theta = n_theta,
        .w_x = w_x,
        .w_y = w_y,
        .w_probe = w_probe
    };

    // Assign the structs to the pointers (by copying the memory)
    molecule = m;
    field = f;
    params = p;
    avg = a;
}


void get_custom_pulse_spline(char *filename, solver_params params)
{
    // Load cross correlation data into two arrays
    double *pulse_times = (double *) malloc(MAX_CC_LENGTH * sizeof(double));
    double *pulse_envel = (double *) malloc(MAX_CC_LENGTH * sizeof(double));
    size_t n_cross_correlation = 0;
    size_t n_steps = params->n_steps_field + params->n_steps_field_free;
    FILE *pulse_file = fopen(filename, "r");
    while (!feof(pulse_file))
    {
        if (fscanf(pulse_file, "%lf %lf", &pulse_times[n_cross_correlation], &pulse_envel[n_cross_correlation]) == 2)
        {
            printf("t = %f and I = %f\n", pulse_times[n_cross_correlation], pulse_envel[n_cross_correlation]);
            n_cross_correlation++;
        }
    }
    fclose(pulse_file);

    // Then interpolate the data (only withtin the time range of the data)
    gsl_interp_accel *_acc = gsl_interp_accel_alloc();
    gsl_spline *_spline = gsl_spline_alloc(gsl_interp_linear, n_cross_correlation);
    gsl_spline_init(_spline, pulse_times, pulse_envel, n_cross_correlation);

    // Interpolate the data to the time grid of the simulation
    double *cc_spline_times = malloc(n_steps * sizeof(double));
    double *cc_spline_envel = malloc(n_steps * sizeof(double));
    if (!cc_spline_times || !cc_spline_envel) 
    {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    for (size_t i = 0; i < n_steps; i++)
    {
        double t_i = params->t_start + i * params->dt;
        cc_spline_times[i] = t_i;
        if (t_i > pulse_times[0] * au_per_ps && t_i < pulse_times[n_cross_correlation - 1] * au_per_ps)
        {
            printf("t = %f\n", t_i / au_per_ps);
            printf("times %f and %f\n", t_i, pulse_times[0] * au_per_ps); 
            cc_spline_envel[i] = gsl_spline_eval(_spline, t_i / au_per_ps, _acc);
        }
        else
        {
            cc_spline_envel[i] = 0.0;
        }
    }
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_linear, n_steps);
    gsl_spline_init(spline, cc_spline_times, cc_spline_envel, n_steps);
    params->field->spline = spline;
    params->field->acc = acc;
    gsl_spline_free(_spline);
    gsl_interp_accel_free(_acc);
}