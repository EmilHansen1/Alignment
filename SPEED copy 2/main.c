//TODO: Maybe try to run in WSL file system for better performace?
//TODO: Custom pulses + interpolaiton
//TODO: Multiple pulses
//TODO: Implement other rotor types

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "units.h"
#include "utils.h"
#include <omp.h>
#include "planar_rotor.h"
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define TWO_SQRT_TWO_LOG_TWO 2.35482004503094938202
#define BUF_SIZE 1024 //65536
#define MAX_LINE_LENGTH 1024
#define MAX_COLS 1024
#define MAX_CC_LENGTH 200

typedef double _Complex dcmplx;

// Prototypes
void get_focal_average_weights(size_t n_focal, const double I0, const double w_pump, const double w_probe, double weights[n_focal], double intensities[n_focal]);
void get_rot_const_avg_weights(const size_t n, const double Delta_B, const double B, double weights[2*n+1], double B_values[2*n+1]);
double e_field_from_intensity(const double intensity);
int count_lines(FILE* file);

int main(int argc, char *argv[])
{
    // Constants
    size_t j_max = 25;
    size_t dim = 2*j_max + 1; 
    double B = 2.0 * M_PI * 1.33 * 1e9 * au_frequency_per_SI;
    double delta_alpha = 31.2 * pow(au_per_angstrom, 3.0);
    double intensity_SI =  15.6e9 * 1e4;
    double e_field_amplitude = sqrt(2.0 * intensity_SI / (speed_of_light_SI * vaccuum_permitivity_SI)) * au_e_field_strength_per_au;
    double fwhm = 171. * au_per_ps;
    double dt = 0.5 * au_per_ps;
    double even_abundance = 5.0;
    double odd_abundance = 3.0;
    double temperature = 0.37;
    double t_start = -2 * fwhm;
    double t_end = 1250 * au_per_ps;
    size_t n_lw = 10;
    double Delta_B = 200.0 * 1.e6 * 2.0 * M_PI * au_frequency_per_SI;
    size_t n_focal = 15;
    double w_pump = 70;
    double w_probe = 26;
    bool custom_pulse_flag = true;
    char *custom_pulse_fname = "Pulses/Na_pulse.dat";
    
    // Find the number of steps for during and after the pulse
    double field_duration = 4.0 * fwhm;
    size_t n_field = (size_t) ceil(field_duration / dt);
    size_t n_field_free = (size_t) ceil((t_end + t_start) / dt);
    size_t n_steps = n_field + n_field_free;

    // The type of rotor used in the simulations
    rotor_type run_type = PLANAR_ROTOR;

    // Create rotational energy list
    double E_rot[dim];
    switch(run_type)
    {
    case PLANAR_ROTOR:
        get_rot_energies_2D(j_max, B, E_rot);
        break;

    case LINEAR_ROTOR:
        break;

    case SYMMETRIC_TOP:
        break;

    default:
        break;
    }

    // Set up the molecule parameter struct
    molecule_params mol_params = {
        .B = B, 
        .delta_alpha = delta_alpha,
        .even_abundance = even_abundance,
        .odd_abundance = odd_abundance,
        .delta_B = 0.0,
        .E_rot = E_rot,
        .j_max = j_max
    };

    // Set up laser parameters/interaction parameters
    field_params laser_params = {
        .n_pulses = 1,
        .custom_pulse_flag = custom_pulse_flag,
        .e_field_squared = pow(e_field_amplitude, 2),
        .fwhm = fwhm
    };

    // Set up all the parameters in a nice struct
    solver_params params = {
        .dim = 2 * j_max + 1,
        .field = &laser_params,
        .molecule = &mol_params,
        .temperature = temperature,
        .beta = 1.0/(boltzmann_constant_SI * temperature * au_energy_per_SI),
        .type = run_type,
        .n_steps_field = n_field,
        .n_steps_field_free = n_field_free,
        .dt = dt,
        .t_start = t_start,
        .t_end = t_end
    };

    // Load cross correlation and interpolate
    double *pulse_times = (double *) malloc(MAX_CC_LENGTH * sizeof(double));
    double *pulse_envel = (double *) malloc(MAX_CC_LENGTH * sizeof(double));
    size_t n_cross_correlation = 0;
    FILE *pulse_file = fopen(custom_pulse_fname, "r");
    while(!feof(pulse_file))
    {
        if(fscanf(pulse_file, "%lf %lf", &pulse_times[n_cross_correlation], &pulse_envel[n_cross_correlation]) == 2)
            n_cross_correlation++;
    }
    fclose(pulse_file);
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_linear, n_cross_correlation);
    gsl_spline_init(spline, pulse_times, pulse_envel, n_cross_correlation);

    // Containers for times and cos² expectation values
    // NOTE: This chuck HAS to come after the loading of the pulse otherwise
    // NOTE: the program breaks down, and the output goes to absolute shit
    double times[n_steps];
    double cos2_exp[n_steps];
    for(int i = 0; i < n_steps; i++) 
    {
        times[i] = (t_start + i * dt) / au_per_ps;
        cos2_exp[i] = 0.;
    }

    // Calculate thermal weights
    double weights[dim]; 
    size_t n_states_thermal = get_planar_thermal_weights(&params, weights);
    n_states_thermal = (size_t) ceil(n_states_thermal / 2.);

    // Calculate lineshape average weights 
    double lw_weights[2 * n_lw + 1];
    double B_values[2 * n_lw + 1];
    get_rot_const_avg_weights(n_lw, Delta_B, B, lw_weights, B_values);
    
    // Calculate focal volume average weights
    double fva_weights[n_focal];
    double intensities[n_focal];
    get_focal_average_weights(n_focal, intensity_SI, w_pump, w_probe, fva_weights, intensities);
    for(int i = 0; i < n_focal; i++)
    {
        printf("Weight: %f\n", fva_weights[i]);
        printf("Intensity: %f\n", intensities[i]/(1e9 * 1e4));
    }
        

    // Perform averaging
    double timer_start = clock();
    size_t n_simulations = (2 * n_lw + 1) * n_states_thermal * n_focal;
    size_t counter = 0;
    double process_bar = 0.0, run_weight = 0.0;
    dcmplx psi[dim];
    printf("Number of jobs: %li\n", n_simulations);
    for(size_t m = 0; m < n_focal; m++)
    {
        for(size_t k = 0; k < 2*n_lw+1; k++)
        {
            for(size_t i = 0; i < n_states_thermal; i++)
            {
                switch(run_type)
                {
                case PLANAR_ROTOR:
                    // Initial state
                    for(size_t j = 0; j < dim; j++) psi[j] = 0.0 + 0.0 * I;
                    psi[j_max + i] = 1.0 + 0.0 * I;

                    // Change values for the current run
                    params.molecule->B = B_values[k];
                    params.field->e_field_squared = pow(e_field_from_intensity(intensities[m]), 2);
                    printf("%f\n", params.field->e_field_squared);
                    run_weight =  weights[j_max + i] * lw_weights[k]* fva_weights[m];

                    // Propagate single state
                    double t1 = clock();
                    planar_rotor_propagation(&params, psi, cos2_exp, run_weight);
                    double t2 = clock();

                    // Log process
                    counter++;
                    process_bar = (double) counter / n_simulations;
                    printf("J = %li completed in %f s\t\t %.1f%%  \n", i, (t2 - t1) / CLOCKS_PER_SEC, process_bar*100);
                    break;

                case LINEAR_ROTOR:
                    break;

                case SYMMETRIC_TOP:
                    break;

                default:
                    break;
                }
                
            }
        }
    }
    double timer_end = clock();
    printf("\n All %li simulations completed in %f seconds.\n\n", n_simulations, (timer_end - timer_start) / CLOCKS_PER_SEC);

    // Write out the degree of alignment to file
    FILE *cos2_file = fopen("cos2.out", "w");
    for(int i = 0; i < n_steps; i++)
    {
        double field_intensity = 0.0;
        if(pulse_times[0] <= times[i] && times[i] <= pulse_times[n_cross_correlation-1])
            field_intensity = gsl_spline_eval(spline, times[i], acc);
        fprintf(cos2_file, "%f, %f, %f\n", times[i], cos2_exp[i], field_intensity);
    }
    fclose(cos2_file);

    // Free the GSL data
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);

    return EXIT_SUCCESS;
}


// ------------------------------------------------------------------------------------------------------------------------------- //

double gaussian_profile(const double r, const double I0, const double w)
{
    return I0 * exp(-2.0 * pow(r / w, 2));
}

void get_focal_average_weights(size_t n_focal, const double I0, const double w_pump, const double w_probe, double weights[n_focal], double intensities[n_focal])
{
    double r_max = 1.5 * w_probe;           // Maximum radius (Areg used 1.7? But 1.5 is more than 99% of the power)
    double dr = (double) r_max / n_focal;   // The radial step size
    double radii[n_focal];                  // The radii used form [dr; 1.5w+dr]
    double normalization = 0.0;             // The total sum of all unnormalized weights

    if(n_focal == 1)
    {
        intensities[0] = I0;
        weights[0] = 1.0;
    }
    else
    {
        for(size_t i = 0; i < n_focal; i++)
        {
            // Add 1 to avoid weight[0] = 0
            radii[i] = (double) i * dr + 0.1;
            intensities[i] = gaussian_profile(radii[i], I0, w_pump);
            weights[i] = gaussian_profile(radii[i], 1.0, w_probe);
            normalization += weights[i];
        }

        // Normalize the weights so they sum to unit
        for(size_t i = 0; i < n_focal; i++) 
            weights[i] /= normalization;
    }

}


double gaussian_lineshape(const double B, const double B0, const double Delta_B)
{
    double sigma_B = Delta_B / TWO_SQRT_TWO_LOG_TWO;
    return exp(-0.5 * pow((B - B0) / sigma_B, 2)) / (sqrt(2.0 * M_PI) * sigma_B);
}

void get_rot_const_avg_weights(const size_t n, const double Delta_B, const double B, double weights[2*n+1], double B_values[2*n+1])
{
    // Find the B-values
    size_t n_B = 2 * n + 1;
    double delta = 3.0 * Delta_B / ((double) n_B);
    for(size_t i = 0; i < n_B; i++)
        B_values[i] = B - 3.0 * Delta_B / 2.0 + delta * i;

    // Find the weights
    for(size_t i = 0; i < n_B; i++)
        weights[i] = gaussian_lineshape(B_values[i], B, Delta_B);
    
    // Then normalize so they sum to unity
    double weight_sum = 0.0;
    for(size_t i = 0; i < n_B; i++) weight_sum += weights[i];
    for(size_t i = 0; i < n_B; i++) weights[i] /= weight_sum;
}

double e_field_from_intensity(const double intensity)
{
    return sqrt(2.0 * intensity / (speed_of_light_SI * vaccuum_permitivity_SI)) * au_e_field_strength_per_au;
}


int count_lines(FILE* file)
{
    char buf[BUF_SIZE];
    int counter = 0;
    for(;;)
    {
        size_t res = fread(buf, 1, BUF_SIZE, file);
        if (ferror(file))
            return -1;

        int i;
        for(i = 0; i < res; i++)
            if (buf[i] == '\n')
                counter++;

        if (feof(file))
            break;
    }

    rewind(file);

    return counter;
}