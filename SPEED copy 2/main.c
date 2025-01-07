//TODO: Maybe try to run in WSL file system for better performace?
//TODO: Custom pulses + interpolaiton
//TODO: Multiple pulses
//TODO: Implement other rotor types
//TODO: Add imaginary value to energy instead of averaging

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <stdlib.h>
#include <stdbool.h>
#include <omp.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include "units.h"
#include "utils.h"
#include "planar_rotor.h"
#include "averaging.h"
#include "input_loader.h"

#define MAX_THERMAL_STATES 1024

typedef double _Complex dcmplx;

// Prototypes
double e_field_from_intensity(const double intensity);


int main(int argc, char *argv[])
{
    omp_set_num_threads(1);

    // --- Load the input parameters --- //


    solver_params params = load_input(argc, argv);
    print_solver_params(&params);
    get_rot_energies_2D(params.molecule.j_max, params.molecule.B, params.molecule.E_rot);


    // --- Prepare observables and time --- //


    size_t n_steps = params.n_steps_field + params.n_steps_field_free;
    double times[n_steps];
    double cos2_exp[n_steps];
    for (int i = 0; i < n_steps; i++) 
    {
        // Set the time and zero-initialize the degree of alignment
        times[i] = (params.t_start + i * params.dt) / au_per_ps;
        cos2_exp[i] = 0.0;
    }

    
    // --- Calculate thermal weights --- //


    double *weights = (double *) calloc(MAX_THERMAL_STATES, sizeof(double *)); 
    size_t n_states_thermal = get_planar_thermal_weights(&params, weights);
    printf("n_states_thermal: %li\n", n_states_thermal);
    const size_t n_homo = n_states_thermal;
    if(params.molecule.homonuclear_abundance != 1.0)
    {
        // For heteronuclear molecules, there is no even-odd asymmetry
        params.molecule.even_abundance = params.molecule.odd_abundance = 1.0;
        n_states_thermal += get_planar_thermal_weights(&params, weights + n_states_thermal);

        // Renormalize weights
        double total_weight = 0.0;
        for(size_t i = 0; i < n_states_thermal; i++) total_weight += weights[i];
        
    }
    printf("n_states_thermal: %li\n", n_states_thermal);


    // --- Calculate lineshape average weights  --- //


    double lw_weights[2 * params.avg.n_lw + 1];
    double B_values[2 * params.avg.n_lw + 1];
    get_rot_const_avg_weights(params.avg.n_lw, params.avg.delta_B, params.molecule.B, lw_weights, B_values);
    
    
    // --- Calculate focal volume average weights --- //


    size_t n_focal = params.avg.n_r * params.avg.n_theta;
    double ellip_weights[n_focal]; 
    double ellip_int[n_focal];
    get_elliptical_focal_average_weights(params.avg.n_r, params.avg.n_theta, params.field.intensity_SI, 
                                         params.avg.w_x, params.avg.w_y, params.avg.w_probe, ellip_weights, ellip_int);

    
    // --- Perform averaging --- //


    double timer_start = omp_get_wtime();
    size_t n_simulations = (2 * params.avg.n_lw + 1) * n_states_thermal * n_focal;
    size_t counter = 0;
    double process_bar = 0.0, run_weight = 0.0;
    dcmplx psi[params.dim];
    double t1, t2;
    printf("Number of jobs: %li\n", n_simulations);
    
    #pragma omp parallel for collapse(2) private(psi, run_weight, t1, t2) shared(counter, process_bar, cos2_exp)
    for (size_t m = 0; m < n_focal; m++)
    {
        for (size_t k = 0; k < 2 * params.avg.n_lw + 1; k++)
        {
            for (size_t i = 0; i < n_states_thermal; i++)
            {
                // Change values for the current run 
                params.molecule.B = B_values[k];
                params.field.e_field_squared = pow(e_field_from_intensity(ellip_int[m]), 2);
                run_weight =  weights[params.molecule.j_max + i] * lw_weights[k] * ellip_weights[m] * (
                    (i < n_homo) ? params.molecule.homonuclear_abundance : 1.0 - params.molecule.homonuclear_abundance
                );

                switch (params.type)
                {
                    case PLANAR_ROTOR:

                        // Initial state
                        for (size_t j = 0; j < params.dim; j++) psi[j] = 0.0 + 0.0 * I;
                        psi[params.molecule.j_max + i] = 1.0 + 0.0 * I;
                        get_rot_energies_2D(params.molecule.j_max, params.molecule.B, params.molecule.E_rot);

                        // Propagate single state
                        t1 = omp_get_wtime();
                        planar_rotor_propagation(&params, psi, cos2_exp, run_weight);
                        t2 = omp_get_wtime();

                        // Log process
                        #pragma omp critical
                        {
                            counter++;
                            process_bar = (double) counter / n_simulations;
                            printf("J = %i completed in %.5f s\t\t %.1f%%  \n", i, (t2 - t1), process_bar*100);
                        }

                        break;

                    case LINEAR_ROTOR:
                        // Initial state
                        for (size_t j = 0; j < params.dim; j++) psi[j] = 0.0 + 0.0 * I;
                        psi[i] = 1.0 + 0.0 * I;

                        // Propagate single state
                        t1 = clock();
                        //linear_rotor_propagation(&params, psi, cos2_exp, run_weight);
                        t2 = clock();

                        // Log process
                        counter++;
                        process_bar = (double) counter / n_simulations;
                        printf("J = %i completed in %f s\t\t %.1f%%  \n", i, (t2 - t1) / CLOCKS_PER_SEC, process_bar*100);
                        break;

                    case SYMMETRIC_TOP:
                        break;

                    default:
                        break;
                }
            }
        }
    }
    double timer_end = omp_get_wtime();
    printf("\nAll %li simulations completed in %f seconds.\n\n", n_simulations, (timer_end - timer_start));


    // Write out the degree of alignment to file
    FILE *cos2_file = fopen("cos2.out", "w");
    for (int i = 0; i < n_steps; i++)
    {
        double field_intensity = 0.0;
        if (times[i] > -2.0 * params.field.fwhm / au_per_ps && times[i] < 2.0 * params.field.fwhm / au_per_ps)
        {
            field_intensity = params.field.custom_pulse_flag 
                ? gsl_spline_eval(params.field.spline, times[i] * au_per_ps, params.field.acc) 
                : e_field_squared(times[i] * au_per_ps, 1.0, params.field.fwhm);
        }
        fprintf(cos2_file, "%f, %f, %f\n", times[i], cos2_exp[i], field_intensity);
    }
    fclose(cos2_file);


    // Free allocated memory
    gsl_spline_free(params.field.spline);
    gsl_interp_accel_free(params.field.acc);
    free(weights);
    free(params.molecule.E_rot);

    return EXIT_SUCCESS;
}