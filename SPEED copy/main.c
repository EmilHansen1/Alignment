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

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef double _Complex dcmplx;



int main(int argc, char *argv[])
{
    // Constants
    size_t j_max = 140;
    size_t dim = 2*j_max + 1; 
    double B = 2.0 * M_PI * 1.33 * 1e9 * au_frequency_per_SI;
    double delta_alpha = 31.2 * pow(au_per_angstrom, 3.0); // 31.2
    double intensity_SI =  50.e9 * 1e4; //1e12 * 1e4;  // W/cm²
    double e_field_amplitude = sqrt(2.0 * intensity_SI / (speed_of_light_SI * vaccuum_permitivity_SI)) * au_e_field_strength_per_au;
    double fwhm = 1. * au_per_ps; //300 * au_per_fs; 
    size_t n_steps_field = 150; //200;
    size_t n_steps_field_free = 5000;
    double dt_field_free = 0.25 * au_per_ps;
    double even_abundance = 3.0;
    double odd_abundance = 5.0;
    bool debug_temp = true;

    // Create rotational energy list
    double E_rot[dim];
    printf("E= %f", E_rot[10]);
    get_rot_energies_2D(j_max, B, E_rot);

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

    // The type of rotor used in the simulations
    rotor_type run_type = PLANAR_ROTOR;

    // Set up laser parameters/interaction parameters
    field_params laser_params = {
        .n_pulses = 1,
        .custom_pulse_flag = false,
        .e_field_squared = &e_field_squared,
        .fwhm = &fwhm
    };

    // Set up all the parameters in a nice struct
    solver_params params = {
        .dim = 2 * j_max + 1,
        .field = &laser_params,
        .molecule = &mol_params,
        .temperature = 0.37,
        .type = run_type,
        .n_steps_field = n_steps_field,
        .n_steps_field_free = n_steps_field_free
    };

    double temperature = 0.37;
    double beta = 1.0/(boltzmann_constant_SI * temperature * au_energy_per_SI);

    // Containers for times and cos² expectation values
    double field_times[n_steps_field];
    double field_free_times[n_steps_field_free];
    double cos2_exp_field[n_steps_field];
    double cos2_exp_field_free[n_steps_field_free];

    // TODO: calculate cos2 before propagation, and then simply just calculate AFTER each step
    // in each solver (with field and field free)  

    double temp_start_time = (double) clock();
    if(!debug_temp)
    {
        // Single j-state propagation
        double cos2_field[n_steps_field];
        double cos2_field_free[n_steps_field_free];

        // Create coefficient matrix
        dcmplx psi[dim];
        psi[j_max] = 1.0 + 0.0 * I;

        // Propagate from -2fwhm to 2fwhm
        double start_time = (double) clock();
        field_times[0] = -2 * fwhm;
        field_propagation(j_max, n_steps_field, (double) 4.0*fwhm/n_steps_field, B, fwhm, pow(e_field_amplitude, 2.0), 
                        delta_alpha, field_times, E_rot, psi, cos2_field);
        double end_time = (double) clock();
        printf("\nODE solved in %f seconds\n", (end_time - start_time)/CLOCKS_PER_SEC);

        // Field-free propagation
        start_time = (double) clock();
        field_free_times[0] = field_times[n_steps_field - 1];
        field_free_propagation(j_max, dt_field_free, n_steps_field_free, field_free_times, B, psi, cos2_field_free);
        end_time = (double) clock();
        printf("\nFF solved in %f seconds\n", (end_time - start_time)/CLOCKS_PER_SEC);

        // Write cos2 expectation values to file
        FILE* cos2_file = fopen("cos2.out", "w");
        for(size_t i = 0; i < n_steps_field; i++) fprintf(cos2_file, "%f %f\n", field_times[i]/au_per_ps, cos2_field[i]);
        for(size_t i = 0; i < n_steps_field_free; i++) fprintf(cos2_file, "%f %f\n", field_free_times[i]/au_per_ps, cos2_field_free[i]);
        printf("I'm done!\n");
        
        double weights[params.dim];
        size_t n_states = get_thermal_weights_funciton(&params, weights);
        printf("N states: %d\n", n_states);
        return 0;
    }


    // Set up for thermal averaging
    size_t n_states_thermal = 0;
    double thermal_weights[dim];
    double partition_function = 0;
    if(temperature != 0.0)
    {
        int j = -j_max;
        for(size_t i = 0; i < dim; i++)
        {
            thermal_weights[i] = ((j % 2) ? even_abundance : odd_abundance) * exp(-beta * E_rot[i]);
            partition_function += thermal_weights[i];
            j++;
        }

        double thermal_contribution = 1;
        size_t counter = j_max;
        while(thermal_contribution > 1e-5)
        {
            thermal_contribution = thermal_weights[counter];
            counter++;
        }   
        n_states_thermal = counter - j_max;
    }

    printf("THERMAL STATES: %i\n", n_states_thermal);
    // Start doing the stuffs
    size_t counter = 0;
    double cos2_field[2*n_states_thermal+1][n_steps_field];
    double cos2_field_free[2*n_states_thermal+1][n_steps_field_free];
    for(int j = j_max - n_states_thermal; j <= j_max + n_states_thermal; j++)
    {
        counter++;
        // Create coefficient matrix
        dcmplx psi[dim];
        for(size_t i = 0; i < dim; i++)
                psi[i] = (i == j) ? 1.0 : 0.0;

        // Propagate from -2fwhm to 2fwhm
        double start_time = (double) clock();
        field_times[0] = -3 * fwhm;
        field_propagation(j_max, n_steps_field, (double) 6.0*fwhm/((double)n_steps_field), B, fwhm, pow(e_field_amplitude, 2.0), 
                        delta_alpha, field_times, E_rot, psi, cos2_field[counter]);
        double end_time = (double) clock();
        printf("\nODE solved in %f seconds\n", (end_time - start_time)/CLOCKS_PER_SEC);
        for(int i = 0; i < n_steps_field; i++) cos2_field[counter][i] *= thermal_weights[j]/partition_function;

        // Field-free propagation
        start_time = (double) clock();
        field_free_times[0] = field_times[n_steps_field - 1];
        field_free_propagation(j_max, dt_field_free, n_steps_field_free, field_free_times, B, psi, cos2_field_free[counter]);
        end_time = (double) clock();
        printf("\nFF solved in %f seconds\n", (end_time - start_time)/CLOCKS_PER_SEC);
        for(int i = 0; i < n_steps_field_free; i++) cos2_field_free[counter][i] *= thermal_weights[j]/partition_function;
    }

    double temp_end_time = (double) clock();
    printf("\nALL SOLVED IN in %f SECONDS\n", (temp_end_time - temp_start_time)/CLOCKS_PER_SEC);


    // Write cos2 expectation values to file
    FILE* cos2_file = fopen("cos2.out", "w");
    for(size_t i = 0; i < n_steps_field; i++)
    {
        for(int j = 0; j < 2*n_states_thermal+1; j++)
            cos2_exp_field[i] += cos2_field[j][i];
        fprintf(cos2_file, "%f %f\n", field_times[i]/au_per_ps, cos2_exp_field[i]);
    }
    for(size_t i = 0; i < n_steps_field_free; i++)
    {
        for(int j = 0; j < 2*n_states_thermal+1; j++)
            cos2_exp_field_free[i] += cos2_field_free[j][i];
        fprintf(cos2_file, "%f %f\n", field_free_times[i]/au_per_ps, cos2_exp_field_free[i]);
        
    }

    return 0;
}