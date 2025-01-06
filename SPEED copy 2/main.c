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
#include <time.h>
#include <math.h>
#include <string.h>
#include "units.h"
#include "utils.h"
#include <omp.h>
#include "planar_rotor.h"
#include "input_loader.h"
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>


#define TWO_SQRT_TWO_LOG_TWO 2.35482004503094938202
#define BUF_SIZE 1024 //65536
#define MAX_LINE_LENGTH 1024
#define MAX_COLS 1024
#define MAX_CC_LENGTH 200
#define EPSILON 1.0e-9
#define MAX_THERMAL_STATES 1024

typedef double _Complex dcmplx;

// Prototypes
void get_focal_average_weights(size_t n_focal, const double I0, const double w_pump, const double w_probe, double weights[n_focal], double intensities[n_focal]);
void get_rot_const_avg_weights(const size_t n, const double Delta_B, const double B, double weights[2*n+1], double B_values[2*n+1]);
void get_elliptical_focal_average_weights(const size_t n_r, const size_t n_theta, const double I0, const double w_x, const double w_y, 
                                            const double w_probe, double weights[n_r*n_theta], double intensities[n_r*n_theta]);
double e_field_from_intensity(const double intensity);
int count_lines(FILE* file);

int main(int argc, char *argv[])
{

    printf("argc: %d\n", argc);

    // Print each argument in argv
    printf("argv:\n");
    for (int i = 0; i < argc; i++) 
    {
        printf("\targv[%d]: %s\n", i, argv[i]);
    }

    // --- Load the input parameters --- //


    molecule_params molecule;
    field_params field;
    solver_params params;
    avg_params avg;
    load_input(argc, argv, &molecule, &field, &params, &avg);


    // --- Prepare observables and time --- //


    size_t n_steps = params.n_steps_field + params.n_steps_field_free;
    double times[n_steps];
    double cos2_exp[n_steps];
    for (int i = 0; i < n_steps; i++) 
    {
        times[i] = (params.t_start + i * params.dt) / au_per_ps;
        cos2_exp[i] = 0.;
    }

    
    // --- Calculate thermal weights --- //

    printf("%li\n", params.temperature);
    return 0;

    double *weights = (double *) calloc(MAX_THERMAL_STATES, sizeof(double *)); 
    size_t n_states_thermal = get_planar_thermal_weights(&params, weights);
    const size_t n_homo = n_states_thermal;
    if(molecule.homonuclear_abundance != 1.0)
    {
        params.molecule->even_abundance = params.molecule->odd_abundance = 1.0;
        n_states_thermal += get_planar_thermal_weights(&params, weights + n_states_thermal);
    }
    //n_states_thermal = (size_t) ceil(n_states_thermal / 2.);
    printf("n_states_thermal: %li\n", n_states_thermal);


    // --- Calculate lineshape average weights  --- //


    double lw_weights[2 * avg.n_lw + 1];
    double B_values[2 * avg.n_lw + 1];
    get_rot_const_avg_weights(avg.n_lw, avg.delta_B, molecule.B, lw_weights, B_values);
    

    // --- Calculate focal volume average weights --- //


    size_t n_focal = avg.n_r * avg.n_theta;
    double ellip_weights[n_focal]; 
    double ellip_int[n_focal];
    get_elliptical_focal_average_weights(avg.n_r, avg.n_theta, field.e_field_squared, avg.w_x, avg.w_y, avg.w_probe, ellip_weights, ellip_int);

    /*
    // Constants
    size_t j_max = 100;  
    double B = 2.0 * M_PI * 0.69 * 1e9 * au_frequency_per_SI; // (3Li 5.25) (3Na 1.33) (3K 0.69) (3Rb 0.28) (1Na 4.60) (1.78 1K) (0.69 1Rb)
    double delta_alpha = 68.54 * pow(au_per_angstrom, 3.0); // (3Li 65.73) (3Na 31.06) (3K 68.54) (3Rb 67.82) (1Na 25.8) (1K 49.45) (1Rb 58.31)
    double intensity_SI = 25.0e9 * 1.0e4; // (Li 22.217) (Na 13.435) (K 11.0) (Rb 15.454)
    double e_field_amplitude = sqrt(2.0 * intensity_SI / (speed_of_light_SI * vaccuum_permitivity_SI)) * au_e_field_strength_per_au;
    double fwhm = 0.75 * au_per_ps; // (Li 191.0) (Na 175.0) (K 169.0) (Rb 179.0)
    double dt = 0.1 * au_per_ps;
    double even_abundance = 5.0; // Triplet = highest first, singlet = lowest first. (5, 3) for all except Rb (7, 5)
    double odd_abundance = 3.0;
    double temperature = 0.37;
    double t_start = -2 * fwhm;
    double t_end = 1250 * au_per_ps;
    size_t n_lw = 0;
    double Delta_B = 50.0 * 1.e6 * 2.0 * M_PI * au_frequency_per_SI; 
    size_t n_focal = 20;
    double w_pump = 90.0;
    double w_probe = 25.0;
    size_t n_r = 50, n_theta = 1;
    double w_x = 70.0, w_y = 70.0; // (Li 35 58) (Na 40 90) (K 42.5 113) (Rb 32 90)
    bool custom_pulse_flag = false;
    char *custom_pulse_fname = "Pulses/K_pulse_smooth_taper.dat";
    double homonuclear_abundance = 1.0; //0.933; // The heteronuclear is just 1 - homonuclear_abundance
    
    // Find the number of steps for during and after the pulse
    double field_duration = 4.0 * fwhm; //40.0 * fwhm;
    size_t n_field = (size_t) ceil(field_duration / dt);
    size_t n_field_free = (size_t) ceil((t_end + t_start) / dt);
    size_t n_steps = n_field + n_field_free;

    // The type of rotor used in the simulations
    rotor_type run_type = PLANAR_ROTOR;
    size_t dim = (run_type == LINEAR_ROTOR) ? j_max * j_max : 2 * j_max + 1;

    // Linear stuffs
    free_3d_qn linear_qn[j_max * j_max];
    double cos2_diag[dim];
    double cos2_off_diag[dim - 2];
    tridiag diags = {.diag=cos2_diag, .off_diag=cos2_off_diag};
    // Create rotational energy list
    double E_rot[dim];
    switch (run_type)
    {
    case PLANAR_ROTOR:
        get_rot_energies_2D(j_max, B, E_rot);
        break;

    case LINEAR_ROTOR:
        //create_idx_map(j_max, linear_qn);
        //fill_diagonals(j_max, linear_qn, &diags);
        //get_rot_energies_3D(j_max, B, E_rot, linear_qn);
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
        .dim = dim,
        .field = &laser_params,
        .molecule = &mol_params,
        .temperature = temperature,
        .beta = 1.0/(boltzmann_constant_SI * temperature * au_energy_per_SI),
        .type = run_type,
        .n_steps_field = n_field,
        .n_steps_field_free = n_field_free,
        .dt = dt,
        .t_start = t_start,
        .t_end = t_end,
        .diagonals=&diags
    };

    // Load cross correlation and interpolate
    double *pulse_times = (double *) malloc(MAX_CC_LENGTH * sizeof(double));
    double *pulse_envel = (double *) malloc(MAX_CC_LENGTH * sizeof(double));
    size_t n_cross_correlation = 0;
    FILE *pulse_file = fopen(custom_pulse_fname, "r");
    while (!feof(pulse_file))
    {
        if (fscanf(pulse_file, "%lf %lf", &pulse_times[n_cross_correlation], &pulse_envel[n_cross_correlation]) == 2)
        {
            //printf("t = %f and I = %f\n", pulse_times[n_cross_correlation], pulse_envel[n_cross_correlation]);
            n_cross_correlation++;
        }
    }
    fclose(pulse_file);
    gsl_interp_accel *_acc = gsl_interp_accel_alloc();
    gsl_spline *_spline = gsl_spline_alloc(gsl_interp_linear, n_cross_correlation);
    gsl_spline_init(_spline, pulse_times, pulse_envel, n_cross_correlation);
    double cc_spline_times[n_steps];
    double cc_spline_envel[n_steps];
    for (size_t i = 0; i < n_steps; i++)
    {
        double t_i = t_start + i * dt;
        cc_spline_times[i] = t_i;
        if (t_i > pulse_times[0] * au_per_ps && t_i < pulse_times[n_cross_correlation - 1] * au_per_ps)
        {
            //printf("t = %f\n", t_i / au_per_ps);
            //printf("times %f and %f\n", t_i, pulse_times[0] * au_per_ps); 
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
    laser_params.spline = spline;
    laser_params.acc = acc;

    // Containers for times and cos² expectation values
    // NOTE: This chunk HAS to come after the loading of the pulse otherwise
    // NOTE: the program breaks down, and the output goes to absolute shit
    double times[n_steps];
    double cos2_exp[n_steps];
    for (int i = 0; i < n_steps; i++) 
    {
        times[i] = (t_start + i * dt) / au_per_ps;
        cos2_exp[i] = 0.;
    }

    // Calculate thermal weights
    double *weights = (double *) calloc(MAX_THERMAL_STATES, sizeof(double *)); 
    size_t n_states_thermal = get_planar_thermal_weights(&params, weights);
    const size_t n_homo = n_states_thermal;
    if(homonuclear_abundance != 1.0)
    {
        params.molecule->even_abundance = params.molecule->odd_abundance = 1.0;
        n_states_thermal += get_planar_thermal_weights(&params, weights + n_states_thermal);
    }
    //n_states_thermal = (size_t) ceil(n_states_thermal / 2.);
    printf("n_states_thermal: %li\n", n_states_thermal);

    // Calculate lineshape average weights 
    double lw_weights[2 * n_lw + 1];
    double B_values[2 * n_lw + 1];
    get_rot_const_avg_weights(n_lw, Delta_B, B, lw_weights, B_values);
    
    // Calculate focal volume average weights
    double fva_weights[n_focal];
    double intensities[n_focal];
    get_focal_average_weights(n_focal, intensity_SI, w_pump, w_probe, fva_weights, intensities);
    
    // Calculate focal volume average weights for the elliptical pump
    double ellip_weights[n_r*n_theta]; 
    double ellip_int[n_r*n_theta];
    get_elliptical_focal_average_weights(n_r, n_theta, intensity_SI, w_x, w_y, w_probe, ellip_weights, ellip_int);
    n_focal = n_r * n_theta;


    // Perform averaging
    double timer_start = clock();
    size_t n_simulations = (2 * n_lw + 1) * n_states_thermal * n_focal;
    size_t counter = 0;
    double process_bar = 0.0, run_weight = 0.0;
    dcmplx psi[dim];
    double t1, t2;
    printf("Number of jobs: %li\n", n_simulations);
    for (size_t m = 0; m < n_focal; m++)
    {
        for (size_t k = 0; k < 2*n_lw+1; k++)
        {
            for (size_t i = 0; i < n_states_thermal; i++)
            {
                // Change values for the current run
                params.molecule->B = B_values[k];
                params.field->e_field_squared = pow(e_field_from_intensity(ellip_int[m]), 2);
                run_weight =  weights[j_max + i] * lw_weights[k] * ellip_weights[m] * ((i < n_homo) ? homonuclear_abundance : 1.0 - homonuclear_abundance);

                switch (run_type)
                {
                    case PLANAR_ROTOR:
                        // Initial state
                        for(size_t j = 0; j < dim; j++) psi[j] = 0.0 + 0.0 * I;
                        psi[j_max + i] = 1.0 + 0.0 * I;

                        // Propagate single state
                        t1 = clock();
                        planar_rotor_propagation(&params, psi, cos2_exp, run_weight);
                        t2 = clock();

                        // Log process
                        counter++;
                        process_bar = (double) counter / n_simulations;
                        printf("J = %i completed in %.5f s\t\t %.1f%%  \n", i, (t2 - t1) / CLOCKS_PER_SEC, process_bar*100);
                        break;

                    case LINEAR_ROTOR:
                        // Initial state
                        for(size_t j = 0; j < dim; j++) psi[j] = 0.0 + 0.0 * I;
                        psi[i] = 1.0 + 0.0 * I;

                        // Propagate single state
                        t1 = clock();
                        linear_rotor_propagation(&params, psi, cos2_exp, run_weight);
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
    double timer_end = clock();
    printf("\nAll %li simulations completed in %f seconds.\n\n", n_simulations, (timer_end - timer_start) / CLOCKS_PER_SEC);


    // Write out the degree of alignment to file
    FILE *cos2_file = fopen("cos2.out", "w");
    for (int i = 0; i < n_steps; i++)
    {
        double field_intensity = 0.0;
        if (times[i] > -2.0 * fwhm / au_per_ps && times[i] < 2.0 * fwhm / au_per_ps)
            field_intensity = custom_pulse_flag ? gsl_spline_eval(spline, times[i] * au_per_ps, acc) 
                                                : e_field_squared(times[i] * au_per_ps, 1.0, fwhm);
        fprintf(cos2_file, "%f, %f, %f\n", times[i], cos2_exp[i], field_intensity);
    }
    fclose(cos2_file);


    // Free the GSL data
    free(pulse_times);
    free(pulse_envel);
    gsl_spline_free(_spline);
    gsl_interp_accel_free(_acc);
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    */
    return EXIT_SUCCESS;
}


// ------------------------------------------------------------------------------------------------------------------------------- //

double gaussian_profile(const double r, const double I0, const double w)
{
    return I0 * exp(-2.0 * pow(r / w, 2));
}

void get_focal_average_weights(size_t n_focal, const double I0, const double w_pump, const double w_probe, double weights[n_focal], double intensities[n_focal])
{
    double r_max = 1.7 * w_probe;           // Maximum radius (Areg used 1.7? But 1.5 is more than 99% of the power)
    double dr = (double) r_max / n_focal;   // The radial step size
    double radii[n_focal];                  // The radii used form [dr; 1.5w+dr]
    double normalization = 0.0;             // The total sum of all unnormalized weights

    if (n_focal == 1)
    {
        intensities[0] = I0;
        weights[0] = 1.0;
    }
    else
    {
        for (size_t i = 0; i < n_focal; i++)
        {
            // Add 1 to avoid weight[0] = 0
            radii[i] = (double) i * dr + 0.1;
            intensities[i] = gaussian_profile(radii[i], I0, w_pump);
            weights[i] = radii[i] * gaussian_profile(radii[i], 1.0, w_probe);
            normalization += weights[i];
        }

        // Normalize the weights so they sum to unit
        for(size_t i = 0; i < n_focal; i++) 
            weights[i] /= normalization;
    }

}


double elliptical_gaussian_profile(const double I0, const double r, const double theta, const double w_x, const double w_y)
{
    return I0 * exp(-2.0 * r * r * (pow(cos(theta) / w_x, 2) + pow(sin(theta) / w_y, 2)));
}


void get_elliptical_focal_average_weights(const size_t n_r, const size_t n_theta, const double I0, const double w_x, const double w_y, 
                                            const double w_probe, double weights[n_r*n_theta], double intensities[n_r*n_theta])
{
    const double r_max = 1.7 * w_probe;
    const double dr = (double) r_max / n_r;
    const double dtheta = (double) (M_PI / 4.0) / (n_theta - 1.0);
    double normalization = 0.0;
    const double epsilon = 0.1;

    if (n_r * n_theta <= 1)
    {
        intensities[0] = I0;
        weights[0] = 1.0;
    }
    else
    {
        size_t idx = 0;
        for (size_t i = 0; i < n_r; i++)
        {
            double r = epsilon + i * dr; 
            for (int j = 0; j < n_theta; j++)
            {
                double theta = j * dtheta;
                double symmetry_factor = (theta == 0.0 || theta == M_PI / 4.0) ? 2.0 : 4.0; 
                intensities[idx] = elliptical_gaussian_profile(I0, r, theta, w_x, w_y);
                weights[idx] = symmetry_factor * r * gaussian_profile(r, 1.0, w_probe);
                normalization += weights[idx];
                idx++;
            }
        }

        for (size_t i = 0; i < n_r * n_theta; i++)
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
    // If no linewidth average! Otherwise it underestimates B
    if (n == 0)
    {
        B_values[0] = B;
        weights[0] = 1.0;
        return;
    }

    // Find the B-values
    size_t n_B = 2 * n + 1;
    double delta = 3.0 * Delta_B / n_B;
    for (size_t i = 0; i < n_B; i++)
        B_values[i] = B - 3.0 * Delta_B / 2.0 + delta * i;

    // Find the weights
    for (size_t i = 0; i < n_B; i++)
        weights[i] = gaussian_lineshape(B_values[i], B, Delta_B);
    
    // Then normalize so they sum to unity
    double weight_sum = 0.0;
    for (size_t i = 0; i < n_B; i++) weight_sum += weights[i];
    for (size_t i = 0; i < n_B; i++) weights[i] /= weight_sum;
}

double e_field_from_intensity(const double intensity)
{
    return sqrt(2.0 * intensity / (speed_of_light_SI * vaccuum_permitivity_SI)) * au_e_field_strength_per_au;
}


int count_lines(FILE* file)
{
    char buf[BUF_SIZE];
    int counter = 0;
    for (;;)
    {
        size_t res = fread(buf, 1, BUF_SIZE, file);
        if (ferror(file))
            return -1;

        int i;
        for (i = 0; i < res; i++)
            if (buf[i] == '\n')
                counter++;

        if (feof(file))
            break;
    }

    rewind(file);

    return counter;
}