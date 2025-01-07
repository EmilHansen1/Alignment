#include "utils.h"
#include "units.h"
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

/**
 * @brief prints the solver parameters
 * 
 * @param params the solver parameters
 */
void print_solver_params(const solver_params *params)
{
    printf("dim = %zu\n", params->dim);
    printf("temperature = %.6f K\n", params->temperature);
    printf("t_start = %.6f, t_end = %.6f, dt = %.6f\n", params->t_start, params->t_end, params->dt);
    printf("molecule B = %.6f GHz, delta_alpha = %.6f Å, j_max = %zu\n",
           params->molecule.B / au_frequency_per_SI / 2.0 / M_PI / 1.0e9, params->molecule.delta_alpha / pow(au_per_angstrom, 3.0), params->molecule.j_max);
    printf("field E_field^2 = %.6f, fwhm = %.6f ps\n", params->field.e_field_squared, params->field.fwhm / au_per_ps);
    printf("custom_pulse_flag = %d\n", params->field.custom_pulse_flag);
    printf("avg w_x = %.6f, w_y = %.6f, n_lw = %zu\n", params->avg.w_x, params->avg.w_y, params->avg.n_lw);
}

/**
 * @brief prints a standard 2D array
 * 
 * @param dim 
 * @param matrix 
 */
void print_matrix(int dim, dcmplx matrix[dim][dim])
{
    for(size_t i = 0; i < dim; i++)
    {
        for(size_t j = 0; j < dim; j++)
            printf("%f + I%f\t", creal(matrix[i][j]), cimag(matrix[i][j]));
        printf("\n");
    }
}


/**
 * @brief prints a dcmplx vector
 * 
 * @param dim the dimension of the hilbert space
 * @param vec the dcmplx vector to be printed
*/
void print_vector(int dim, dcmplx vec[dim])
{
    for(size_t i = 0; i < dim; i++)
        printf("%f + I%f\t", creal(vec[i]), cimag(vec[i]));
    printf("\n");
}


/**
 * @brief performs fast 3-band matrix multiplication
 * 
 * @param dim 
 * @param mat 
 * @param vec 
 * @return dcmplx* 
 */
dcmplx* matmul_fast_three_band(size_t dim, dcmplx mat[dim][dim], dcmplx vec[dim])
{
    
    return vec;
}


/**
 * @brief performs the scalar product of two complex-valued arrays
 * 
 * @param dim 
 * @param vec1 
 * @param vec2 
 * @return dcmplx 
 */
dcmplx inline scalar_product(size_t dim, dcmplx vec1[dim], dcmplx vec2[dim])
{
    dcmplx result = 0.0 + 0.0*I;
    for(size_t i = 0; i < dim; i++)
        result += conj(vec1[i]) * vec2[i];
    return result;
}


/**
 * @brief The electric field squared involved in the interaction Hamiltonia
 * 
 * @param t the time in a.u.
 * @param amplitude_squared the square of the electric field amplitude in a.u.
 * @param fwhm the temporal full width at half maximum in a.u.
 * @return double 
 */
double inline e_field_squared(const double t, const double amplitude_squared, const double fwhm)
{   
    return amplitude_squared * exp(-FOUR_LOG_TWO * t*t / (fwhm*fwhm));
}

/**
 * @brief The electric field squared for a custom pulse
 * 
 * @param t the time in a.u.
 * @param amplitude_squared the amplitude of the e-field in a.u.
 * @param spline the spline object used
 * @param acc the spline accelerator
 */
double e_field_squared_custom(const double t, const double amplitude_squared, const gsl_spline *spline, gsl_interp_accel *acc)
{
    return amplitude_squared * gsl_spline_eval(spline, t, acc);
}


/**
 * @brief Get the thermal weights
 * 
 * @param dim the dimension of the Hilbert space 
 * @param hamiltonian the field-free diagonal of the Hamiltonian
 * @param temperature the ensemble temperautre in K
 * @param weights list to be overwritten with the thermal weights of each state
 * @return size_t the the maximum j quantum number to be in the ensemble
 */
size_t get_thermal_weights_funciton(solver_params *solver, double weights[solver->dim])
{
    // The cursed constant convertes the temperature in K to energy in units of kB. I'll unfuck it later...
    double beta = 1.0/(1.380649e-23 / 4.3597447222071e-18 * solver->temperature);
    size_t ensemble_size = 0;
    double partition_funciton = 0.0;

    // Calculate the individual Boltzmann factors and weights based on rotor type
    int j_counter = (solver->type == PLANAR_ROTOR) ? -solver->molecule.j_max : 0;
    for(size_t i = 0; i < solver->dim; i++)
    {  
        double abundance = (j_counter % 2) ? solver->molecule.even_abundance : solver->molecule.odd_abundance;
        weights[i] = abundance * exp(-beta * solver->molecule.E_rot[i]);

        // Calculate partition function directly for non-planar rotors
        if(solver->type != PLANAR_ROTOR)
        {
            partition_funciton += weights[i];
            ensemble_size++;
            if(partition_funciton > 0.99)
                break;
        }
        j_counter++;
    }

    // 'Trim off' the unneccesary parts of the ensemble and calculate the partition function
    if(solver->type == PLANAR_ROTOR)
    {   
        // Start from the middle of the weights (ground state) and work left and right simultaneously
        for(size_t i = solver->molecule.j_max; i < solver->dim; i++)
        {
            size_t idx_left = solver->dim - i - 1;
            size_t idx_right = i;

            if(i == solver->molecule.j_max)
            {
                partition_funciton += weights[i];
            }
            else
            {
                partition_funciton += weights[idx_left] + weights[idx_right];
            }
            ensemble_size++;

            // Set weights to zero if we have 99% of the ensemble
            if(partition_funciton > 0.99)
            {
                weights[idx_left] = 0.0;
                weights[idx_right] = 0.0;
                break;
            }
        }
    }

    // Divide all factors by the partition function to normalize
    for(size_t i = 0; i < solver->dim; i++)
        weights[i] /= partition_funciton;

    return ensemble_size;
}

