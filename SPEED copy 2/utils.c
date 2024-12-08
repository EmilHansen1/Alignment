#include "utils.h"
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

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
    int j_counter = (solver->type == PLANAR_ROTOR) ? -solver->molecule->j_max : 0;
    for(size_t i = 0; i < solver->dim; i++)
    {  
        double abundance = (j_counter % 2) ? solver->molecule->even_abundance : solver->molecule->odd_abundance;
        weights[i] = abundance * exp(-beta * solver->molecule->E_rot[i]);

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
        for(size_t i = solver->molecule->j_max; i < solver->dim; i++)
        {
            size_t idx_left = solver->dim - i - 1;
            size_t idx_right = i;

            if(i == solver->molecule->j_max)
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



