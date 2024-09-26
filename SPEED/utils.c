#include "utils.h"

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
dcmplx scalar_product(size_t dim, dcmplx vec1[dim], dcmplx vec2[dim])
{
    dcmplx result = 0.0 + 0.0*I;
    for(size_t i = 0; i < dim; i++)
        result += conj(vec1[i]) * vec2[i];
    return result;
}


double inline e_field_squared(const double t, const double amplitude_squared, const double fwhm)
{   
    return amplitude_squared * exp(-2.772588722239781 * t*t / (fwhm*fwhm));  // The number is 4*log(2)
}

