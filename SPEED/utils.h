#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <stdbool.h>

typedef double _Complex dcmplx;

typedef struct
{
    const double B;
    const double delta_alpha;
    const double even_abundance;
    const double odd_abundance;
    const double delta_B;
    const double *E_rot;
    const size_t j_max;

} molecule_params;


typedef struct
{
    const size_t n_pulses;
    const double *fwhm;
    const double *e_field_squared;
    const bool custom_pulse_flag;
    
} laser_params;


void print_matrix(int dim, dcmplx matrix[dim][dim]);

void print_vector(int dim, dcmplx vec[dim]);

dcmplx* matmul_fast_three_band(size_t dim, dcmplx mat[dim][dim], dcmplx vec[dim]);

dcmplx scalar_product(size_t dim, dcmplx vec1[dim], dcmplx vec2[dim]);

double e_field_squared(const double t, const double amplitude_squared, const double fwhm);

#endif /* UTILS_H */