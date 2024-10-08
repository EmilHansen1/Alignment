#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <stdbool.h>

#define FOUR_LOG_TWO 2.772588722239781

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
} field_params;

typedef enum
{
    PLANAR_ROTOR,
    LINEAR_ROTOR,
    SYMMETRIC_TOP
} rotor_type;

typedef struct
{
    const molecule_params* molecule; 
    const field_params *field;
    const rotor_type type;
    const double temperature;
    const size_t dim;
    const size_t n_steps_field;
    const size_t n_steps_field_free;
} solver_params;


void print_matrix(int dim, dcmplx matrix[dim][dim]);

void print_vector(int dim, dcmplx vec[dim]);

dcmplx* matmul_fast_three_band(size_t dim, dcmplx mat[dim][dim], dcmplx vec[dim]);

dcmplx scalar_product(size_t dim, dcmplx vec1[dim], dcmplx vec2[dim]);

double extern e_field_squared(const double t, const double amplitude_squared, const double fwhm);

size_t get_thermal_weights_funciton(solver_params *solver, double weights[solver->dim]);

#endif /* UTILS_H */