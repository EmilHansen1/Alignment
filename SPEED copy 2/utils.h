#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <stdbool.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#define FOUR_LOG_TWO 2.772588722239781

typedef double _Complex dcmplx;

typedef struct
{
    double B;
    const double delta_alpha;
    double even_abundance;
    double odd_abundance;
    const double delta_B;
    const double *E_rot;
    const size_t j_max;

} molecule_params;


typedef struct
{
    const size_t n_pulses;
    double fwhm;
    double e_field_squared;

    const bool custom_pulse_flag;
    gsl_interp_accel *acc;
    const gsl_spline *spline;
    
} field_params;


typedef enum
{
    PLANAR_ROTOR,
    LINEAR_ROTOR,
    SYMMETRIC_TOP
} rotor_type;


typedef struct
{
    molecule_params* molecule; 
    field_params *field;
    const rotor_type type;
    const double temperature;
    const double beta;
    const size_t dim;
    const size_t n_steps_field;
    const size_t n_steps_field_free;
    const double dt;
    const double t_start;
    const double t_end;
} solver_params;


void print_matrix(int dim, dcmplx matrix[dim][dim]);

void print_vector(int dim, dcmplx vec[dim]);

dcmplx* matmul_fast_three_band(size_t dim, dcmplx mat[dim][dim], dcmplx vec[dim]);

dcmplx scalar_product(size_t dim, dcmplx vec1[dim], dcmplx vec2[dim]);

double e_field_squared(const double t, const double amplitude_squared, const double fwhm);

double e_field_squared_custom(const double t, const double amplitude_squared, const gsl_spline *spline, gsl_interp_accel *acc);

size_t get_thermal_weights_funciton(solver_params *solver, double weights[solver->dim]);

#endif /* UTILS_H */