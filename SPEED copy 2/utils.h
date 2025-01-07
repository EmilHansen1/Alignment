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

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

typedef double _Complex dcmplx;


typedef struct
{
    size_t N;
    double *diag;
    double *off_diag; 
} tridiag;

typedef struct 
{
    size_t j;
    int m;
} free_3d_qn;

typedef struct
{
    double B;
    const double delta_alpha;
    double even_abundance;
    double odd_abundance;
    const double delta_B;
    double *E_rot;
    const size_t j_max;
    const double homonuclear_abundance;
    const free_3d_qn *qn; 

} molecule_params;


typedef struct
{
    const size_t n_pulses;
    double fwhm;
    double e_field_squared;
    double intensity_SI;

    const bool custom_pulse_flag;
    const char *custom_pulse_fname;
    gsl_interp_accel *acc;
    gsl_spline *spline;
    
} field_params;


typedef struct 
{
    // Linewidth averaging
    const size_t n_lw;
    const double delta_B;

    // Focal averaging
    const size_t n_theta;
    const size_t n_r;
    const double w_x;
    const double w_y;
    const double w_probe;

    // Thermal averaging
    const double T;
    const double beta;
} avg_params;


typedef enum
{
    PLANAR_ROTOR,
    LINEAR_ROTOR,
    SYMMETRIC_TOP,
    CONSTRAINED_ROTOR
} rotor_type;


typedef struct
{
    molecule_params molecule; 
    field_params field;
    avg_params avg;
    const rotor_type type;
    const double temperature;
    const double beta;
    const size_t dim;
    const size_t n_steps_field;
    const size_t n_steps_field_free;
    const double dt;
    const double t_start;
    const double t_end;
    const char *out_fname;
    const tridiag *diagonals;
} solver_params;

void print_solver_params(const solver_params *params);


typedef struct 
{
    const size_t j_max;         // The maximum angular momentum quantum number
    const double e_field_sq;    // The peak E-field squared
    const double fwhm;          // The full width at half maximum of the of the Gaussian pulse 
    const double *E_rot;        // Array with the rotational energy levels
    const double delta_alpha;   // Polarizability anisotropy
    const field_params *laser;  // Laser params
    const tridiag *diagonals;    // The diagonals of the interaction Hamiltonian
} ode_params; 



void print_matrix(int dim, dcmplx matrix[dim][dim]);

void print_vector(int dim, dcmplx vec[dim]);

dcmplx* matmul_fast_three_band(size_t dim, dcmplx mat[dim][dim], dcmplx vec[dim]);

dcmplx scalar_product(size_t dim, dcmplx vec1[dim], dcmplx vec2[dim]);

double e_field_squared(const double t, const double amplitude_squared, const double fwhm);

double e_field_squared_custom(const double t, const double amplitude_squared, const gsl_spline *spline, gsl_interp_accel *acc);

size_t get_thermal_weights_funciton(solver_params *solver, double weights[solver->dim]);

#endif /* UTILS_H */