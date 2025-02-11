#ifndef PLANAR_ROTOR_H
#define PLANAR_ROTOR_H

#include "utils.h"

typedef double _Complex dcmplx;

double E_rot_2D(int j, double B);

void get_rot_energies_2D(size_t j_max, double B, double rot_energies[2*j_max + 1]);

void get_exponentials(int j_max, double B, double dt, dcmplx exponentials[2*j_max+1][2*j_max+1]);

void get_field_free_propagator(int j_max, double B, double dt, dcmplx propagator_diag[2*j_max+1]);

extern void multiply_cos2_2D(const size_t j_max, const dcmplx vec[2*j_max+1], dcmplx result[2*j_max+1]);

void field_free_propagation(const size_t j_max, const double dt, const size_t n_steps, size_t *counter, 
                            const double B, dcmplx psi0[2*j_max+1], double cos2_exp[], double weight);

int planar_rotor_ode(double t, const double _psi[], double _psi_deriv[], void* params);

void field_propagation(const size_t j_max, const size_t n_steps, const double dt, const double B, const double fwhm, 
                       const double e_field_sq, const double delta_alpha, const double E_rot[2*j_max+1], 
                       dcmplx psi0[2*j_max+1], double t0, size_t *counter, double cos2_exp[], double weight, const field_params *fp);

void planar_rotor_propagation(solver_params *params, dcmplx psi0[params->dim], double cos2[], double weight);

size_t get_planar_thermal_weights(solver_params *params, double weights[params->dim]);





#endif  /* PLANAR_ROTOR_H */