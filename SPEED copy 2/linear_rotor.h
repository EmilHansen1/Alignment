#ifndef LINEAR_ROTOR_H
#define LINEAR_ROTOR_H

#include <stdio.h>
#include "utils.h"

typedef struct
{
    double *diag;
    double *offdiag;
} cos2_diags;

void create_idx_map(const size_t j_max, free_3d_qn quantum_numbers[j_max * j_max]);

void fill_diagonals(const size_t j_max, const free_3d_qn quantum_numbers[j_max * j_max], tridiag *diags);

double E_rot_3D(const free_3d_qn qn, const double B);

void get_rot_energies_3D(const size_t j_max, const double B, double energies[j_max * j_max], const free_3d_qn quantum_numbers[j_max * j_max]);

void get_field_free_prpagator_3D(const size_t j_max, const double B, double dt, const free_3d_qn quantum_numbers[j_max * j_max], dcmplx propagator_diag[j_max * j_max]);

void linear_rotor_propagation(solver_params *params, dcmplx psi0[params->dim], double cos2[], double weight);

#endif // LINEAR_ROTOR_H