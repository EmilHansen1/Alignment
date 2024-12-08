#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include "linear_rotor.h"
#include "utils.h"


typedef struct 
{
    size_t j;
    int m;
} free_qn;


/// @brief Function for filling a list with free quantum numbers corresponding to the index
/// @param j_max 
/// @param quantum_numbers A list if free quantum number strutcts
void create_idx_map(const size_t j_max, free_qn quantum_numbers[j_max * j_max])
{
    size_t idx = 0;
    for(size_t j = 0; j < j_max; j++)
    {
        for(int m = -j; m <= j; m++)
        {
            free_qn qn = {j, m};
            quantum_numbers[idx] = qn;
            idx++;
        }
    }
}


/// @brief The energy of a free 3D rotor given the quantum numbers and the rotational constant
/// @param qn The qunatum number struct
/// @param B The rotational constant
/// @return The energy of a free 3D rotor
double E_rot_3D(const free_qn qn, const double B)
{
    return B * qn.j * (qn.j + 1.0);
}


/// @brief Fills an array with the free rotor Hamiltonian eigenvalues (with degeneracy)
/// @param j_max The maximum quantum number (not included)
/// @param B The rotational constant
/// @param energies The list of energies to be filled
void get_rot_energies_3D(const size_t j_max, const double B, double energies[j_max * j_max])
{
    const size_t dim = j_max * j_max;
    
    for(size_t i = 0; i < dim; i++)
    {
        energies[i] = 
    }
    
}