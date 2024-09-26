#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <complex>
#include <functional>

// Type abbreviations
using dcmplx = std::complex<double>;
using dvec = std::vector<double>;
using cvec = std::vector<dcmplx>;

typedef struct 
{
    const size_t n_pulses;
    const std::function<double(double)>& e_field_squared;
    
} field_params;


enum rotor_type 
{
    PLANAR_ROTOR,
    LINEAR_ROTOR,
    SYMMETIC_TOP
};


struct molecule_params
{
    const double rotational_constant;
    const dvec energy_levels;
    const double polarizability_anisotropy;
    const std::pair<double, double> even_odd_abundance;
    const size_t j_max;
};


struct tdse_params
{
    field_params field_paramters;
    rotor_type type;
    molecule_params molecule_parameters;
    T interaction_matrix;
};


#endif // UTILS_H