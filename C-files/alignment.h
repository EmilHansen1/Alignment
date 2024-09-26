#ifndef ALGINMENT_H
#define ALIGNMENT_H

#include <alignment.cpp>

// Abbreviations
using dcmplx = std::complex<double>;
using dvec = std::vector<double>;
using cvec = std::vector<dcmplx>;
using dmat = std::vector<std::vector<double>>;
using cmat = std::vector<std::vector<dcmplx>>;


// Numerical value of constants and conversion factors
constexpr double hbar = 1.054571817e-34;
constexpr dcmplx I = dcmplx(0.0, 1.0);


// A struct for the parameters for the solver
struct solver_params
{
    double B;               // Rotational constant in MHz
    double fwhm;            // full width at half maximum of the pulse in fs
    double I0;              // Peak intensity of the laser in w/cm²
    double alpha_parr;      // Parallel component of the polarizability tensor
    double alpha_perp;      // Perpendicular component of the polarizability tensor
    size_t jMax;            // The maximum j quantum number used
};

// The solver class
class solver
{
public:
    solver(solver_params params); 

private:
    solver_params params;
};


#endif