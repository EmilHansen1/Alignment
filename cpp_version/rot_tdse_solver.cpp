#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <boost/numeric/odeint.hpp>

using dcmplx = std::complex<double>;
using dvec = std::vector<double>;
using cvec = std::vector<dcmplx>;
using dmat = std::vector<std::vector<double>>;
using cmat = std::vector<std::vector<dcmplx>>;

const double H_BAR = 1.054571817e-34;      // The reduced Planck's constant
const double SPEED_OF_LIGHT = 299792458.0; // The speed of light in SI-units
const double TWO_LOG_TWO = 1.38629436112;  // Constant for electric field
const double FOUR_LOG_TWO = 2.77258872224; // Constant for electric field squared
constexpr dcmplx imi = dcmplx(0.0, 1.0);   // The imaginary unit

/**
 * @brief The rotational energy
 *
 * @param j the angular momentum quantum number
 * @param B the rotational constant in rad/s
 * @return double
 */
inline double E_rot(const size_t j, const double B)
{
    return H_BAR * B * j * (j + 1.0);
}

inline double omega_rot(const size_t j, const double B)
{
    return B * j * (j + 1.0);
}

inline double E_field_envelope(const double t, const double amplitude, const double fwhm)
{
    return amplitude * std::exp(-TWO_LOG_TWO * std::pow(t / fwhm, 2.0));
}

inline double E_field_envelope_squared(const double t, const double amplitude, const double fwhm)
{
    return amplitude * std::exp(-FOUR_LOG_TWO * std::pow(t / fwhm, 2.0));
}

inline double cos2(const size_t lp, const int mp, const size_t l, const int m)
{
    auto alpha = [&](int j, int m) -> double
    {
        return std::sqrt((j - m + 1.0) * (j + m + 1.0) / ((2.0 * j + 1.0) * (2.0 * j + 3.0)));
    };

    auto beta = [&](int j, int m) -> double
    {
        return std::sqrt((j - m) * (j + m) / ((2.0 * j - 1.0) * (2.0 * j + 1.0)));
    };

    if (mp == m)
    {
        if (std::abs(m) > l || std::abs(m) > lp)
            return 0.0;

        if (lp == l + 2)
            return alpha(l, m) * alpha(l + 1.0, m);
        else if (lp == l)
            return (alpha(l, m) * beta(l + 1.0, m) + alpha(l - 1.0, m) * beta(l, m));
        else if (lp == l - 2)
            return beta(l, m) * beta(l - 1.0, m);
        else
            return 0.0;
    }
    return 0.0;
}

dmat getPotMatrixElements(const size_t jMax, const int m0, const double alphaPerp, const double alphaParr)
{
    dmat mat(jMax + 1);
    for (int i = 0; i < jMax + 1; i++)
    {
        for (int j = 0; j < jMax + 1; j++)
        {
            double V = -(alphaParr - alphaPerp) / 4.0 * cos2(i, m0, j, m0);// + alphaPerp * (i == j ? 1.0 : 0.0);
            mat[i].push_back(V);
        }
    }
    return mat;
}

struct SolverParams
{
    const size_t jMax;
    const size_t j0;
    const int m0;
    const double B;
    const double E0;
    const dmat& potMel;
    dvec& times;
    std::vector<cvec>& states;

    SolverParams(const size_t jMax, const size_t j0, const int m0, const double B, const double E0, const dmat& potMel, 
                 dvec& times, std::vector<cvec>& states):
        jMax(jMax), j0(j0), m0(m0), B(B), E0(E0), potMel(potMel), times(times), states(states) { /* Empty ctor */ }
};



void interactionProgagation(SolverParams &params)
{

}

void fieldFreePropagation(SolverParams &params)
{
    // To be implementd
}



int main()
{
    // Set up constants
    const size_t j0 = 0;
    const size_t jMax = 50;
    const int m0 = 0;
    const double B = 1.1e9;     // MHz
    const double E0 = 1.1;
    const double alphaPerp = 2.2;
    const double alphaParr = 3.3;

    // State and observble related objects
    cvec state(jMax + 1);
    state[j0] = 1.0;
    std::vector<cvec> states;
    dvec times;
    dvec cos2ExpVal;
    dvec JSquaredExpVal;

    // Calculate matrix elements and create params objet
    dmat potMel = getPotMatrixElements(jMax, m0, alphaPerp, alphaParr);
    SolverParams params(jMax, j0, m0, B, E0, potMel, times, states);

    // Propagation under influence of laser field
    try
    {
        interactionProgagation(params);
    }
    catch(std::exception& e)
    {
        std::cout << e.what() << "\n";
        return -1;
    }

    // Field-free propagation
    fieldFreePropagation(params);

    // Calculate expectation calues

    std::cout << "Done John" << "\n";
    return 0;
}