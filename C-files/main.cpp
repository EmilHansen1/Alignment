//
// Created by emil on 3/12/23.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <boost/numeric/odeint.hpp>
#include <boost/phoenix/core.hpp>
#include <time.h>

/// BOOST ODEINT INFO LINK! NOT EASY TO FIND...
/// https://www.boost.org/doc/libs/1_62_0/libs/numeric/odeint/doc/html/index.html
using namespace boost::numeric::odeint;

using dcmplx = std::complex<double>;
using dmat = std::vector<std::vector<double>>;
using cmat = std::vector<std::vector<dcmplx>>;
using dvec = std::vector<double>;
using cvec = std::vector<dcmplx>;

constexpr double FOUR_LOG_TWO = 2.77258872224;
constexpr double TWO_LOG_TWO = 1.38629436112;
constexpr dcmplx imi = dcmplx(0.0, 1.0);


/**
 * The electric dipole approximation Gaussian electric field envelope
 * @param t time in a.u.
 * @param fwhm full width at half maximum
 * @return value of the electric field envelope
 */
inline double electricField(const double t, const double fwhm, const double E0)
{
    return E0 * std::exp(-TWO_LOG_TWO * std::pow(t/fwhm, 2.0)); 
}

inline double electricField2(const double t, const double fwhm, const double E0, const double delay, const double ratio)
{
    return E0 * std::exp(-TWO_LOG_TWO * std::pow(t/fwhm, 2.0))
           + E0 * std::sqrt(ratio) * std::exp(-TWO_LOG_TWO * std::pow((t - delay)/fwhm, 2.0));
}


inline double E_rot(int j, double B)
{
    return B * j*(j + 1.0);
}


double cos2(int lp, int mp, int l, int m)
{
    auto alpha = [&](int j, int m) -> double
    {
        return std::sqrt((j - m + 1.0)*(j + m + 1.0)/((2.0*j + 1.0)*(2.0*j + 3.0)));
    };

    auto beta = [&](int j, int m) -> double
    {
        return std::sqrt((j - m)*(j + m)/((2.0*j - 1.0)*(2.0*j + 1.0)));
    };

    if(mp == m)
    {
        if(std::abs(m) > l || std::abs(m) > lp)
            return 0.0;

        if(lp == l + 2)
            return alpha(l, m)*alpha(l + 1.0, m);
        else if(lp == l)
            return (alpha(l, m)*beta(l + 1.0, m) + alpha(l - 1.0, m)*beta(l, m));
        else if(lp == l - 2)
            return beta(l, m)*beta(l - 1.0, m);
        else
            return 0.0;
    }
    return 0.0;
}

struct TDSE
{
    double offset, fwhm, E0, ratio;
    int jMax;
    dmat& potMel;

    TDSE(double offset, double fwhm, double E0, double ratio, int jMax, dmat& potMel):
        offset(offset),
        fwhm(fwhm),
        E0(E0),
        ratio(ratio),
        jMax(jMax),
        potMel(potMel)
    {
        //  --- EMPTY --- //
    }

    void operator()(const std::vector<dcmplx>& c, std::vector<dcmplx>& dcdt, double t)
    {
        //double eFieldSq = std::pow(electricField(t, fwhm, E0), 2);
        double eFieldSq = std::pow(electricField(t, fwhm, E0), 2);
        for(int i = 0; i < jMax + 1; i++)
        {
            dcdt[i] = -dcmplx(0.0, 1.0) * eFieldSq * ( potMel[i][i]*c[i]
                                                       + (i - 2 < 0 ? 0.0 : potMel[i][i-2]*c[i-2])
                                                       + (i + 2 > jMax ? 0.0 : potMel[i][i+2]*c[i+2]) );
        }
    }
};


struct Observer
{
    dmat& cos2Mel;
    dvec& cos2, times;
    double B;
    int jMax;

    Observer(dmat& cos2Mel, dvec& cos2, dvec& times, double B, int jMax):
        cos2Mel(cos2Mel), cos2(cos2), times(times), B(B), jMax(jMax)
    {
        // --- EMPTY --- //
        //std::cout << "I'm made!" << "\n";
    }

    void operator()(const std::vector<dcmplx>& state, const double t)
    {
        //std::cout << "I'm called!" << "\n";
        times.push_back(t);
        double cos2_exp = 0.0;
        for(int i = 0; i < jMax + 1; i++)
        {
            for(int j = 0; j < jMax + 1; j++)
            {
                dcmplx wiggle = std::exp(imi*dcmplx{E_rot(i, B) - E_rot(j, B), 0.0});
                cos2_exp += std::real(std::conj(state[i])*state[j] * wiggle * dcmplx(cos2Mel[i][j], 0.0));
            }
        }
        cos2.push_back(cos2_exp);
    }
};


int main(int argc, char* argv[])
{
    clock_t start = clock();

    // Conversion factors and constants
    const double c_SI = 299792458.0;
    const double epsilon_0_SI = 8.8541878128e-12;
    const double au_per_s = 2.4188843265857e17;
    const double au_per_ps = 1.0e-12 * au_per_s;
    const double au_per_fs = 1.0e-15 * au_per_s;
    const double cm_m1_per_au = 219474.63068;

    // Parameters for the simulation
    const int j_max = std::atoi(argv[5]);//80;           // Maximum value for orbital angular momentum
    const int m = std::atoi(argv[4]);                // Magnetic quantum number (conserved!!)
    const int J = std::atoi(argv[3]);                // OAM quantum number
    const double I = std::atof(argv[2]) * 1.0e4;    // W/m²  6.0e12
    const double E_0_SI = std::sqrt(2*I/(c_SI*epsilon_0_SI));
    const double E_0 = E_0_SI/(5.14220674763*1.0e11);
    const double alphaParr = std::atof(argv[6]); //98.05;  //791.57057201;
    const double alphaPerp = std::atof(argv[7]); //52.779; // 406.33350678;
    const double B = std::atof(argv[1]); //1.6991682638e-7; //1.9474/cm_m1_per_au;
    const double fwhm = std::atof(argv[8]) * au_per_fs; //200.0 * au_per_fs;
    const double offset = std::atof(argv[9]) * au_per_ps;
    const double intensityRatio = std::atof(argv[10]);
    const double Trev = 2.0 * M_PI / B;

    // Temporal parameters
    double dt = 200.0; // Step size for the ODE
    double t = -1.5 * fwhm;
    int nSteps = (int) (3.0*fwhm + offset)/dt;
    int nSteps2 = 10000; //(int) 1.1*Trev/dt;
    double dt2 = (double) (1.1*Trev - nSteps*dt)/nSteps2; // Step size for the field-free propagation

    std::cout << "Number of steps: " << nSteps << " " << nSteps2 << "\n";

    // Calculate interaction matrix elements for given m and all possible js.
    dmat Vmat(j_max + 1);
    dmat cos2Mat(j_max + 1);
    std::ofstream cos2File{"cos2.out"};
    std::ofstream VFile{"matElem.out"};
    for(int i = 0; i < j_max + 1; i++)
    {
        for(int j = 0; j < j_max + 1; j++)
        {
            double V = -(alphaParr - alphaPerp)/4.0 * cos2(i, m, j, m) + alphaPerp * (i == j ? 1.0 : 0.0);
            Vmat[i].push_back(V);
            cos2Mat[i].push_back(cos2(i, m, j, m));

            cos2File << cos2Mat[i][j] << " ";
            VFile << Vmat[i][j] << " ";
        }
        cos2File << "\n";
        VFile << "\n";
    }

    // Calculate the wiggle factors for ODE and field-free prop.
    cmat wiggleODE(j_max + 1);
    cmat wiggleProp(j_max + 1);
    for(int i = 0; i < j_max + 1; i++)
    {
        for(int j = 0; j < j_max + 1; j++)
        {
            wiggleODE[i].push_back(std::exp(imi * dcmplx(E_rot(i, B) - E_rot(j, B)) * dt));
            wiggleProp[i].push_back(std::exp(imi * dcmplx(E_rot(i, B) - E_rot(j, B)) * dt2));
        }
    }

    // Initialize initial state
    std::vector<dcmplx> state(j_max + 1, dcmplx{0.0, 0.0});
    state[J] = 1.0;

    // Solve the differential equation
    dvec times, cos2_exp;
    std::vector<cvec> states;

    //Observer obs(cos2Mat, cos2_exp, times, B, j_max);   // TODO: ..fix this? Maybe not..
    TDSE tdse(offset, fwhm, E_0, intensityRatio, j_max, Vmat);
    runge_kutta_fehlberg78<cvec> stepper;

    try
    {
        //steps = integrate_adaptive(make_controlled<stepper_78>(1.0e-12, 1.0e-12), tdse, initialState, 0.0, 15*au_per_ps, 1.0, obs);
        //steps = integrate(tdse, initialState, 0.0, 10.0*au_per_ps, 10.0, obs);
        for(size_t i = 0; i < nSteps; i++)
        {
            times.push_back(t);
            double cos2_val = 0.0;
            for(int i = 0; i < j_max + 1; i++)
            {
                for(int j = 0; j < j_max + 1; j++)
                {
                    dcmplx wiggle = std::exp(imi*dcmplx{(E_rot(i, B) - E_rot(j, B))*t, 0.0});
                    cos2_val += std::real(std::conj(state[i])*state[j] * wiggle * dcmplx(cos2Mat[i][j], 0.0));
                }
            }
            cos2_exp.push_back(cos2_val);
            states.push_back(state);
            stepper.do_step(tdse, state, t, dt);
            t += dt;
        }

        std::cout << "ODE done" << "\n";
        // Field-free propagation after the pulse: state is constant! Only wiggle factors
        for(int i = 0; i < nSteps2; i++)  //(int i = 0; i < nSteps2 - nSteps; i++)
        {
            double cos2_val = 0.0;
            times.push_back(t);
            for(int j = 0; j < j_max + 1; j++)
            {
                dcmplx wigglePlus2 = std::exp(imi*dcmplx{B*(4.0*j + 6.0)*t, 0.0});
                dcmplx wiggleMinus2 = std::exp(imi*dcmplx{B*(2.0 - 4.0*j)*t, 0.0});
                dcmplx term0 = std::pow(std::abs(state[j]), 2) * dcmplx(cos2Mat[j][j], 0.0);
                dcmplx termPlus2 = (j + 2 > j_max ? 0.0 : std::conj(state[j+2])*state[j] * dcmplx(cos2Mat[j+2][j], 0.0) * wigglePlus2);
                dcmplx termMinus2 = (j - 2 < 0 ? 0.0 : std::conj(state[j-2])*state[j] * dcmplx(cos2Mat[j-2][j], 0.0) * wiggleMinus2);
                cos2_val += std::real(termMinus2 + term0 + termPlus2);
            }
            cos2_exp.push_back(cos2_val);
            states.push_back(state);
            t += dt2;
        }
        std::cout << "Propagation done" << "\n";
    }
    catch(std::exception& e)
    {
        std::cerr << e.what();
    }

    clock_t stop = clock();
    double elapsed = (double) (stop - start) / CLOCKS_PER_SEC;
    printf("\nTime elapsed: %.5f\n", elapsed);

    // Print results to file
    std::ofstream outFile{"dat.out"};
    std::ofstream stateFile{"states.out"};
    for(int i = 0; i < times.size(); i++)
    {
        outFile << times[i] << " " << cos2_exp[i] << "\n";
        for(int j = 0; j < j_max + 1; j++)
        {
            stateFile << (double) std::pow(std::abs(states[i][j]), 2) << " ";
        }
        stateFile << "\n";
    }
    stateFile.close();
    outFile.close();

    return 0;   // 'Grrreat succes' - Borat
}
