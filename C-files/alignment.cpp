#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include "alignment.h"

solver::solver(solver_params params)
{
    this->params = params;
}

int main()
{
    solver_params params;
    params.B = 1;
    params.I0 = 2;

    return 0;
}