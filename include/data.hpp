/**
 * \file find_roots.hpp
 * \brief This file defines different root finding methods applied
 * to find the minimas in the smallest singular value
 * of the Helmholtz transmission problem solutions operator.
 *
 * This file is a part of the simpleTBEM library.
 * It was adapted from the HelmholtzTransmissionProblem library.
 */
#ifndef DATA_HPP
#define DATA_HPP

#include <vector>
#include <iomanip>
#include <ostream>

using namespace std;

enum flag {
    active,
    nozero,
    zerofound
};

struct data {
public:
    double grid_point;
    double value;
    double derivative;
    flag flag_val;

};
#endif //OPERATORS
