/**
 * \file mass.cpp
 * \brief This target builds a script that computes the mass matrix
 * for the Helmholtz transmission problem using
 * second-kind direct BIEs and Galerkin BEM.
 * The results are written to file.
 * The script can be run as follows:
 * <tt>
 * /path/to/mass \<shape flag\>
 *      \<shape parameters\> \<number of panels\>
 *      \<order of quadrature rule\> \<outputfile\>
 * </tt>
 * SHAPE FLAG may be "circle", "square", or "star".
 *   for "circle", the SHAPE PARAMETERS are the radius;
 *   for "square", the SHAPE PARAMETERS are the half side length;
 *   for "star", the SHAPE PARAMETERS are the inner and outer
 *               radius, separated by a "_";
 *
 * This file is a part of the simpleTBEM library.
 */
#include <iostream>
#include <fstream>
#include "math.h"
#include "build_shapes.cpp"
#include "continuous_space.hpp"
#include "mass_matrix.hpp"

int main(int argc, char** argv) {

    // define number of panels and order of quadrature rule
    unsigned numpanels = atoi(argv[3]);
    unsigned order = atoi(argv[4]);

    // define shape
    std::string shape = argv[1];
    std::string parameters = argv[2];
    std::vector<double> params;
    PanelVector panels;
    if (shape == "circle") {
        params.push_back(atof(parameters.c_str()));
        panels = buildCircle(params[0], numpanels);
    } else if (shape == "square") {
        params.push_back(atof(parameters.c_str()));
        panels = buildSquare(params[0], numpanels);
    } else if (shape == "star") {
        std::istringstream parameterstream(parameters);
        std::string parameter;
        std::getline(parameterstream, parameter, '_');
        params.push_back(atof(parameter.c_str()));
        std::getline(parameterstream, parameter, '_');
        params.push_back(atof(parameter.c_str()));
        panels = buildStar(params[0], params[1], numpanels);
    }

    // set FEM-sapces of lowest order for validation
    ContinuousSpace<1> cont_space;

    // generate mesh on boudary
    ParametrizedMesh mesh(panels);

    // compute mass matrix
    Eigen::MatrixXcd M_cont = mass_matrix::GalerkinMatrix(mesh,cont_space,cont_space,order);

    // write mass matrix to file
    std::ofstream file_matrix;
    file_matrix.open(argv[5], std::ofstream::out);
    file_matrix.precision(10);
    file_matrix << std::scientific;
    for (unsigned i=0; i < numpanels; i++){
        for (unsigned j=0; j < numpanels; j++){
            file_matrix << M_cont(i,j).real();
            if (j < numpanels-1) file_matrix << ",";
        }
        if (i < numpanels-1) file_matrix << std::endl;
    }
    file_matrix.close();

    return 0;
}
