/**
 * \file assembler.cpp
 * \brief This target builds a script that assembles the system
 * of the Helmholtz transmission problem using
 * second-kind direct BIEs and Galerkin BEM.
 * The results are written to file.
 * The script can be run as follows:
 * <tt>
 * /path/to/assembler \<shape flag\>
 *      \<shape parameters\> \<refraction inside\>
 *      \<wavenumber\> \<source angle\> \<number of panels\>
 *      \<order of quadrature rule\>
 *      \<whether to rescale neumann part by wavenumber\>
 *      \<matrix1outputfile\> \<matrix2outputfile\>
 *      \<vectoroutputfile\>
 * </tt>
 * SHAPE FLAG may be "circle", "square", "star", or "cshape".
 *   for "circle", the SHAPE PARAMETERS are the radius;
 *   for "square", the SHAPE PARAMETERS are the half side length;
 *   for "star", the SHAPE PARAMETERS are the inner and outer
 *               radius, separated by a "_";
 *   for "cshape", the SHAPE PARAMETERS are the side length and
 *                 inner radius, separated by a "_";
 *
 * This file is a part of the simpleTBEM library.
 */
#include <complex>
#include <iostream>
#include <fstream>
#include "math.h"
#include "build_shapes.cpp"
#include "solvers.hpp"
#include "continuous_space.hpp"

// define shorthand for time benchmarking tools, complex data type and immaginary unit
typedef std::complex<double> complex_t;
complex_t ii = complex_t(0,1.);

int main(int argc, char** argv) {
    
    // define refraction index and wavenumber
    double c_i = atof(argv[3]);
    double c_o = 1.;
    double k = atof(argv[4]);
    double theta = atof(argv[5]);

    // define number of panels and order of quadrature rule
    unsigned numpanels = atoi(argv[6]);
    unsigned order = atoi(argv[7]);

    // define whether to rescale neumann part by wavenumber
    bool rescale_neumann = atoi(argv[8]) != 0;

    // define output filenames
    std::string filename_matrix1 = argv[9];
    std::string filename_matrix2 = argv[10];
    bool also_matrices = (filename_matrix1 != filename_matrix2);

    // define shape
    std::string shape = argv[1];
    std::string parameters = argv[2];
    std::vector<double> params;
    PanelVector panels;
    std::function<double(double, double)> normalAngle;
    if (shape == "circle" || shape == "square") { // single parameter
        params.push_back(atof(parameters.c_str()));
        if (shape == "circle") {
            panels = buildCircle(params[0], numpanels);
            normalAngle = normalCircle;
        } else if (shape == "square") {
            panels = buildSquare(params[0], numpanels);
            normalAngle = normalSquare;
        }
    } else if (shape == "star" || shape == "cshape") { // double parameter
        std::istringstream parameterstream(parameters);
        std::string parameter;
        std::getline(parameterstream, parameter, '_');
        params.push_back(atof(parameter.c_str()));
        std::getline(parameterstream, parameter, '_');
        params.push_back(atof(parameter.c_str()));
        if (shape == "star") {
            panels = buildStar(params[0], params[1], numpanels);
            normalAngle = [&] (double x1, double x2) { return normalStar(params[0], params[1], x1, x2); };
        } else if (shape == "cshape") {
            panels = buildCShape(params[0], params[1], numpanels);
            normalAngle = [&] (double x1, double x2) { return normalCShape(params[0], params[1], x1, x2); };
        }
    }

    // define incoming wave exp(1i * k * ((c, s) . (x1, x2)))
    auto u_i_dir = [&] (double x1, double x2) {
        // simplify parameters
        double x = k * (cos(theta) * x1 + sin(theta) * x2);
        complex_t result = complex_t(cos(x), sin(x));
        return result;
    };
    // define incoming wave normal derivative
    // 1i * k * exp(1i * k * ((c, s) . (x1, x2))) * ((c, s) . (nu1, nu2))
    auto u_i_neu = [&] (double x1, double x2) {
        double x = k * (cos(theta) * x1 + sin(theta) * x2);
        double normal_angle = normalAngle(x1, x2);
        double w_dir = cos(theta) * cos(normal_angle) + sin(theta) * sin(normal_angle);
        complex_t result = ii * k * complex_t(cos(x), sin(x)) * w_dir;
        return result;
    };

    // set FEM-sapces of lowest order for validation
    ContinuousSpace<1> cont_space;

    // generate mesh on boudary
    ParametrizedMesh mesh(panels);

    if (also_matrices) {
        // compute FEM matrices
        std::tuple<Eigen::MatrixXcd, Eigen::MatrixXcd> matrices = tp::direct_second_kind::matrix(mesh, order, k, c_o, c_i, rescale_neumann);

        // generate output files
        std::ofstream file_matrix1;
        std::ofstream file_matrix2;
        file_matrix1.open(argv[9], std::ofstream::out);
        file_matrix1.precision(10);
        file_matrix1 << std::scientific;
        file_matrix2.open(argv[10], std::ofstream::out);
        file_matrix2.precision(10);
        file_matrix2 << std::scientific;

        // write to file
        for (unsigned i=0; i < 2*numpanels; i++){
            for (unsigned j=0; j < 2*numpanels; j++){
                complex_t M1ij = std::get<0>(matrices)(i,j);
                file_matrix1 << M1ij.real() << "," << M1ij.imag();
                complex_t M2ij = std::get<1>(matrices)(i,j);
                file_matrix2 << M2ij.real() << "," << M2ij.imag();
                if (j < 2*numpanels-1) {
                    file_matrix1 << ",";
                    file_matrix2 << ",";
                }
            }
            if (i < 2*numpanels-1) {
                file_matrix1 << std::endl;
                file_matrix2 << std::endl;
            }
        }
        file_matrix1.close();
        file_matrix2.close();
    }

    // compute FEM vector
    Eigen::VectorXcd vector = tp::direct_second_kind::vector(mesh, u_i_dir, u_i_neu, k, rescale_neumann);

    // generate output file
    std::ofstream file_vector;
    file_vector.open(argv[11], std::ofstream::out);
    file_vector.precision(10);
    file_vector << std::scientific;

    // write to file
    for (unsigned i=0; i < 2*numpanels; i++){
        complex_t Vi = vector[i];
        file_vector << Vi.real() << "," << Vi.imag();
        if (i < 2*numpanels-1) file_vector << std::endl;
    }
    file_vector.close();

    return 0;
}
