/**
 * \file solvers.hpp
 * \brief This file defines lowest order solvers for direct first kind BIEs
 * for the Dirichlet and the Neumann problem as well as a direct second kind
 * BIE for the Helmholtz transmission problem.
 *
 * This file is a part of the simpleTBEM library.
 * It was adapted from the HelmholtzTransmissionProblem library.
 */

#ifndef SOLVERSHPP
#define SOLVERSHPP

#include "parametrized_mesh.hpp"
/**
 * \namespace bvp
 * \brief This namespace contains solvers for boundary value problems.
 */
namespace bvp {

/**
 * \namespace direct_first_kind
 * \brief This namespace contains solvers using direct first kind BIE
 * for the Dirichlet and the Neumann problem.
 * The solvers use the lowest order BEM spaces for computation.
 */
    namespace direct_first_kind {

        /**
         * This function solves the Dirichlet problem given a mesh \p mesh, the Dirichlet data of u
         * \p u_dir, the order of the quadrature ruile used to compute the Galerkin matrix entries
         * \p order and the wavenumber \p k.
         * @param mesh mesh of the boundary on which to compute BIOs
         * @param u_dir Dirichlet data
         * @param order order of qudrature rule for matrix entries
         * @param k wavenumber
         * @return Neumann data of u
         */
        Eigen::VectorXcd solve_dirichlet(const ParametrizedMesh &mesh,
                                         const std::function<std::complex<double>(double, double)> u_dir,
                                         const unsigned order,
                                         const double k);

        /**
         * This function solves the Neumann problem given a mesh \p mesh, the Neumann data of u
         * \p u_dir, the order of the quadrature ruile used to compute the Galerkin matrix entries
         * \p order and the wavenumber \p k.
         * @param mesh mesh of the boundary on which to compute BIOs
         * @param u_neu Dirichlet data
         * @param order order of qudrature rule for matrix entries
         * @param k wavenumber
         * @return Dirichlet data of u
         */
        Eigen::VectorXcd solve_neumann(const ParametrizedMesh &mesh,
                                       const std::function<std::complex<double>(double, double)> u_neu,
                                       const unsigned order,
                                       const double k);
    } // namespace direct_first_kind
} // namespace bvp

/**
 * \namespace tp
 * \brief This namespace contains the solver for the Helmholtz transmission problem.
 */
namespace tp {

    /**
	 * \namespace direct_second_kind
	 * \brief This namespace contains the solver using direct second kind BIEs
     * for the Helmholtz transmission problem.
     * The solver uses the lowest order BEM spaces for computation.
     */
    namespace direct_second_kind {
        /**
         * This function returns the matrices of the Helmholtz transmission problem
         * on boundary given by \p mesh. The wavenumber is set by \p k and th refraction
         * indices by \p c_o and \p c_i. The Galerkin matrix entries are computed
         * with a quadrature rule defined by the parameter \p order.
         * @param mesh mesh of the boundary on which to compute BIOs
         * @param order order of qudrature rule for matrix entries
         * @param k wavenumber
         * @param c_o refraction index outer domain
         * @param c_i refraction index on inner domain
         * @param rescale_neumann whether to rescale Neumann portions by wavenumber
         * @return tuple of matrices
         */
        std::tuple<Eigen::MatrixXcd, Eigen::MatrixXcd> matrix(const ParametrizedMesh &mesh,
                                                              const unsigned order,
                                                              const double k,
                                                              const double c_o,
                                                              const double c_i,
                                                              bool rescale_neumann = false);
        /**
         * This function returns the vector of the Helmholtz transmission problem
         * on boundary given by \p mesh for an incoming wave defined by \p u_inc_dir and
         * \p u_inc_neu. The wavenumber is set by \p k and th refraction indeces by
         * \p c_o and \p c_i. The Galerkin matrix entries are compute with a quadrature rule
         * defined by the parameter \p order.
         * @param mesh mesh of the boundary on which to compute BIOs
         * @param u_inc_dir Dirichlet data of incoming wave
         * @param u_inc_neu Neumann data of incoming wave
         * @param k wavenumber
         * @param rescale_neumann whether to rescale Neumann portions by wavenumber
         * @return Dirichlet and Neumann data of incoming wave
         */
        Eigen::VectorXcd vector(const ParametrizedMesh &mesh,
                                const std::function<std::complex<double>(double, double)> u_inc_dir,
                                const std::function<std::complex<double>(double, double)> u_inc_neu,
                                const double k,
                                bool rescale_neumann = false);
        /**
         * This function returns the solution to the Helmholtz transmission problem
         * on boundary given by \p mesh for an incoming wave defined by \p u_inc_dir and
         * \p u_inc_neu. The wavenumber is set by \p k and th refraction indeces by
         * \p c_o and \p c_i. The Galerkin matrix entries are compute with a quadrature rule
         * defined by the parameter \p order.
         * @param mesh mesh of the boundary on which to compute BIOs
         * @param u_inc_dir Dirichlet data of incoming wave
         * @param u_inc_neu Neumann data of incoming wave
         * @param order order of qudrature rule for matrix entries
         * @param k wavenumber
         * @param c_o refraction index outer domain
         * @param c_i refraction index on inner domain
         * @param rescale_neumann whether to rescale Neumann portions by wavenumber
         * @return Dirichlet and Neumann data of resulting wave
         */
        Eigen::VectorXcd solve(const ParametrizedMesh &mesh,
                               const std::function<std::complex<double>(double, double)> u_inc_dir,
                               const std::function<std::complex<double>(double, double)> u_inc_neu,
                               const unsigned order,
                               const double k,
                               const double c_o,
                               const double c_i,
                               bool rescale_neumann = false);
    } // namespace direct_second_kind
} // namespace tp
#endif // DIRICHLETHPP
