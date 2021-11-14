/*
 * ibex_to_eigen.h
 *
 ** Lagrangian Reachtubes: The next Generation
 *
 *      Authors: Sophie Gruenbacher and Md Ariful Islam
 *      Contact: sophie.gruenbacher@tuwien.ac.at
 */

#ifndef IBEX_TO_EIGEN_H
#define IBEX_TO_EIGEN_H

/* Eigen Library*/
#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>
#include <Eigen/Eigenvalues>
#include <Eigen/Dense>

/* ibex Library*/
#include "ibex.h"

/*Define namespace*/
using namespace Eigen;
using namespace ibex;

////////////// BEGIN: API between IBEX and EIGEN ////////////////////////////////
/* IBEX to Eigen */
MatrixXd ibex_to_eigen(ibex::Matrix ibex_mat, int dim) {
    MatrixXd eig_mat(dim, dim);
    for(int i=0; i<dim; i++){
        for(int j=0; j<dim; j++){
            eig_mat(i,j) = ibex_mat[i][j];
        }
    }
    return eig_mat;
}

VectorXd ibex_to_eigen(ibex::Vector ibex_vec, int dim){
    VectorXd eig_vec(dim);
    for(int i=0; i<dim; i++){
        eig_vec(i) = ibex_vec[i];
    }
    return eig_vec;
}

/* Eigen to IBEX */
// MatrixXd -> Matrix
ibex::Matrix eigen_to_ibex(MatrixXd eig_mat, int dim) {
    ibex::Matrix ibex_mat(dim, dim);
    for(int i=0; i<dim; i++){
        for(int j=0; j<dim; j++){
            ibex_mat[i][j] = eig_mat(i,j);
        }
    }
    return ibex_mat;
}

// VectorXd -> IntervalVector
ibex::IntervalVector eigen_to_ibex(VectorXd eig_vec, int dim){
    ibex::IntervalVector ibex_vec(dim);
    for(int i=0; i<dim; i++){
        ibex_vec[i] = eig_vec(i);
    }
    return ibex_vec;
}

#endif /* IBEX_TO_EIGEN_H */
