/*
 * computeNextStep.h
 *
 ** Lagrangian Reachtubes: The next Generation
 *
 *      Authors: Sophie Gruenbacher
 *      Contact: sophie.gruenbacher@tuwien.ac.at
 */

#ifndef COMPUTE_NORM_H
#define COMPUTE_NORM_H

/* Standard C++ library*/
#include <iostream>

/* Eigen Library*/
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

/* IBEX */
#include "ibex.h"

/* Libraries defined by myself */
#include "ibex_matrixInv.h"

/*Define namespace*/
using namespace Eigen;
using namespace std;
using namespace ibex;


/*function computes the new best metric and returns usedC and usedCi as references*/
// get Metric parallel to axis of ellipse, returning also length of semi-axes
bool getMetric(IntervalMatrix Fmid, IntervalMatrix oldC, IntervalMatrix oldCi,
		IntervalMatrix& usedC, IntervalMatrix& usedCi, int dim, ibex::Vector& semiAxis, int variablesDim)
{

	usedCi.put(0,0,Fmid * oldCi.submatrix(0,variablesDim-1,0,variablesDim-1)); //take submatrix if the model has time as the last variable

	usedC = gaussInverseMatrix(usedCi);

	// The following part is important to get xbox with the axes parallel to the major and minor axes of the ellipse
	IntervalMatrix Mmat = usedC.transpose() * usedC;

    /* Compute Eigenvalue and Eigenvector in Eigen */
    MatrixXd A = ibex_to_eigen(Mmat.mid(), dim);
    //VectorXd rowj_eigen = VectorXd::Zero(dim);
    EigenSolver<MatrixXd> es(A);

    MatrixXd rVec = es.eigenvectors().real();
    VectorXd rowj_eigen = es.eigenvalues().real();


	for (int j=0;j<dim;j++){
		semiAxis[j] = 1./(sqrt(rowj_eigen(j)));
		VectorXd rVecCol_j = rVec.col(j);
		for (int i = 0; i<dim; i++){
			usedCi[i][j] = semiAxis[j] * rVecCol_j(i);
		}
	}

	usedC = gaussInverseMatrix(usedCi);


	bool normChanged = true;

	return normChanged;
}

#endif /* COMPUTE_NORM_H */
