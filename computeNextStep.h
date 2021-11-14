/*
 * computeNextStep.h
 *
 ** Lagrangian Reachtubes: The next Generation
 *
 *      Authors: Sophie Gruenbacher
 *      Contact: sophie.gruenbacher@tuwien.ac.at
 */

#ifndef COMPUTE_NEXT_STEP_H
#define COMPUTE_NEXT_STEP_H

#include <iostream>
#include <fstream>
#include <string>
#include "stdlib.h"
#include <math.h>

/* EIGEN */
#include <Eigen/SVD>
#include <Eigen/QR>

/* IBEX */
#include "ibex.h"

/* Libraries defined by myself */
#include "computeIntEig.h"
#include "computeNorm.h"
#include "ibex_to_eigen_api.h"
#include "ibex_matrixInv.h"

using namespace std;
using namespace ibex;

//Takes the intersection between the interval arithmetics evaluation and the mean value form evaluation of Function f on the interval xbox
IntervalVector evalMeanValue(Function &f, IntervalMatrix &fjacob_xbox, IntervalVector &xbox)
{
	IntervalVector returnVec = f.eval_vector(xbox.mid()) + fjacob_xbox * (xbox - xbox.mid());
	returnVec &= f.eval_vector(xbox);
	return returnVec;
}

/* d^2/dt x(t) = d/dt f(x(t)) */
IntervalVector evalSecondTaylor(Function &f, Function &fjacob, IntervalVector roughEnc)
{
	IntervalMatrix fjacob_roughEnc = fjacob.eval_matrix(roughEnc);
	return fjacob_roughEnc * evalMeanValue(f, fjacob_roughEnc, roughEnc);
}

/* Functions for tensor-vector multiplications */
IntervalMatrix mul2TensorVector(std::vector<ibex::Function> &d2f, IntervalVector x, IntervalVector vec, int dim)
{
	IntervalMatrix d2f_vec = ibex::Matrix::zeros(dim);

	for (int i = 0; i < dim; i++)
	{
		d2f_vec[i] = d2f[i].eval_matrix(x) * vec; // Jacobian of f_i evaluated at x times vector vec
	}

	return d2f_vec;
}

IntervalMatrix mul3TensorVector(std::vector<std::vector<ibex::Function>> &d3f, IntervalVector x,
								IntervalVector vec1, IntervalVector vec2, int dim)
{
	IntervalMatrix d3f_vec1_i = ibex::Matrix::zeros(dim);
	IntervalMatrix d3f_vec1_vec2 = ibex::Matrix::zeros(dim);

	for (int i = 0; i < dim; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			d3f_vec1_i[j] = d3f[i][j].eval_matrix(x) * vec1; //Hessian of f_i_j evaluated at x times vector vec1
		}
		d3f_vec1_vec2[i] = d3f_vec1_i * vec2;
	}

	return d3f_vec1_vec2;
}

IntervalMatrix mul4TensorVector(std::vector<std::vector<std::vector<ibex::Function>>> &d4f, IntervalVector x,
								IntervalVector vec1, IntervalVector vec2, IntervalVector vec3, int dim)
{
	IntervalMatrix d4f_vec1_i_j = ibex::Matrix::zeros(dim);
	IntervalMatrix d4f_vec1_vec2_i = ibex::Matrix::zeros(dim);
	IntervalMatrix d4f_vec1_vec2_vec3 = ibex::Matrix::zeros(dim);

	for (int i = 0; i < dim; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			for (int k = 0; k < dim; k++)
			{
				d4f_vec1_i_j[k] = d4f[i][j][k].eval_matrix(x) * vec1;
			}
			d4f_vec1_vec2_i[j] = d4f_vec1_i_j * vec2;
		}
		d4f_vec1_vec2_vec3[i] = d4f_vec1_vec2_i * vec3;
	}

	return d4f_vec1_vec2_vec3;
}

IntervalMatrix mul5TensorVector(std::vector<std::vector<std::vector<std::vector<ibex::Function>>>> &d5f, IntervalVector x,
								IntervalVector vec1, IntervalVector vec2, IntervalVector vec3, IntervalVector vec4, int dim)
{
	IntervalMatrix d5f_vec1_i_j_k = ibex::Matrix::zeros(dim);
	IntervalMatrix d5f_vec1_vec2_i_j = ibex::Matrix::zeros(dim);
	IntervalMatrix d5f_vec1_vec2_vec3_i = ibex::Matrix::zeros(dim);
	IntervalMatrix d5f_vec1_vec2_vec3_vec4 = ibex::Matrix::zeros(dim);

	for (int i = 0; i < dim; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			for (int k = 0; k < dim; k++)
			{
				for (int l = 0; l < dim; l++)
				{
					d5f_vec1_i_j_k[l] = d5f[i][j][k][l].eval_matrix(x) * vec1;
				}
				d5f_vec1_vec2_i_j[k] = d5f_vec1_i_j_k * vec2;
			}
			d5f_vec1_vec2_vec3_i[j] = d5f_vec1_vec2_i_j * vec3;
		}
		d5f_vec1_vec2_vec3_vec4[i] = d5f_vec1_vec2_vec3_i * vec4;
	}

	return d5f_vec1_vec2_vec3_vec4;
}

/* eval LTE of RK1 computation */
IntervalMatrix evalLTE_RK1(Function &f, Function &fjacob, std::vector<ibex::Function> &d2f,
		IntervalVector &f_roughEnc, IntervalMatrix &fjacob_roughEnc,
		IntervalVector &roughEnc, IntervalMatrix gbox, int dim)
{

    IntervalMatrix d2f_f = mul2TensorVector(d2f, roughEnc, f_roughEnc, dim);
    IntervalMatrix df_df = fjacob_roughEnc * fjacob_roughEnc;

	IntervalMatrix LTE = d2f_f + df_df;
	LTE *= gbox; //C1-Lohner algorithm paper by Zgliczynski: Lemma 1 and especially equation [17]

	return LTE;
}

/* eval LTE of RK2 computation */
IntervalMatrix evalLTE_RK2(Function &f, Function &fjacob, std::vector<ibex::Function> &d2f,
						   std::vector<std::vector<ibex::Function>> &d3f, IntervalVector &f_roughEnc, IntervalMatrix &fjacob_roughEnc,
						   IntervalVector &roughEnc, IntervalMatrix gbox, int dim)
{

	//term 1
	IntervalMatrix d3f_f_f = mul3TensorVector(d3f, roughEnc, f_roughEnc, f_roughEnc, dim);

	//term 2
	IntervalMatrix d2f_f = mul2TensorVector(d2f, roughEnc, f_roughEnc, dim);

	//term 3
	IntervalVector df_f = fjacob_roughEnc * f_roughEnc;
	IntervalMatrix d2f_df_f = mul2TensorVector(d2f, roughEnc, df_f, dim);

	//term 4
	IntervalMatrix df_d2f_f = fjacob_roughEnc * d2f_f;

	//term 5
	IntervalMatrix df_df_df = fjacob_roughEnc * fjacob_roughEnc * fjacob_roughEnc;

	IntervalMatrix LTE = 0.25 * d3f_f_f;
	LTE += 0.5 * d2f_f * fjacob_roughEnc;
	LTE += d2f_df_f;
	LTE += df_d2f_f;
	LTE += df_df_df;
	LTE *= gbox; //C1-Lohner algorithm paper by Zgliczynski: Lemma 1 and especially equation [17]

	return LTE;
}

/* eval LTE of RK4 computation */
IntervalMatrix evalLTE_RK4 (Function &f, Function &fjacob, std::vector<ibex::Function> &d2f,
		std::vector<std::vector<ibex::Function>> &d3f, std::vector<std::vector<std::vector<ibex::Function>>> &d4f,
		std::vector<std::vector<std::vector<std::vector<ibex::Function>>>> &d5f,
		IntervalVector &f_xRoughEnc, IntervalMatrix &fjacob_xRoughEnc,
		IntervalVector &xRoughEnc, IntervalMatrix gbox, ibex::Matrix &idMatrix, int &dim)
{
	//LTE for x = xRoughEnc
	IntervalMatrix df_LTE = fjacob_xRoughEnc;
	IntervalMatrix F_LTE = idMatrix;

	IntervalMatrix d2f_f_LTE = mul2TensorVector(d2f, xRoughEnc, f_xRoughEnc, dim);
	IntervalMatrix d2f_df_f_LTE = mul2TensorVector(d2f, xRoughEnc, df_LTE*f_xRoughEnc, dim);
	IntervalMatrix d2f_d2f_f_f_LTE = mul2TensorVector(d2f, xRoughEnc, d2f_f_LTE*f_xRoughEnc, dim);
	IntervalMatrix d2f_d2f_df_f_f_LTE = mul2TensorVector(d2f, xRoughEnc, d2f_df_f_LTE*f_xRoughEnc, dim);
	IntervalMatrix d2f_d2f_f_df_f_LTE = mul2TensorVector(d2f, xRoughEnc, (d2f_f_LTE*df_LTE)*f_xRoughEnc, dim);
	IntervalMatrix d2f_df_df_f_LTE = mul2TensorVector(d2f, xRoughEnc, (df_LTE*df_LTE)*f_xRoughEnc, dim);
	IntervalMatrix d2f_df_df_df_f_LTE = mul2TensorVector(d2f, xRoughEnc, (df_LTE*df_LTE*df_LTE)*f_xRoughEnc, dim);
	IntervalMatrix d2f_df_d2f_f_f_LTE = mul2TensorVector(d2f, xRoughEnc, (df_LTE*d2f_f_LTE)*f_xRoughEnc, dim);

	IntervalMatrix d3f_f_f_LTE = mul3TensorVector(d3f, xRoughEnc, f_xRoughEnc, f_xRoughEnc, dim);
	IntervalMatrix d3f_df_f_f_LTE = mul3TensorVector(d3f, xRoughEnc, df_LTE*f_xRoughEnc, f_xRoughEnc, dim);
	IntervalMatrix d3f_df_f_df_f_LTE = mul3TensorVector(d3f, xRoughEnc, df_LTE*f_xRoughEnc, df_LTE*f_xRoughEnc, dim);
	IntervalMatrix d3f_df_df_f_f_LTE = mul3TensorVector(d3f, xRoughEnc, (df_LTE*df_LTE)*f_xRoughEnc, f_xRoughEnc, dim);
	IntervalMatrix d3f_f_df_f_LTE = mul3TensorVector(d3f, xRoughEnc, f_xRoughEnc, df_LTE*f_xRoughEnc, dim);
	IntervalMatrix d3f_f_df_df_f_LTE = mul3TensorVector(d3f, xRoughEnc, f_xRoughEnc, (df_LTE*df_LTE)*f_xRoughEnc, dim);
	IntervalMatrix d3f_d2f_f_f_f_LTE = mul3TensorVector(d3f, xRoughEnc, d2f_f_LTE * f_xRoughEnc, f_xRoughEnc, dim);
	IntervalMatrix d3f_f_d2f_f_f_LTE = mul3TensorVector(d3f, xRoughEnc, f_xRoughEnc, d2f_f_LTE * f_xRoughEnc, dim);

	IntervalMatrix d2f_d3f_f_f_f_LTE = mul2TensorVector(d2f, xRoughEnc, d3f_f_f_LTE*f_xRoughEnc, dim);

	IntervalMatrix d4f_f_f_f_LTE = mul4TensorVector(d4f, xRoughEnc, f_xRoughEnc, f_xRoughEnc, f_xRoughEnc, dim);
	IntervalMatrix d4f_df_f_f_f_LTE = mul4TensorVector(d4f, xRoughEnc, df_LTE*f_xRoughEnc, f_xRoughEnc, f_xRoughEnc, dim);
	IntervalMatrix d4f_f_df_f_f_LTE = mul4TensorVector(d4f, xRoughEnc, f_xRoughEnc, df_LTE*f_xRoughEnc, f_xRoughEnc, dim);
	IntervalMatrix d4f_f_f_df_f_LTE = mul4TensorVector(d4f, xRoughEnc, f_xRoughEnc, f_xRoughEnc, df_LTE*f_xRoughEnc, dim);

	IntervalMatrix d5f_f_f_f_f_LTE = mul5TensorVector(d5f, xRoughEnc, f_xRoughEnc, f_xRoughEnc, f_xRoughEnc, f_xRoughEnc, dim);


	// derivatives of df and F
	IntervalMatrix deriv1_df_LTE = d2f_f_LTE;
	IntervalMatrix deriv2_df_LTE = d3f_f_f_LTE + d2f_df_f_LTE;
	IntervalMatrix deriv3_df_LTE = d4f_f_f_f_LTE + d3f_df_f_f_LTE + d3f_f_df_f_LTE + d2f_d2f_f_f_LTE + d3f_f_df_f_LTE + d2f_df_df_f_LTE;
	IntervalMatrix deriv4_df_LTE = d5f_f_f_f_f_LTE + d4f_df_f_f_f_LTE + d4f_f_df_f_f_LTE + d3f_d2f_f_f_f_LTE + d3f_f_d2f_f_f_LTE + d2f_d3f_f_f_f_LTE;
	deriv4_df_LTE += d4f_f_df_f_f_LTE + d3f_df_df_f_f_LTE + d2f_d2f_df_f_f_LTE + d4f_f_f_df_f_LTE + d3f_df_f_df_f_LTE + d2f_d2f_f_df_f_LTE;
	deriv4_df_LTE += 2*d4f_f_f_df_f_LTE + 2*d3f_df_f_df_f_LTE + d3f_f_df_df_f_LTE + d2f_d2f_f_df_f_LTE;
	deriv4_df_LTE += 2*d3f_f_d2f_f_f_LTE + d2f_df_d2f_f_f_LTE + 2*d3f_f_df_df_f_LTE + d2f_df_df_df_f_LTE;

	IntervalMatrix deriv1_F_LTE = df_LTE;
	IntervalMatrix deriv2_F_LTE = deriv1_df_LTE + df_LTE*df_LTE;
	IntervalMatrix deriv3_F_LTE = deriv2_df_LTE + 2*deriv1_df_LTE*df_LTE + df_LTE * deriv2_F_LTE;
	IntervalMatrix deriv4_F_LTE = deriv3_df_LTE + 3*deriv2_df_LTE*df_LTE + 3*deriv1_df_LTE*deriv2_F_LTE + df_LTE*deriv3_F_LTE;


	//term LTE with xRoughEnc
	IntervalMatrix LTE = -1./2880 * deriv4_df_LTE * F_LTE;
	LTE += -1./720 * deriv3_df_LTE * deriv1_F_LTE;
	LTE += -1./480 * deriv2_df_LTE * deriv2_F_LTE;
	LTE += 1./480 * deriv1_df_LTE * deriv3_F_LTE;
	LTE += 1./720 * df_LTE * deriv4_F_LTE;
	LTE += -1./96 * deriv1_df_LTE * df_LTE * deriv2_F_LTE;
	LTE += -1./288 * df_LTE * df_LTE * deriv3_F_LTE;
	LTE += 1./96 * df_LTE * df_LTE * df_LTE * deriv2_F_LTE;
	LTE *= gbox; //C1-Lohner algorithm paper by Zgliczynski: Lemma 1 and especially equation [17]

	return LTE;
}

/* Algorithm for rough enclosure for ODE in section 2.4 of C1 Lohner algorithm paper
 * of P. Zgliczynski - similar to enclosure in FirstOrderEnclosure.hpp of CAPD */
IntervalVector computeRoughEnc(Function &f, Function &fjacob, IntervalVector &f_xbox,
		IntervalVector &xbox, double h, int dim, double eps = 1e-05)
{
	Interval trialStep(-0.01 * h, 1.01 * h);
	IntervalVector newEnc(dim), currentEnc(dim); //new enclosure and current enclosure of iterative procedure for equation (22)

	Interval step(0, h);
	//typename ScalarType::BoundType multf = 1.5;   // factor to multiply coordinates if inclusion fails

	IntervalMatrix fjacob_currentEnc(dim, dim);
	IntervalVector f_currentEnc(dim);

	currentEnc = xbox + trialStep * f_xbox + IntervalVector(dim, Interval(-eps, eps));

	bool found = false;
	int counter = 0,
		limit = 10 + 2 * dim; // maximum numbers of attempts to find enclosure

	while ((!found) && (counter < limit))
	{
		counter++;

		f_currentEnc = f.eval_vector(currentEnc);

		/* Tighter but slower approach:
		 *
		 * fjacob_currentEnc = fjacob.eval_matrix(currentEnc);
		 * f_currentEnc = evalMeanValue(f, fjacob_currentEnc, currentEnc);
		 */

		newEnc = xbox + step * f_currentEnc;
		if (newEnc.is_interior_subset(currentEnc))
		{
			found = true;
		}
		else
		{
			currentEnc = newEnc.inflate(0.01);
		}
	}

	if (!found)
	{
		throw std::runtime_error("Solver error: cannot find rough enclosure guaranteeing bounds");
	}
	return newEnc;
}

/* Algorithm for rough enclosure for the variational part in section 2.5 of C1 Lohner
 * algorithm paper of P. Zgliczynski - similar to jacEnclosure in FirstOrderEnclosure.hpp of CAPD */
IntervalMatrix computeRoughEncJac(IntervalMatrix &fjacob_roughEnc, ibex::Matrix &idMatrix, double h, int dim)
{
	Interval step(0, h);

	double l = computeIntEig(1. / 2 * (fjacob_roughEnc + fjacob_roughEnc.transpose()), dim); //computation of logarithmic norm for euclidian norm (Theorem 3 of Zgliczynski)

	/* Less tight but faster approaches:
	 *
	 * double l = rohnEigMax(1. / 2 * (fjacob_roughEnc + fjacob_roughEnc.transpose()), dim);
	 *
	 * or
	 *
	 * double l = cfEigMax(1. / 2 * (fjacob_roughEnc + fjacob_roughEnc.transpose()), dim);
	 */

	IntervalMatrix W(dim, dim, Interval(-1, 1) * exp(step * l)); // W in paper "C^1 - Lohner algorithm"

	IntervalMatrix W3(dim, dim); // W_3 in paper "C^1 - Lohner algorithm"

	W3 = idMatrix + step * fjacob_roughEnc * W;

	W3 &= W; // Intersection operation of ibex for IntervalMatrix

	return W3;
}

/* Evaluation #3 in the paper: choose a matrix U in A*B and perform QR decomposition of U */
void getQnextStep(ibex::Matrix U, ibex::Vector l, IntervalMatrix &Q, IntervalMatrix &Q_inv, int dim)
{

	VectorXd l_eigen = ibex_to_eigen(l, dim); //l is a vector containing the lengths of the edges of the solution set

	// copy row_sum to std::vector (vec_sum)
	std::vector<double> l_stdvec;
	l_stdvec.resize(dim);
	VectorXd::Map(&l_stdvec[0], dim) = l_eigen;

	// sorting
	std::vector<int> sort_idx;			 // Index: std::vector
	std::vector<double> l_stdvec_sorted; // N
	sort(l_stdvec, l_stdvec_sorted, sort_idx); //Sorting that gives Index
	// converting std::vec to Eigen::Vector
	VectorXi sort_idx_eigen = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(sort_idx.data(), dim);

	// change columns of U in decreasing order of the lengths in vector l
	int i, j;
	ibex::Matrix U_copy = U;
	for (i = 0; i < dim; i++)
	{
		for (j = 0; j < dim; j++)
		{
			U[i][j] = U_copy[i][sort_idx_eigen[dim - 1 - j]]; //because sort_idx is in increasing order
		}
	}

	MatrixXd U_eigen = ibex_to_eigen(U, dim);

	HouseholderQR<MatrixXd> qr(U_eigen);

	MatrixXd Q_eigen = qr.householderQ();

	Q = eigen_to_ibex(Q_eigen, dim);

	Q_inv = gaussInverseMatrix(Q); //TODO compare with Inverse in the Eigen LU module
}


//Euler method to integrate an interval
IntervalVector RK_center(Function &f, Function &fjacob, IntervalVector &f_cx, IntervalVector &cx, double h, int dim)
{
    IntervalVector cxRoughEnc = computeRoughEnc(f, fjacob, f_cx, cx, h, dim, 1e-16);

    IntervalVector cx_new = cx + Interval(h) * f_cx + Interval(h*h) / 2 * evalSecondTaylor(f, fjacob, cxRoughEnc);

	return cx_new;
}

//RK2 integration of the variational equation
IntervalMatrix computeRK2variationalQR(Function &f, Function &fjacob, std::vector<ibex::Function> &d2f,
									   std::vector<std::vector<ibex::Function>> &d3f, IntervalMatrix &fjacob_xbox,
									   IntervalVector &f_xbox, IntervalVector &xbox, IntervalMatrix &gbox, double h, int dim,
									   IntervalMatrix &gDeltaQ, ibex::Matrix &gPointQ, IntervalMatrix &Qg,
									   IntervalMatrix &Qg_inv, IntervalMatrix usedC, ibex::Matrix &idMatrix, int variablesDim)
{
	IntervalVector xRoughEnc = computeRoughEnc(f, fjacob, f_xbox, xbox, h, dim); //W1 in the paper

	IntervalMatrix fjacob_xRoughEnc = fjacob.eval_matrix(xRoughEnc);
	IntervalVector f_xRoughEnc = evalMeanValue(f, fjacob_xRoughEnc, xRoughEnc);

	IntervalMatrix gRoughEnc = computeRoughEncJac(fjacob_xRoughEnc, idMatrix, h, dim).submatrix(0, variablesDim - 1, 0, variablesDim - 1); //W3 in the paper

	IntervalVector xbox_05h = RK_center(f, fjacob, f_xbox, xbox, 0.5 * h, dim); //xbox + 0.5 * h * f_xbox + LTE

	IntervalMatrix k1 = fjacob_xbox;
	IntervalMatrix k2 = fjacob.eval_matrix(xbox_05h) * (idMatrix + 0.5 * h * k1);

	IntervalMatrix gEncJk2 = (idMatrix + h * k2).submatrix(0, variablesDim - 1, 0, variablesDim - 1);

	IntervalMatrix LTE = evalLTE_RK2(f, fjacob, d2f, d3f, f_xRoughEnc, fjacob_xRoughEnc, xRoughEnc, idMatrix, dim);

	IntervalMatrix gEncJ = gEncJk2 + h * h * h * 1. / 6. * LTE.submatrix(0, variablesDim - 1, 0, variablesDim - 1)
			* gRoughEnc; //J_k in the paper - rough enclosure of the variational part with starting value gbox = Id

	IntervalMatrix Q(variablesDim, variablesDim), Q_inv(variablesDim, variablesDim);

	ibex::Matrix U = (gEncJ * Qg).mid();

	ibex::Vector length_edges(variablesDim); //used to reorder U before QR decomposition

	int i;
	for (i = 0; i < variablesDim; i++)
	{
		length_edges[i] = norm(U.col(i)) * norm(gDeltaQ.row(i).diam());
	}

	getQnextStep(U, length_edges, Q, Q_inv, variablesDim);

	IntervalMatrix J_gPointQ = gEncJ * IntervalMatrix(gPointQ); // J_k times gPoint
	gPointQ = J_gPointQ.mid();									//this is V_{k+1}

	//Delta\tilde{V_{k+1}}          tilde{Delta V_k}        Z_{k+1}
	gDeltaQ = (Q_inv * gEncJ * Qg) * gDeltaQ + Q_inv * (J_gPointQ - gPointQ);

	Qg = Q;

	Qg_inv = Q_inv;

	IntervalMatrix usedC_temp = usedC.submatrix(0, variablesDim - 1, 0, variablesDim - 1);

	// simpler formula for [ V_{k+1} = gbox_new]
	IntervalMatrix gbox_usedC = (usedC_temp * gEncJ) * gbox;
	gbox_usedC &= (usedC_temp * Qg) * gDeltaQ + usedC_temp * gPointQ; //intersection between directly computed version and version with QR method

	gbox = gEncJ * gbox;			// simpler formula for next step gbox
	gbox &= Qg * gDeltaQ + gPointQ; //intersection between directly computed version and version with QR method

	return gbox_usedC; //V_{k+1}
}

//RK4 integration of the variational equation
IntervalMatrix computeRK4variationalQR(Function &f, Function &fjacob, std::vector<ibex::Function> &d2f,
									   std::vector<std::vector<ibex::Function>> &d3f,
									   std::vector<std::vector<std::vector<ibex::Function>>> &d4f,
									   std::vector<std::vector<std::vector<std::vector<ibex::Function>>>> &d5f,
									   IntervalMatrix &fjacob_xbox, IntervalVector &f_xbox, IntervalVector &xbox,
									   IntervalMatrix &gbox, double h, int dim,
									   IntervalMatrix &gDeltaQ, ibex::Matrix &gPointQ, IntervalMatrix &Qg, IntervalMatrix &Qg_inv,
									   IntervalMatrix usedC, ibex::Matrix &idMatrix, int variablesDim)
{
	IntervalVector xRoughEnc = computeRoughEnc(f, fjacob, f_xbox, xbox, h, dim); //W1 in the paper

	IntervalMatrix fjacob_xRoughEnc = fjacob.eval_matrix(xRoughEnc);
	IntervalVector f_xRoughEnc = evalMeanValue(f,fjacob_xRoughEnc,xRoughEnc);

	IntervalMatrix gRoughEnc = computeRoughEncJac(fjacob_xRoughEnc, idMatrix, h, dim).submatrix(0,variablesDim-1,0,variablesDim-1); //W3 in the paper

	IntervalVector xbox_05h(dim), xbox_h(dim);
	xbox_05h = RK_center(f, fjacob, f_xbox, xbox, 0.5 * h, dim); //xbox + 0.5 * h * f_xbox + LTE
	xbox_h = RK_center(f, fjacob, f_xbox, xbox, h, dim); //xbox + h * f_xbox + LTE

	IntervalMatrix k1 = fjacob_xbox;
    IntervalMatrix k2 = fjacob.eval_matrix(xbox_05h) * (idMatrix + 0.5 * h * k1);
    IntervalMatrix k3 = fjacob.eval_matrix(xbox_05h) * (idMatrix + 0.5 * h * k2);
    IntervalMatrix k4 = fjacob.eval_matrix(xbox_h) * (idMatrix + h * k3);

    IntervalMatrix gEncJk2 = (idMatrix + h * 1./6 * (k1 + 2*k2 + 2*k3 + k4)).submatrix(0,variablesDim-1,0,variablesDim-1);

    IntervalMatrix LTE = evalLTE_RK4(f, fjacob, d2f, d3f, d4f, d5f, f_xRoughEnc, fjacob_xRoughEnc, xRoughEnc, idMatrix, idMatrix, dim);

    IntervalMatrix gEncJ = gEncJk2 + pow((double)h,5.) * LTE.submatrix(0,variablesDim-1,0,variablesDim-1)
    		* gRoughEnc; //J_k in the paper - rough enclosure of the variational part with starting value gbox = Id

	IntervalMatrix Q(variablesDim, variablesDim), Q_inv(variablesDim, variablesDim);

	ibex::Matrix U = (gEncJ * Qg).mid();

	ibex::Vector length_edges(variablesDim);//used to reorder U before QR decomposition

	int i;
	for(i=0; i<variablesDim; i++){
		length_edges[i] = norm(U.col(i)) * norm(gDeltaQ.row(i).diam());
	}

    getQnextStep(U, length_edges, Q, Q_inv, variablesDim);

	IntervalMatrix J_gPointQ = gEncJ * IntervalMatrix(gPointQ); // J_k times gPoint
    gPointQ = J_gPointQ.mid(); //this is V_{k+1}

    //Delta\tilde{V_{k+1}}          tilde{Delta V_k}        Z_{k+1}
	gDeltaQ = (Q_inv * gEncJ * Qg) * gDeltaQ + Q_inv * (J_gPointQ - gPointQ);

	Qg = Q;

	Qg_inv = Q_inv;

	IntervalMatrix usedC_temp = usedC.submatrix(0,variablesDim-1,0,variablesDim-1);

    // simpler formula for [ V_{k+1} = gbox_new]
    IntervalMatrix gbox_usedC = (usedC_temp * gEncJ) * gbox;
	gbox_usedC &= (usedC_temp * Qg) * gDeltaQ + usedC_temp * gPointQ;//intersection between directly computed version and version with QR method

	gbox = gEncJ * gbox; // simpler formula for next step gbox
	gbox &= Qg * gDeltaQ + gPointQ; //intersection between directly computed version and version with QR method

	return gbox_usedC; //V_{k+1}
}

//RK1 (=Euler) integration of the variational equation
IntervalMatrix computeRK1variationalQR(Function& f, Function& fjacob, std::vector<ibex::Function>& d2f,
		IntervalMatrix& fjacob_xbox, IntervalVector& f_xbox, IntervalVector& xbox, IntervalMatrix& gbox, double h, int dim,
		IntervalMatrix& gDeltaQ, ibex::Matrix& gPointQ, IntervalMatrix& Qg, IntervalMatrix& Qg_inv, IntervalMatrix usedC, ibex::Matrix& idMatrix, int variablesDim){

	IntervalVector xRoughEnc = computeRoughEnc(f, fjacob, f_xbox, xbox, h, dim); //W1 in the paper

	IntervalMatrix fjacob_xRoughEnc = fjacob.eval_matrix(xRoughEnc);
	IntervalVector f_xRoughEnc = evalMeanValue(f,fjacob_xRoughEnc,xRoughEnc);

	IntervalMatrix gRoughEnc = computeRoughEncJac(fjacob_xRoughEnc, idMatrix, h, dim).submatrix(0,variablesDim-1,0,variablesDim-1); //W3 in the paper

	IntervalMatrix k1 = fjacob_xbox;

    IntervalMatrix gEncJk2 = (idMatrix + h * k1).submatrix(0,variablesDim-1,0,variablesDim-1);

    IntervalMatrix LTE = evalLTE_RK1(f, fjacob, d2f, f_xRoughEnc, fjacob_xRoughEnc, xRoughEnc, idMatrix, dim);

    IntervalMatrix gEncJ = gEncJk2 + 0.5 * h*h * LTE.submatrix(0,variablesDim-1,0,variablesDim-1)
    		* gRoughEnc; //J_k in the paper - rough enclosure of the variational part with starting value gbox = Id

	IntervalMatrix Q(variablesDim, variablesDim), Q_inv(variablesDim, variablesDim);

	ibex::Matrix U = (gEncJ * Qg).mid();

	ibex::Vector length_edges(variablesDim);//used to reorder U before QR decomposition

	for(int i=0; i<variablesDim; i++){
		length_edges[i] = norm(U.col(i)) * norm(gDeltaQ.row(i).diam());
	}

    getQnextStep(U, length_edges, Q, Q_inv, variablesDim);

	IntervalMatrix J_gPointQ = gEncJ * IntervalMatrix(gPointQ); // J_k times gPoint
    gPointQ = J_gPointQ.mid(); //this is V_{k+1}

    //Delta\tilde{V_{k+1}}          tilde{Delta V_k}        Z_{k+1}
	gDeltaQ = (Q_inv * gEncJ * Qg) * gDeltaQ + Q_inv * (J_gPointQ - gPointQ);

	Qg = Q;

	Qg_inv = Q_inv;

	IntervalMatrix usedC_temp = usedC.submatrix(0,variablesDim-1,0,variablesDim-1);

    // simpler formula for [ V_{k+1} = gbox_new]
    IntervalMatrix gbox_usedC = (usedC_temp * gEncJ) * gbox;
	gbox_usedC &= (usedC_temp * Qg) * gDeltaQ + usedC_temp * gPointQ;//intersection between directly computed version and version with QR method

	gbox = gEncJ * gbox; // simpler formula for next step gbox
	gbox &= Qg * gDeltaQ + gPointQ; //intersection between directly computed version and version with QR method

	return gbox_usedC; //V_{k+1}
}



IntervalMatrix computeRK1variationalQR_reinitGbox(Function& f, Function& fjacob, std::vector<ibex::Function>& d2f,
		IntervalMatrix& fjacob_xbox, IntervalVector& f_xbox, IntervalVector& xbox, IntervalMatrix &gbox, double h, int dim,
		IntervalMatrix usedC, ibex::Matrix& idMatrix, int variablesDim){

	ibex::Matrix idMatrixVariables = ibex::Matrix::eye(variablesDim);

	IntervalMatrix gDeltaC_sigma = ibex::Matrix::zeros(variablesDim);
    ibex::Matrix gPointC_sigma = idMatrixVariables;
    IntervalMatrix Qg_sigma = idMatrixVariables;
    IntervalMatrix Qg_inv_sigma = idMatrixVariables;

    return computeRK1variationalQR(f, fjacob, d2f, fjacob_xbox, f_xbox, xbox, gbox, h, dim,
    		gDeltaC_sigma, gPointC_sigma, Qg_sigma, Qg_inv_sigma, usedC, idMatrix, variablesDim);
}



IntervalMatrix computeRK2variationalQR_reinitGbox(Function& f, Function& fjacob, std::vector<ibex::Function>& d2f,
		std::vector<std::vector<ibex::Function>>& d3f,
		IntervalMatrix& fjacob_xbox, IntervalVector& f_xbox, IntervalVector& xbox, IntervalMatrix &gbox, double h, int dim,
		IntervalMatrix usedC, ibex::Matrix& idMatrix, int variablesDim){

	ibex::Matrix idMatrixVariables = ibex::Matrix::eye(variablesDim);

	IntervalMatrix gDeltaC_sigma = ibex::Matrix::zeros(variablesDim);
    ibex::Matrix gPointC_sigma = idMatrixVariables;
    IntervalMatrix Qg_sigma = idMatrixVariables;
    IntervalMatrix Qg_inv_sigma = idMatrixVariables;

    return computeRK2variationalQR(f, fjacob, d2f, d3f, fjacob_xbox, f_xbox, xbox, gbox, h, dim,
    		gDeltaC_sigma, gPointC_sigma, Qg_sigma, Qg_inv_sigma, usedC, idMatrix, variablesDim);
}



IntervalMatrix computeRK4variationalQR_reinitGbox(Function& f, Function& fjacob, std::vector<ibex::Function>& d2f,
		std::vector<std::vector<ibex::Function>>& d3f, std::vector<std::vector<std::vector<ibex::Function>>>& d4f,
		std::vector<std::vector<std::vector<std::vector<ibex::Function>>>>& d5f,
		IntervalMatrix& fjacob_xbox, IntervalVector& f_xbox, IntervalVector& xbox, IntervalMatrix &gbox, double h, int dim,
		IntervalMatrix usedC, ibex::Matrix& idMatrix, int variablesDim){

	ibex::Matrix idMatrixVariables = ibex::Matrix::eye(variablesDim);

	IntervalMatrix gDeltaC_sigma = ibex::Matrix::zeros(variablesDim);
    ibex::Matrix gPointC_sigma = idMatrixVariables;
    IntervalMatrix Qg_sigma = idMatrixVariables;
    IntervalMatrix Qg_inv_sigma = idMatrixVariables;

    return computeRK4variationalQR(f, fjacob, d2f, d3f, d4f, d5f, fjacob_xbox, f_xbox, xbox, gbox, h, dim,
    		gDeltaC_sigma, gPointC_sigma, Qg_sigma, Qg_inv_sigma, usedC, idMatrix, variablesDim);
}

IntervalMatrix linear_grad(IntervalMatrix& fjacob_cx, double h, int dim){

	//e^{At}*F0
	//F0 = id because for the metric we suppose that we are reinitializing every timestep
	MatrixXd fJacEigen = ibex_to_eigen(fjacob_cx.mid(), dim);
	MatrixXd exp_fJacEigen = (fJacEigen * h).exp();
	IntervalMatrix exp_fJac = eigen_to_ibex(exp_fJacEigen, dim);

	return exp_fJac;
}

ibex::Interval get_SF(IntervalMatrix CgCi, int dim)
{
	IntervalMatrix sq_CgCi = CgCi.transpose() * CgCi;
	double eig_CgCi = computeIntEig(sq_CgCi, dim);
	return sqrt(eig_CgCi);
}

#endif /* COMPUTE_NEXT_STEP_H */
