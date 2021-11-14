/*
 * ibex_api.h
 *
 ** Lagrangian Reachtubes: The next Generation
 *
 *      Authors: Sophie Gruenbacher, Md Ariful Islam and Jacek Cyranka
 *      Contact: sophie.gruenbacher@tuwien.ac.at
 */

#ifndef IBEX_API_H
#define IBEX_API_H

/* Standard C++ library*/
#include <unistd.h>
#include <iomanip>
#include <string.h>
#include <stdio.h>
#include <sstream>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <iostream>
#include <vector>
#include <functional>
#include <numeric>
#include <set>
#include <iterator>

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
using namespace std;

/* std::vector<double> -> ibex::Vector*/
ibex::Vector std2ibexVec(std::vector<double> std_vec, int begin, int end){

    int dim = end - begin;
    ibex::Vector ibex_vec(dim);
    for(int i=begin; i<end; i++)
        ibex_vec[i] = std_vec[i];
    return ibex_vec;    
}

/* ibex::Vector -> std::vector<double>*/
std::vector<double> ibex2stdVec(ibex::Vector ibex_vec){
    int dim = ibex_vec.size();
    std::vector<double> std_vec(dim);
    for(int i=0; i<dim; i++)
        std_vec[i] = ibex_vec[i];
    return std_vec;    
}

/* IBEX: Matrix -> IntervalMatrix */
ibex::IntervalMatrix d2iMatrix(ibex::Matrix &m){
    const int R = m.nb_rows();
    const int C = m.nb_cols();
    
    ibex::IntervalMatrix r(R,C);
    for(int i=0; i < R; i++){
        for(int j=0; j < C; j++){
            r[i][j] = ibex::Interval(m[i][j],m[i][j]);
        }
    }
    return r;
}

/* IBEX: IntervalMatrix -> ibex::Matrix  */
ibex::Matrix ibex_lMatrix(ibex::IntervalMatrix iA){
    int R = iA.nb_rows();
    int C = iA.nb_cols();
    ibex::Matrix r(R,C);
    for(int i=0; i<R; i++){
        for(int j=0; j<C; j++)
            r[i][j] = iA[i][j].lb();
    }
    return r;
}

/* IBEX: IntervalMatrix -> ibex::Matrix  */
ibex::Matrix ibex_rMatrix(ibex::IntervalMatrix iA){
    int R = iA.nb_rows();
    int C = iA.nb_cols();
    ibex::Matrix r(R,C);
    for(int i=0; i<R; i++){
        for(int j=0; j<C; j++)
            r[i][j] = iA[i][j].ub();
    }
    return r;
}

/* Matlab's all(in0(A,B)): checking inclusion of iA in iB componentwise*/
bool  all_in(ibex::IntervalVector iA, ibex::IntervalVector iB, int dim){
    for(int i=0; i< dim; i++){
        if(!(iB[i].lb() < iA[i].lb() && iA[i].ub() < iB[i].ub()))
            return false;
    }
    return true;
}

/*Matlab's B = flipud(A) method -- non-rigorous: flip upside down or reverse in case of vector */
ibex::IntervalVector flipud(ibex::IntervalVector vecA, int dim){
    ibex::IntervalVector vecB(dim);
    for(int i = 0; i<dim; i++){
        vecB[i] = vecA[dim-1-i];
    }
   return vecB;
}

ibex::IntervalMatrix getIMatrix(ibex::Matrix lm, ibex::Matrix rm){
    int R = rm.nb_rows();
    int C = rm.nb_cols();

    ibex::IntervalMatrix im(R,C);
    for(int i=0; i<R; i++){
        for(int j=0; j<C; j++)
            im[i][j] = ibex::Interval(lm[i][j], rm[i][j]);
    }
    return im;
}


ibex::IntervalMatrix deleteRowCol(ibex::IntervalMatrix& m, int idx){
    const int R = m.nb_rows();
    const int C = m.nb_cols();
    
    ibex::IntervalMatrix r(R-1,C-1);
    int cCount, rCount;
    rCount = 0;
    for(int i=0; i < R; i++){
        if(i == idx){
            continue;
        }
        cCount = 0;
        for(int j=0; j < C; j++){
            if(j == idx){
                continue;
            }
            r[rCount][cCount] = m[i][j];
            cCount++;
        }
        rCount++;
    }
    return r;      
}

ibex::IntervalMatrix getRowsCols(ibex::IntervalMatrix& m, vector<int> idxVec){
    //const int R = m.nb_rows();
    //const int C = m.nb_cols();
    
    int nDim = idxVec.size();
    ibex::IntervalMatrix r(nDim,nDim);
    int cCount, rCount;
    rCount = 0;
    for(int i: idxVec){
        cCount = 0;
        for(int j: idxVec){
            r[rCount][cCount] = m[i][j];
            cCount++;
        }
        rCount++;
    }
    return r;
        
}
// ivec = ivec * intv
ibex::IntervalVector ibex_mult(ibex::IntervalVector ivec, ibex::Interval intv){
    int dim = ivec.size();
    IntervalVector result(dim);
    for(int i=0; i<dim; i++)
        result[i] = ivec[i]*intv;
    return result;
}

// ivec = ivec + intv
ibex::IntervalVector ibex_sum(ibex::IntervalVector ivec, ibex::Interval intv){
    int dim = ivec.size();
    IntervalVector result(dim);
    for(int i=0; i<dim; i++)
        result[i] = ivec[i]+intv;
    return result;
}

ibex::Matrix deleteRowCol(ibex::Matrix& m, int idx){
    const int R = m.nb_rows();
    const int C = m.nb_cols();
    
    ibex::Matrix r(R-1,C-1);
    int cCount, rCount;
    rCount = 0;
    for(int i=0; i < R; i++){
        if(i == idx){
            continue;
        }
        cCount = 0;
        for(int j=0; j < C; j++){
            if(j == idx){
                continue;
            }
            r[rCount][cCount] = m[i][j];
            cCount++;
        }
        rCount++;
    }
    return r;    
}

void printIMatrix(IntervalMatrix A, std::ostream &fout, bool isDim = false){
	int rows = A.nb_rows();
	int cols = A.nb_cols();

	// print lower part
	for(int i=0; i<rows; i++){
		fout << "\t";
		for(int j=0; j<cols; j++){
			fout << A[i][j].lb() << " ";
		}
		fout << endl;
	}
	fout << endl;
	// print upper part
	for(int i=0; i<rows; i++){
		fout << "\t";
		for(int j=0; j<cols; j++){
			fout << A[i][j].ub()<< " ";
		}
		fout << endl;
	}

	// print diameter
	if(isDim) {
		fout << endl << "\tDiameter: " << endl;
		for(int i=0; i<rows; i++){
			fout << "\t";
			for(int j=0; j<cols; j++){
				fout  << A[i][j].ub() - A[i][j].lb()<< " ";
			}
			fout << endl;
		}
	}
	fout << endl;
}


void printDMatrix(ibex::Matrix A, std::ostream &fout){
	int rows = A.nb_rows();
	int cols = A.nb_cols();

	// print lower part
	for(int i=0; i<rows; i++){
		fout << "\t";
		for(int j=0; j<cols; j++){
			fout << A[i][j] << " ";
		}
		fout << endl;
	}
	fout << endl;
}

#endif /* IBEX_API_H */
