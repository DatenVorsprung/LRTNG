/*
 * computeIntEig.h
 *
 ** Lagrangian Reachtubes: The next Generation
 *
 *      Authors: Sophie Gruenbacher and Md Ariful Islam
 *      Contact: sophie.gruenbacher@tuwien.ac.at
 */

#ifndef COMPUTE_INTEIG_H
#define COMPUTE_INTEIG_H

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

/* Library defined by myself*/
#include "ibex_api.h"
#include "ibex_to_eigen_api.h"

/*Define namespace*/
using namespace Eigen;
using namespace ibex;
using namespace std;

/* Define some GLOBAL VARIABLES*/
#define REAL_MIN 2.225073858507201e-308
#define MMAX 100
#define INF 1.7976931348623157E+308

bool _DEBUG; // only important
bool _DEBUG_; // many more

////////////////// WRAPPER for STD::IOTA //////////////////////////////////////
template < class T > struct  IotaWrapper
{
    typedef T type;
    typedef std::function<type(const type&)> IncrFunction;

    type value;
    IncrFunction incrFunction;
    IotaWrapper () = delete;
    IotaWrapper(const type& n, const IncrFunction incrFunction): value(n), incrFunction(incrFunction){};
    operator type() { return value; }
    IotaWrapper& operator++ () { value = incrFunction(value); return *this;}
};

////////////////// Begin of sorting with index map /////////////////////////////////////
/* Comparison struct used by sort method*/
template <class T> struct index_cmp
{
    index_cmp(const T arr) : arr(arr) {}
    bool operator()(const int a, const int b) const
    {
        return arr[a] < arr[b];
    }
    const T arr;
};

/* This implementation is O(n), but also uses O(n) extra memory*/
template< class T >
void reorder(
        std::vector<T> & unordered,
        std::vector<int> const & index_map,
        std::vector<T> & ordered)
{
    // copy for the reorder according to index_map, because unsorted may also be
    // sorted
    std::vector<T> copy = unordered;
    ordered.resize(index_map.size());
    for(int i = 0; i<index_map.size();i++)
    {
        ordered[i] = copy[index_map[i]];
    }
}

/* Sorting that gives Index
Same as Matlab's: [B, I] = sort(A)
*/
template <class T>
void sort(
        std::vector<T> &unsorted,
        std::vector<T> &sorted,
        std::vector<int> &index_map)
{
    // Original unsorted index map
    index_map.resize(unsorted.size());
    for(size_t i=0;i<unsorted.size();i++)
    {
        index_map[i] = i;
    }
    // Sort the index map, using unsorted for comparison
    sort(
            index_map.begin(),
            index_map.end(),
            index_cmp<std::vector<T>& >(unsorted));

    sorted.resize(unsorted.size());
    reorder(unsorted,index_map,sorted);
}

/*Matlab's [eigVal, eigVec] = eig(A) method -- non-rigorous*/
void getEigen(MatrixXd A, VectorXd& eigVal, MatrixXd& eigVec){

    SelfAdjointEigenSolver<MatrixXd> eigensolver(A); //TO DO compute only eigenvalues
    if (eigensolver.info() != Success) abort();
    eigVal = eigensolver.eigenvalues();
    eigVec = eigensolver.eigenvectors();
}

////////////// BEGIN: Rump's Method for point Matrix -- Rigorous ////////////////////
void verifyeig(MatrixXd A, double lam, VectorXd rvec, ibex::Interval &iLam, ibex::IntervalVector &iVec, int dim){
    // consider rvecs contains single eigenvector corresponding to eigen value lam
    /* sorting the element in rvecs */
    int rows = dim;
    int cols = 1;
    VectorXd row_sum = rvec;

    double initial_lam = lam;

    // copy row_sum to std::vector (vec_sum)
    std::vector<double> vec_sum;
    vec_sum.resize(row_sum.size());
    VectorXd::Map(&vec_sum[0], row_sum.size()) = row_sum;

    // sorting [N,I] = sort(sum(abs(xs),2));
    std::vector<int> sort_idx; // Index: std::vector
    std::vector<double> sorted_sum; // N
    sort(vec_sum,sorted_sum,sort_idx);
    // converting std::vec to Eigen::Vector
    VectorXi idx_vec = Eigen::Map<Eigen::VectorXi,Eigen::Unaligned>(sort_idx.data(), sort_idx.size());

    // u = I(1:n-k)
    VectorXi u_idx = idx_vec.head(rows-cols);
    // v = I(n-k+1:n)
    VectorXi v_idx = idx_vec.tail(cols);
    MatrixXd midA = A; // A is point matrix, so midA is same as A

    // R = midA - lambda*speye(n);
    MatrixXd R = midA - lam*MatrixXd::Identity(dim, dim);

    //R(:,v) = -xs;
    R.col(v_idx(0)) = -rvec;

    //y = R\(midA*xs-lambda*xs);
    VectorXd y = R.colPivHouseholderQr().solve(midA*rvec-lam*rvec);

    //xs(u,:) = xs(u,:) - y(u,:);
    for(int i=0; i< rows-cols; i++){
        rvec(u_idx(i)) = rvec(u_idx(i)) - y(u_idx(i));
    }

    //lambda = lambda - sum(diag(y(v,:)))/k;
    /******** From here, we assume k = 1 */
    lam = lam - y(v_idx(0));
    //R = midA - lambda*speye(n);
    R = midA - lam*MatrixXd::Identity(dim, dim);
    //R(:,v) = -xs;
    R.col(v_idx(0)) = -rvec;
    R = R.inverse();

    //C = A - intval(lambda)*speye(n);
    iLam = ibex::Interval(lam,lam);

    // convert A to iA
    ibex::Matrix dA = eigen_to_ibex(A, dim);
    ibex::IntervalMatrix iA = d2iMatrix(dA);
    ibex::IntervalMatrix iC = iA - iLam*ibex::Matrix::eye(dim);

    //Z = - R * ( C * xs ); rowsXcols matrix
    // convert R to IMatrix
    ibex::Matrix dR = eigen_to_ibex(R, dim);
    ibex::IntervalMatrix iR = d2iMatrix(dR);

    // rvecs to IMatrix
    ibex::IntervalVector iRvec = eigen_to_ibex(rvec, dim);
    ibex::IntervalVector iZ = - iR * (iC * iRvec);

    // C(:,v) = -xs;
    for(int i=0; i < rows; i++){
        iC[i][v_idx[0]] = -iRvec[i];
    }

    //  C = speye(n) - R * C;
    iC = ibex::Matrix::eye(rows) - iR * iC;

    //  Y = Z;
    ibex::IntervalVector iY = iZ;

    //Eps = 0.1*mag(Y)*hull(-1,1) + midrad(0,realmin);
    ibex::IntervalVector iYmag = iY.mag();
    ibex::IntervalVector iEps = ibex_sum(ibex_mult(iYmag, ibex::Interval(-0.1,0.1)), ibex::Interval(-REAL_MIN, REAL_MIN));

    bool ready = false;
    int m = 0;
    ibex::IntervalVector iX(dim);
    ibex::IntervalVector iXX(dim);
    std::vector<bool> isSubset;
    while(!ready && m < MMAX){
        m++;
        iX = iY + iEps;
        iXX = iX;
        iXX[v_idx[0]] = 0.0;
        // Y = Z + C*X + R*(XX*X(v,:));
        //ibex::Interval iCX = iC*iX;
        ibex::IntervalVector iXX_temp = ibex_mult(iXX, iX[v_idx[0]]);
        iY = iZ + iC*iX + iR*iXX_temp;
        //ready = all(all(in0(Y,X)));
        ready = all_in(iY,iX, dim);
    }

    if(ready){ // inclusion found
        ibex::Interval rad = iY[v_idx[0]].mag(); // Eigenvalue Correction
        iLam = iLam + ibex::Interval(-rad.ub(), rad.ub());
        //Y(v,:) = 0;
        iY[v_idx[0]] = 0;
        iVec = iRvec + iY;

    } else {
        cout << "Verified eigenvalue inclusion not found using Rump's method!!" << endl;
        /*Test to use an alternative inclusion if Rump's method is not found, just to be able to run the code*/
        double eps = 1e-06;
        iVec = ibex_sum(eigen_to_ibex(rvec, dim),ibex::Interval(-eps,eps));
        iLam = ibex::Interval(initial_lam-eps, initial_lam+eps);
        cout << "initial_lam: " << initial_lam << endl;
        cout << "iLam: " << iLam << endl;
        //exit(1);
    }
}

// Courant Fischer Theorem
// λ_m(iA) ≤ λ_m(mag(iA)).
double cfEigMax(ibex::IntervalMatrix iA, int dim) {//TODO why do we need eigenvectors?
    double lamMax;
    ibex::Matrix magA = iA.mag();
    MatrixXd A = ibex_to_eigen(magA, dim);

    // non-rigourous eigenvalue and eigenvector
    VectorXd evalA;
    MatrixXd evecA;
    getEigen(A,evalA, evecA);

    ibex::Interval i_evalA;
    ibex::IntervalVector i_evecA(dim);
    //TODO do we need pairs of eigenvectors for complex eigenvalues?
    if(evalA(dim-1) != 0)
        verifyeig(A, evalA(dim-1), evecA.col(dim-1), i_evalA, i_evecA, dim);
    else
        i_evalA = 0.0;

    if(_DEBUG_) {
        cout << endl << "cf w.o. verify:	" << evalA(dim-1) << endl << endl;
    }

    lamMax = i_evalA.ub();
    return lamMax;

    //return evalA(dim-1); //non-rigorous
}


////////////// BEGIN: Rohn's Method for Interval Matrix -- Rigorous ////////////////////
Interval rohnEig(MatrixXd Ac, MatrixXd Ad, int dim){
    // compute eigenvalue of Ad
    VectorXd evalAc;
    MatrixXd evecAc;
    getEigen(Ac,evalAc, evecAc); // non-rigourous eigenvalue and eigenvector [V,D] = eig(A)
    ibex::Interval i_evalAc; // need to compute spectral radius, so all eigenvalues are needed
    ibex::IntervalVector i_evecAc(dim);

    // only maximum rigorous value is computed: max is in (dim-1)-th index
    if(evalAc(dim-1) != 0)
        verifyeig(Ac, evalAc(dim-1), evecAc.col(dim-1), i_evalAc, i_evecAc, dim);
    else
        i_evalAc = 0.0;

    if(_DEBUG_){
        cout << endl << "Ac: " << endl << Ac << endl;
        cout << "Max Eigenvalue of Ac: " << endl;
        cout << "[" << i_evalAc.lb() << ", " << i_evalAc.ub() << "]" << endl;
    }

    // compute eigenvalue of Ad
    VectorXd evalAd;
    MatrixXd evecAd;
    getEigen(Ad,evalAd, evecAd); // non-rigourous eigenvalue and eigenvector
    // Maximum eigVal =
    ibex::Interval i_evalAd; // need to compute spectral radius, so all eigenvalues are needed
    ibex::IntervalVector i_evecAd(dim);

    if(evalAd(dim-1) != 0)
        verifyeig(Ad, evalAd(dim-1), evecAd.col(dim-1), i_evalAd, i_evecAd, dim);
    else
        i_evalAd = 0.0;

    // Minimum eigVal =
    ibex::Interval i_evalAd_min; // need to compute spectral radius, so all eigenvalues are needed
    ibex::IntervalVector i_evecAd_min(dim);

    if(evalAd(0) != 0)
        verifyeig(Ad, evalAd(0), evecAd.col(0), i_evalAd_min, i_evecAd_min, dim);
    else
        i_evalAd_min = 0.0;

    if(_DEBUG_) {
        cout << endl << "Ad: " << endl << Ad << endl;
        cout << "Max Eigenvalue of Ad: " << endl;
        for(int i=0; i<dim; i++){
            cout << "[" << i_evalAd.lb() << ", " << i_evalAd.ub() << "]" << endl;
        }
    }

    // compute spectral of i_evalAd_vec
    // spectral = max([abs(inf(iEigsAd)); abs(sup(iEigsAd))]);
    double spectral_ub = max(abs(i_evalAd).ub(), abs(i_evalAd_min).ub()); // for symmetric matrix it is upper bound of max Eig
    //Apply the theorem
    //ienc = infsup( inf(iEigsAc) - spectral, sup(iEigsAc) + spectral );
    ibex::Interval iLam = ibex::Interval(i_evalAc.lb()-spectral_ub,i_evalAc.ub()+spectral_ub);

    if(_DEBUG_) {
        cout << endl << "rohn w.o. verify:	" << evalAc(dim-1) + max(abs(evalAd(dim-1)), abs(evalAd(0))) << endl << endl;
        cout << endl << "iA: " << endl << Ad << endl;
        cout << "Max Eigenvalue of iA: " << endl;
        cout << "[" << iLam.lb() << ", " << iLam.ub() << "]" << endl;
    }
    return iLam;

}


//Rohn Algorithm for interval matrix
double rohnEigMax(ibex::IntervalMatrix iA, int dim) {
    double lamMax;
    ibex::Matrix dAc = iA.mid(); // center
    ibex::Matrix dAd = iA.rad(); // radius

    MatrixXd Ac = ibex_to_eigen(dAc, dim);
    MatrixXd Ad = ibex_to_eigen(dAd, dim);
    if(_DEBUG_){
        cout << endl << "Inside rohnEig: " << endl;
        cout << "Ac: " << endl << Ac << endl;
        cout << "Ad: " << endl << Ad << endl;
    }
    ibex::Interval iLam = rohnEig(Ac, Ad, dim); // returns ibex::Interval
    lamMax = iLam.ub();
    return lamMax;
}

// Combine All Methods
double effEigMax(IntervalMatrix iA, int dim){
    double best = INF;
    double lam = rohnEigMax(iA, dim);

    if(lam < best)
        best = lam;

    lam = cfEigMax(iA, dim);

    if(lam < best)
        best = lam;

    return best;
}

double hladikDiagMax(ibex::IntervalMatrix iA, int dim){
    IntervalMatrix iB(iA);
    for(int i=0; i < dim; i++){
        double diagMax = iA[i][i].ub();
        iB[i][i] = ibex::Interval(diagMax, diagMax);
    }
    double lam = effEigMax(iB, dim);
    return lam;
}

// CombinedApproach
double computeIntEig(ibex::IntervalMatrix iA, int dim){
    double best = INF;
    best = hladikDiagMax(iA, dim);
    return best;
}

#endif /* COMPUTE_INTEIG_H */
