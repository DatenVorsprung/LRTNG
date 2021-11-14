/*
 * ibex_matrixInv.h
 *
 ** Lagrangian Reachtubes: The next Generation
 *
 *      Authors: Sophie Gruenbacher and Md Ariful Islam
 *      Contact: sophie.gruenbacher@tuwien.ac.at
 */

#ifndef IBEX_MATRIX_INV_H
#define IBEX_MATRIX_INV_H

#include <vector>
#include <stdexcept>
#include <iostream>
#include <math.h>
#include "ibex.h"

using namespace ibex;
using namespace std;

// division operation for vector: result[i] = A[i]/b
ibex::Vector divVector(ibex::Vector& A, double b){
	int dim = A.size();
	ibex::Vector result(dim);
	for(int i=0; i < dim; i++){
		result[i] = A[i]/b;
	}
	return result;
}

// division operation for interval vector: result[i] = A[i]/b
ibex::IntervalVector divIVector(ibex::IntervalVector& A, ibex::Interval b){
  int dim = A.size();
  ibex::IntervalVector result(dim);
  for(int i=0; i < dim; i++){
    result[i] = A[i]/b;
  }
  return result;
}

// need for interval
bool is_gt(ibex::Interval ii1, ibex::Interval ii2){
  return (ii1.lb()> ii2.ub());
}

// for point 
bool isSingular(double A) {
    return (A == 0);
}

// for interval  
bool isSingular(const ibex::Interval& A_x) {
    return ((A_x.lb()<=0) && (A_x.ub()>=0));
}


/********************************************************************************/
/***************** Matrix Inverse for Point Matrix ***************************/
/********************************************************************************/

void gauss(ibex::Matrix& a, ibex::Matrix& b, ibex::Matrix& result)
{
  //typedef typename MatrixType::RowVectorType VectorType;
  //typedef typename MatrixType::ScalarType ScalarType;
  // Scalartype is double
  int i,j,k;
  int dimension = a.nb_rows();
  int *p = new int[dimension];
  for(i=0;i<dimension;++i) 
    p[i]=i;
	
  for(j=0;j<dimension-1;++j)
  {
  	/* Find the maximum magnitude in j-th col */
    double a_max,a_max_temp; // need to change for interval
    int i_max,temp;
    a_max = std::abs(a[p[j]][j]); // need to change for interval: ibex::abs()
    i_max = j;
    temp = p[j];
    for(i=j+1;i<dimension;++i)
    {
      a_max_temp = std::abs(a[p[i]][j]); // need to change for interval: ibex::abs()
      if(a_max_temp > a_max)
      {
        a_max = a_max_temp;
        i_max = i;
      }
    }
    /* swapping the index between (j,j) entry with (jmax, j)*/ 
    p[j] = p[i_max]; // p[j]: index of maximum element
    p[i_max] = temp; // p[i_max]: index of diagonal element

    double divisor = a[p[j]][j]; // maximum magnitude in i-th row
    //cout << divisor << endl;
    if( isSingular(divisor))
    {
      throw std::runtime_error( "Gauss elimination: singular matrix\n");
    }

    // 
    for(i=j+1;i<dimension;++i)
    {
      double factor = a[p[i]][j]/a[p[j]][j];
      a[p[i]][j] = 0;
      b.row(p[i]) -= factor*b.row(p[j]);
      for(k=j+1;k<dimension;++k)
      {
        a[p[i]][k] -= factor*a[p[j]][k];
      }
    }
    
  } // End of outer loop

  if( isSingular(a[p[dimension-1]][dimension-1]) )
  {
    throw std::runtime_error( "Gauss elimination: singular matrix\n");
  }

  for(int i=dimension-1;i>=0;--i)
  {
    //result[i] = b[p[i]];
    result.row(i) = b.row(p[i]);
    for(int j=i+1;j<dimension;++j)
      result.row(i) -= a[p[i]][j]*result.row(j);
    result.row(i) = divVector(result.row(i),a[p[i]][i]);
  }
  delete[] p;
}

// for point matrix
ibex::Matrix gaussInverseMatrix(ibex::Matrix a)
{
   //int rows = A.nb_rows();
   //int cols = A.nb_cols();	
  if(a.nb_rows()!=a.nb_cols())
  {
    throw std::runtime_error("Cannot inverse nonsquare matrix!");
  }
  int dim = a.nb_rows();
  ibex::Matrix result(dim, dim);
  ibex::Matrix b = ibex::Matrix::eye(dim);

  gauss(a,b,result);

  return result;
}

/********************************************************************************/
/***************** Matrix Inverse for Interval Matrix ***************************/
/********************************************************************************/
void gauss(ibex::IntervalMatrix& a, ibex::IntervalMatrix& b, ibex::IntervalMatrix& result)
{
  int i,j,k;
  int dimension = a.nb_rows();
  int *p = new int[dimension];
  for(i=0;i<dimension;++i) 
    p[i]=i;
  
  for(j=0;j<dimension-1;++j)
  {
    /* Find the maximum magnitude in j-th col */
    ibex::Interval a_max,a_max_temp; // need to change for interval
    int i_max,temp;
    a_max = ibex::abs(a[p[j]][j]); // need to change for interval: ibex::abs()
    i_max = j;
    temp = p[j];
    for(i=j+1;i<dimension;++i)
    {
      a_max_temp = ibex::abs(a[p[i]][j]); // need to change for interval: ibex::abs()
      if(is_gt(a_max_temp, a_max)) // a_max_temp > a_max
      {
        a_max = a_max_temp;
        i_max = i;
      }
    }
    /* swapping the index between (j,j) entry with (jmax, j)*/ 
    p[j] = p[i_max]; // p[j]: index of maximum element
    p[i_max] = temp; // p[i_max]: index of diagonal element

    ibex::Interval divisor = a[p[j]][j]; // maximum magnitude in i-th row
    //cout << divisor << endl;
    if( isSingular(divisor))
    {
      throw std::runtime_error( "Gauss elimination: singular matrix\n");
    }

    // 
    for(i=j+1;i<dimension;++i)
    {
      ibex::Interval factor = a[p[i]][j]/a[p[j]][j];
      a[p[i]][j] = 0;
      b.row(p[i]) -= factor*b.row(p[j]);
      for(k=j+1;k<dimension;++k)
      {
        a[p[i]][k] -= factor*a[p[j]][k];
      }
    }
    
  } // End of outer loop

  if( isSingular(a[p[dimension-1]][dimension-1]) )
  {
    throw std::runtime_error( "Gauss elimination: singular matrix\n");
  }

  for(int i=dimension-1;i>=0;--i)
  {
    //result[i] = b[p[i]];
    result.row(i) = b.row(p[i]);
    for(int j=i+1;j<dimension;++j)
      result.row(i) -= a[p[i]][j]*result.row(j);
    result.row(i) = divIVector(result.row(i),a[p[i]][i]);
  }
  delete[] p;
}

//Inverse of Interval Matrix
ibex::IntervalMatrix gaussInverseMatrix(ibex::IntervalMatrix a)
{
   //int rows = A.nb_rows();
   //int cols = A.nb_cols();  
  if(a.nb_rows()!=a.nb_cols())
  {
    throw std::runtime_error("Cannot inverse nonsquare matrix!");
  }
  int dim = a.nb_rows();
  ibex::IntervalMatrix result(dim, dim);
  ibex::IntervalMatrix b = ibex::Matrix::eye(dim);

  gauss(a,b,result);

  return result;
}

#endif /* IBEX_MATRIX_INV_H */
