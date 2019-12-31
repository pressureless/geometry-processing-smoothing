#include "massmatrix.h"
#include <cmath>
#include <iostream>

void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{
  // Add your code here
	int num = F.maxCoeff()+1;
	Eigen::VectorXd areas(num);
    areas.setZero();
	for (int i = 0; i < F.rows(); ++i)
    {
    	double s = (l(i,0) + l(i,1) + l(i,2))/2;
    	double area = std::sqrt(s * (s-l(i,0)) * (s-l(i,1)) * (s-l(i,2)));
    	areas(F(i,0)) += area;
    	areas(F(i,1)) += area;
    	areas(F(i,2)) += area; 
    }
    areas = areas / 3;
    M = areas.asDiagonal();
}

