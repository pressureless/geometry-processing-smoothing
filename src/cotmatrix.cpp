#include "cotmatrix.h"
#include <iostream>

void cotmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L)
{
  // Add your code here
	int num = F.maxCoeff()+1;
	L.resize(num, num);
	L.setZero();
	std::vector<Eigen::Triplet<double>> tripletList; 
    Eigen::MatrixXd angle = Eigen::MatrixXd(num, num);
    angle.setZero();
    Eigen::MatrixXd tag = Eigen::MatrixXd(num, num);
    tag.setZero();

    for (int i = 0; i < F.rows(); ++i)
    {
    	double a2 = (l(i,0)*l(i,0) + l(i,1)*l(i,1) - l(i,2)*l(i,2)) / (2 * l(i,0) * l(i,1)); //cos
    	a2 = a2/std::sqrt(1 - a2*a2); //cot
    	double a3 = (l(i,1)*l(i,1) + l(i,2)*l(i,2) - l(i,0)*l(i,0)) / (2 * l(i,1) * l(i,2)); //cos
    	a3 = a3/std::sqrt(1 - a3*a3); //cot
    	double a1 = (l(i,0)*l(i,0) + l(i,2)*l(i,2) - l(i,1)*l(i,1)) / (2 * l(i,0) * l(i,2));  //cos
    	a1 = a1/std::sqrt(1 - a1*a1); //cot

    	angle(F(i,0), F(i,1)) += a2/2;
    	angle(F(i,1), F(i,2)) += a3/2;
    	angle(F(i,2), F(i,0)) += a1/2;

    	angle(F(i,1), F(i,0)) = angle(F(i,0), F(i,1));
    	angle(F(i,2), F(i,1)) = angle(F(i,1), F(i,2));
    	angle(F(i,0), F(i,2)) = angle(F(i,2), F(i,0));

    	angle(F(i,0), F(i,0)) -= a2/2;
    	angle(F(i,1), F(i,1)) -= a2/2;

    	angle(F(i,1), F(i,1)) -= a3/2;
    	angle(F(i,2), F(i,2)) -= a3/2;

    	angle(F(i,2), F(i,2)) -= a1/2;
    	angle(F(i,0), F(i,0)) -= a1/2;

    	tag(F(i,0), F(i,1)) = 1;
    	tag(F(i,1), F(i,2)) = 1;
    	tag(F(i,2), F(i,0)) = 1;

    	tag(F(i,1), F(i,0)) = 1;
    	tag(F(i,2), F(i,1)) = 1;
    	tag(F(i,0), F(i,2)) = 1;

    	tag(F(i,0), F(i,0)) = 1;
    	tag(F(i,1), F(i,1)) = 1;
    	tag(F(i,2), F(i,2)) = 1; 
    } 
    //add to SparseMatrix L
    for (int i = 0; i < num; ++i)
    {
    	for (int j = i; j < num; ++j)
    	{
    		if (tag(i,j) > 0)
    		{
	    		tripletList.push_back(Eigen::Triplet<double>(i, j, angle(i,j) ));
	    		if (i != j)
	    		{
	    			tripletList.push_back(Eigen::Triplet<double>(j, i, angle(i,j) ));
	    		} 
    		}
    	}
    } 
	L.setFromTriplets(tripletList.begin(), tripletList.end());
}

