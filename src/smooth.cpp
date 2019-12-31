#include "smooth.h"
#include "cotmatrix.h"
#include "massmatrix.h"
#include <igl/edge_lengths.h>
#include <igl/cotmatrix.h> 
#include <igl/massmatrix.h>
#include <Eigen/Core>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseCore>
#include <iostream>
#include <Eigen/Dense>


void smooth(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & G,
    double lambda,
    Eigen::MatrixXd & U)
{
  // Replace with your code 
	Eigen::MatrixXd l;
	igl::edge_lengths(V,F,l); 

	Eigen::SparseMatrix<double> L;
	cotmatrix(l, F, L); //the order of l is not easy to understand, watch out

    Eigen::DiagonalMatrix<double, Eigen::Dynamic> M;
    massmatrix(l, F, M); 

    Eigen::SparseMatrix<double> MM(M.rows(), M.cols());
    MM.setZero();

	std::vector<Eigen::Triplet<double>> tripletList;
    for (int i = 0; i < M.rows(); ++i)
    {
    	tripletList.push_back(Eigen::Triplet<double>(i, i, M.diagonal()(i) ));
    } 
	MM.setFromTriplets(tripletList.begin(), tripletList.end());

    Eigen::SparseMatrix<double> A = MM-lambda*L;
    // Eigen::SparseMatrix<double> A = MM - 0.001*L;

    // Eigen::SimplicialLLT <Eigen::SparseMatrix<double> > cholesky;
    // cholesky.analyzePattern(A);
    // cholesky.factorize(A);
    // Eigen::MatrixXd matrixL = cholesky.matrixL();
    // Eigen::MatrixXd y = matrixL.colPivHouseholderQr().solve( M*G );
    // U = matrixL.transpose().colPivHouseholderQr().solve(y);
 
	Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >   solver;
	solver.analyzePattern(A); 
	solver.factorize(A);  
	U = solver.solve(MM*G);   
}
