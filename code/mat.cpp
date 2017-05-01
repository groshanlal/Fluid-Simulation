#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace Eigen;
using namespace std;


int main()
{
/*
  MatrixXd m(2,2);
  m(0,0) = 3;
  m(1,0) = 2.5;
  m(0,1) = -1;
  m(1,1) = m(1,0) + m(0,1);
  std::cout << m << std::endl;
*/
  int n = 2;
  VectorXd x(n), b(n);
  SparseMatrix<double> A(n,n);

  typedef Eigen::Triplet<double> T;
  vector<T> tripletList;
  tripletList.push_back(T(0,0,1));
  tripletList.push_back(T(0,0,10));
  tripletList.push_back(T(0,1,2));
  tripletList.push_back(T(1,0,3));
  tripletList.push_back(T(1,1,50));


  A.setFromTriplets(tripletList.begin(), tripletList.end(), [] (const float&,const float &b) { return b; });
  cout << "A: \n" << A << endl;
  A.setFromTriplets(tripletList.begin(), tripletList.end());
  cout << "A: \n" << A << endl;
  

  b(0) = 4;
  b(1) = 10;
  ConjugateGradient<SparseMatrix<double>, Lower|Upper> cg;
  cg.setMaxIterations(10);
  cg.compute(A);
  x = cg.solve(b);
  cout << "#iterations:     " << cg.iterations() << endl;
  cout << "estimated error: " << cg.error()      << endl;
  cout << "x: \n" << x << endl;
  cout << "b: \n" << b << endl;
}
