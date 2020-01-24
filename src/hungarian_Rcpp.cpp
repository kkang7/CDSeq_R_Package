#include <Rcpp.h>
#include "Hungarian.h"

using namespace std;
using namespace Rcpp;
//' This is the Hungarian algorithm wrapper for cell type assignment
//' \code{hungarian_Rcpp} returns cell type assignment given reference GEPs
//' @param costMat correlation matrix
//' @return cost for the assignment and cell type assignment
// [[Rcpp::export]]
List hungarian_Rcpp(NumericMatrix costMat) {
  int nrows=costMat.rows(), ncols=costMat.cols();
  int i,j;

  vector< vector<double> > costMatrix;
  // the following step is critical to make it work
  costMatrix.resize(nrows, vector<double>(ncols, 0));

  for(i=0;i<nrows;i++)
  {
    for(j=0;j<ncols;j++)
    {
      costMatrix[i][j] = costMat(i,j);
      //Rprintf("costMat(%d,%d)=%f\n",i,j,costMatrix[i][j]);
    }
  }

  HungarianAlgorithm HungAlgo;
  vector<int> assignment;
  
  double cost = HungAlgo.Solve(costMatrix, assignment);
  
  //for (int x = 0; x < costMatrix.size(); x++)
  //  Rprintf("%d,%f\t",x,assignment[x]);
    //std::cout << x << "," << assignment[x] << "\t";
  
  //std::cout << "\ncost: " << cost << std::endl;
  //Rprintf("\ncost:%f\n",cost);
  
  List result;
  result["cost"] = cost;
  result["cost_assignment"] = assignment;
  return result;
}

