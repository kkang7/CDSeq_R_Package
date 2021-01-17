// [[Rcpp::depends(RcppArmadillo)]]
# include <RcppArmadillo.h>

#include <Rcpp.h>
#include <stdio.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <chrono>
#include <ctime> 
//#include <numeric> // will use std::accumulate function, for sum up the elements in mixtureSamples

#include "cokus.h"

using namespace std;
using namespace Rcpp;

//' This is the Gibbs sampler for CDSeq.
//' \code{GibbsSampler} returns estimated GEPs and cell type proportions.
//' @param ALPHA hyperparameter for cell type proportion.
//' @param BETA hyperparameter for cell-type-specific GEPs.
//' @param mixtureSamples bulk RNA-seq data in form of read counts.
//' @param T number of cell types.
//' @param NN number of MCMC iteration.
//' @param OUTPUT MCMC progress output control.
//' @param processID worker process ID when using parallel computing.
//' @param data_block_idx index for data blocks from bulk RNA-seq input.
//' @param CDSeq_tmp_log temporary log file recording the workers' jobs.
//' @param write_2_file print to progress msg to CDSeq_tmp_log if it is 1, not printing otherwise.
//' @return random integers  uniformly distributed in 0..(2^32 - 1).
// [[Rcpp::export]]
List gibbsSampler(double ALPHA, std::vector<double> BETA, NumericMatrix mixtureSamples, int T, int NN, int OUTPUT, int processID, int data_block_idx, std::string CDSeq_tmp_log, int write_2_file)
{
  //clock_t start,finish; 
  //start=clock();
  /*
   *  initialize variables
   */
  int genes = mixtureSamples.rows(), samples = mixtureSamples.cols(); 
  int n = 0;
  //n = std::accumulate(mixtureSamples.begin(),mixtureSamples.end(),0); 
  vector<unsigned int> sid, gid;// unsigned int for saving space
  for(int i = 0; i < genes; i ++)
  {
    for(int j = 0; j < samples; j ++)
    {
      n += mixtureSamples(i, j);
      for(int k = 0; k < mixtureSamples(i, j); k ++)
      {
        sid.push_back(j);
        gid.push_back(i);
      }
    }
  }
  
  vector<int>  cellTypeTot(T),csGEP_vec(genes*T), SSP_vec(T*samples);
  vector<double> probs(T);
  vector<unsigned int> order(n), cellTypeBag(n);// unsigned int for saving space
  
  // this matrix stores reads-cell-type-assignment for each sample
  // use RcppArmadillo cube data structure as 3D matrix container
  arma::cube cellTypeAssignSplit(genes,samples,T);
  cellTypeAssignSplit.fill(0);// initialize with zeros
  
  int wi,di,i,ii,j,cellType, rp, iter, wioffset, dioffset;
  double totprob, WBETA, r, max;
  ofstream logfile (CDSeq_tmp_log, std::ios_base::app);
  auto time_now = std::chrono::system_clock::now();
  std::time_t cdseq_time = std::chrono::system_clock::to_time_t(time_now);
  // seeding
  int SEED = 3;
  seedMT( 1 + SEED * 2 ); // seeding only works on uneven numbers
  for (i = 0; i < n; i ++)
    {
      //cout << i << "/" << n << endl;
      wi = gid[i];
      di = sid[i];
      // pick a random cell type 0..T-1
       cellType = (int) ( (double) randomMT() * (double) T / (double) (4294967296.0 + 1.0) ); // To prevent overflow
      cellTypeBag[i] = cellType; // assign this gene read to this cell type
      csGEP_vec[wi*T + cellType]++;
      SSP_vec[di*T + cellType]++;
      cellTypeTot[cellType]++; // increment cellTypeTot matrix
      cellTypeAssignSplit(wi,di,cellType)++;
    }
  for (i = 0; i < n; i ++) order[i]=i; // fill with increasing series
  
  for (i = 0; i < (n-1); i ++) {
    // pick a random integer between i and nw
    rp = i + (int) ((double) (n-i) * (double) randomMT() / (double) (4294967296.0 + 1.0));
    
    // switch contents on position i and position rp
    swap(order[rp], order[i]);
  }
  
  WBETA = (double) (accumulate(BETA.begin(), BETA.end(), 0));
  for (iter = 0; iter <= NN; iter ++) 
  {
    Rcpp::checkUserInterrupt();
    if (OUTPUT >= 1) {
      if ((iter % 10) == 0) Rprintf( "%-4d of %-4d MCMC runs (finished %2d%%, press esc at any time to terminate the MCMC iteration)\r" , iter , NN , 100*iter/NN);
      if (iter == NN) Rprintf( "%-4d of %-4d MCMC runs (finished %2d%%)\nGibbs sampler completed succefully\n" , iter , NN, 100*iter/NN );
    }
    if(OUTPUT==0 && write_2_file==1){
      //if (iter % 10 == 0){cout<<"Running in parallel:"<< 100*iter/NN<<"% of the job is finished\r";}
      //if (iter % 10 == 0){Rcpp::Rcout<<"Cell type number is"<<setw(4)<<T<<", data block ["<<setw(4)<<data_block_idx<<"] finished "<<setw(4)<<100*iter/NN<<"% of MCMC iterations (worker process ID "<<setw(7)<<processID<<")\r";}
      //if(iter == NN){Rcpp::Rcout<<"Cell type number is"<<setw(4)<<T<<", data block ["<<setw(4)<<data_block_idx<<"] finished "<<setw(4)<<100*iter/NN<<"% of MCMC iterations (worker process ID "<<setw(7)<<processID<<")"<<std::endl;}
      if (logfile.is_open()){
        if (iter == 0){logfile<<std::ctime(&cdseq_time);}
        if (iter % 10 == 0){logfile<<"Cell type number is"<<setw(4)<<T<<", data block ["<<setw(4)<<data_block_idx<<"] finished "<<setw(4)<<100*iter/NN<<"% of MCMC iterations (worker process ID "<<setw(7)<<processID<<")\n";}
        if (iter == NN){logfile << "\n\n";}
        //if(iter == NN){logfile<<"Cell type number is"<<setw(4)<<T<<", data block ["<<setw(4)<<data_block_idx<<"] finished "<<setw(4)<<100*iter/NN<<"% of MCMC iterations (worker process ID "<<setw(7)<<processID<<")\n";}
      }  
    }
    
    for (ii = 0; ii < n; ii++) {
      i = order[ii]; // current gene read to assess
      wi  = gid[i]; // current gene index
      di  = sid[i]; // current sample index  
      cellType = cellTypeBag[i]; // current cell type assignment to gene read
      cellTypeTot[cellType]--;  // substract this from counts
      cellTypeAssignSplit(wi,di,cellType)--;
      
      wioffset = wi*T;
      dioffset = di*T;
      
      csGEP_vec[wioffset+cellType]--;
      SSP_vec[dioffset+cellType]--;
      
      totprob = (double) 0;
      for (j = 0; j < T; j++) {
        probs[j] = ((double) csGEP_vec[wioffset + j] + (double) BETA[wi])/( (double) cellTypeTot[j]+ (double) WBETA)*( (double) SSP_vec[dioffset+ j] + (double) ALPHA);
        totprob += probs[j];
      }
      // sample a cell type from the distribution
      r = (double) totprob * (double) randomMT() / (double) 4294967296.0;
      max = probs[0];
      cellType = 0;
      while (r > max) {
        cellType ++;
        max += probs[cellType];
      }
      cellTypeBag[i] = cellType; // assign current gene read i to cell type j
      csGEP_vec[wioffset + cellType ]++; // and update counts
      SSP_vec[dioffset + cellType ]++;
      cellTypeTot[cellType] ++;
      cellTypeAssignSplit(wi,di,cellType)++;//keep read-cell-type-assignment for each cell type in each sample
    }
  }
  logfile.close();
  List result;
  result["csGEP_vec"] = csGEP_vec;
  result["SSP_vec"] = SSP_vec;
  result["cellTypeTot"] = cellTypeTot;
  result["cellTypeAssignSplit"] = cellTypeAssignSplit;
  return result;
}



