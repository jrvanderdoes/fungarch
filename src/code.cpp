#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
double compute_sigma_inp_k_q_cpp(unsigned int k, arma::mat d, arma::mat As,
                                 arma::mat prev_ys, unsigned int q,int M){

  double sigma;
  arma::mat A, tmp;

  sigma=d(k,0);

  for(unsigned int i=0;i<q;i++){
    A=As.submat(k+i*M,0,k+i*M,M-1);
    tmp=(A*prev_ys.submat(0,i,M-1,i));
    sigma=sigma+tmp(0,0);
  }

  return sigma;
}

// [[Rcpp::export]]
double obj_fn_q_cpp(arma::mat d, arma::mat As,
                    arma::mat ysq,
                    unsigned int q){

  //initialize
  arma::mat sigma0,y0,prev_sigma,prev_ys,new_sigma, sigmat;
  unsigned int cap=q;
  double obj;

  //ysq=NxM
  unsigned int N=size(ysq)(0);
  unsigned int M=size(ysq)(1);


  sigma0=sum(ysq.submat(0,0,4,M-1),0)/5;
  y0=sigma0;


  //compute initial values to start recursion

  prev_ys.set_size(q,M);

  for(unsigned int i=0;i<q;i++)
    prev_ys.row(i)=sigma0;

  //get initial q sigma inner products:

  new_sigma.set_size(1,M);


  cap=q;

  sigmat.set_size(cap,M);

  for(unsigned int i=0;i<cap;i++){


    for(unsigned int j=0;j<M;j++){
      new_sigma(0,j)=compute_sigma_inp_k_q_cpp(j,d,As,prev_ys.t(),q,M);
    }
    sigmat.row(i)=new_sigma;

    if(q>1)
      prev_ys.submat(0,0,q-2,M-1)=prev_ys.submat(1,0,q-1,M-1);
    prev_ys.row(q-1)=ysq.row(i);

  }


  for(unsigned int i=cap;i<N;i++){
    //just do this matrix type in case q=1
    prev_ys=ysq.submat(i-q,0,(i-1),M-1);
    //recursively compute sigma inner products

    for(unsigned int j=0;j<M;j++){
      new_sigma(0,j)=compute_sigma_inp_k_q_cpp(j,d,As,prev_ys.t(),q,M);
    }

    sigmat=join_vert(sigmat,new_sigma);
  }


  obj=mean(sum(ysq/sigmat+log(sigmat),1));


  return obj;
}
