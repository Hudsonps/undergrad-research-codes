#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include "erro.h"
#include "lanczos.h"
#include "indata.h"

double fdot(double* u, double* v)
// returns dot product of vectors u and v, of dimension dim;
// to improve accuracy, run outwards from Fermi level 
{
  double ret = 0.0;
  int ellfermi = NBAND / 2;
  for(int ell = 0; ell < ellfermi; ell++){
    ret += u[ellfermi - ell - 1] * v[ellfermi - ell -1]; // below Fermi level
    ret += u[ellfermi + ell] * v[ellfermi + ell];	 // above
  }
  return ret;
}

double gram_schmidt_ortho(double** alpha, int nvect)
// Gram_Schmidt orthonormalize vector alpha[nvect] to preceding vectors
// alpha[n] (n=0,..,nvect-1];
// Each vector has NBAND components alpha[n][ell] (ell=0,..,NBAND-1);
// Norm of resulting vector, prior to normalization, is returned.
{
  for(int n = 0; n < nvect; n++){
    double proj = fdot(alpha[nvect], alpha[n]) / fdot(alpha[n], alpha[n]);
    for(int ell = 0; ell < NBAND; ell++)
      alpha[nvect][ell] -= proj * alpha[n][ell];
  }
  double norm = 0.0;
  int fermi = NBAND / 2;

  for(int ell = 0; ell < fermi; ell++){
    norm += pow(alpha[nvect][fermi - ell -1], 2.0);
    norm += pow(alpha[nvect][fermi + ell], 2.0);
  }

  for(int ell = 0; ell < NBAND; ell++)
    alpha[nvect][ell] /= sqrt(norm);
  
  return sqrt(norm);
}

double codiags(double *erg, double** alpha, double* t, int dim)
// computes codiagonal matrix elements of tridiagonal Hamiltonian
// H=\sum_n t_n f_n^\dagger f_{n+1} + Hc,
// equivalent to H=\sum_\ell erg[\ell]c_\ell^\dagger c_\ell; 
///////////////////////////////////////////////////////////////////
// upon calling:
// the vector erg contains the energies erg[\ell] (\ell=0,..,dim-1);
// the vector alpha[0][\ell] (\ell=0,..,dim-1) contains the anticommutators
// {f_0^\dagger, c_{\ell}
//////////////////////////////////////////////////////////////////
// upon returning:
// alpha contains the elements defined by
// f_n=\sum_n\alpha[n][\ell]c_\ell, including the \alpha[0][\ell];
// the vector t contains the codiagonal elements t[n]
//////////////////////////////////////////////////////////////////
{
  double norm = 0.0;
  int n_max = MAX_ITER;
  for(int ell = 0; ell < NBAND; ell++)
    alpha[1][ell] = erg[ell] * alpha[0][ell];
  
  t[0] = gram_schmidt_ortho(alpha, 1);
  
  int asymptotic = 0;
  double decay = 0.0;
  for(int n = 2; n < n_max; n++){
    switch(asymptotic){
    case 0:
      if(n > 3 && (fabs(t[n-2] / t[n-3] - t[n-3] / t[n-4]) < THRESHOLD)){
	asymptotic = n - 2;
	decay = t[n-2] / t[n-3];
      }
      
      for(int ell = 0; ell < NBAND; ell++)
	alpha[n][ell] = erg[ell] * alpha[n-1][ell] - t[n-2] * alpha[n-2][ell];
      t[n-1] = gram_schmidt_ortho(alpha, n);
      for(int prior_n = 0; prior_n < n; prior_n++){
	double error = 0.0;
	if((error = fabs(fdot(alpha[n], alpha[prior_n]))) > TOLERANCE){
	  std::cerr << "\n***Nonrthogonal vectors in codiags: (" 
		    << n << "|" << prior_n << ") = " 
		    << error << std::endl;
	  exit(1);
	}
      }
      break;

    default:
      t[n-1] = t[n-2] * decay;
      break;
    }
  }
}

void ksis(double lambda, double z, double* ksi, int dim)
// computes first dim ksis and store them in vector ksis
{
  std::ofstream fksis("ksis.s");

  double* cband = new double[dim];
  double** alpha = new double*[dim];
  double* t = new double[dim];

  for(int n = 0; n < dim; n++){
    alpha[n] = new double[dim];
  }
  alpha[0][0] = alpha[0][NBAND -1] 
    = sqrt((1.0 - pow(lambda, -z))/2.0);
  cband[0] = -(1.0 + pow(lambda, -z))/(1.0 + pow(lambda, -1.0));
    cband[NBAND - 1] = -cband[0];
  
  for(int ell = 1; ell < NBAND / 2; ell++){
    alpha[0][ell] = alpha[0][NBAND - ell - 1] 
      = sqrt((1.0 - pow(lambda, -1.0))/2.0) * pow(lambda, double(1.0-z-ell) / 2.0);
    cband[ell] = -pow(lambda, double(1.0-z-ell));
    cband[NBAND - ell - 1] = -cband[ell];
  }

  codiags(cband, alpha, t, dim);
  for(int n = 0; n < dim; n++)
    ksi[n] = t[n] * pow(lambda, n / 2.0);
  // for comparison, compute Wilson's coefficients (z==1)
  for(int n = 0; n < dim; n++){
    double lambda_p_npp = pow(lambda, n + 1.);
    double ksi_wils_n = lambda_p_npp * (lambda_p_npp - 1.)
	/ sqrt(lambda_p_npp * lambda_p_npp - lambda)
	/ sqrt(lambda_p_npp * lambda_p_npp - 1. / lambda);
    
    fksis << n << "\t" << std::setw(10) << std::setprecision(13) << ksi[n] 
	  << "\t" << ksi_wils_n << std::endl;
  }
}

double ksi(rundata* rundat, int n, int dim)
// returns codiagonal element ksi_n;
// first time called, computes ksi's and stores them;
// if not given, dim takes default value NBAND
{
  if(n >= dim)
    erro("***n>=dim calling ksi\n");
  static double *_ksis = 0;
  if(!_ksis){
    _ksis = new double [dim];
    ksis(rundat->lambda(), rundat->z(), _ksis, dim);
  }
  return _ksis[n];
}
