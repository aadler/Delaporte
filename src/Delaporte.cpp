#include <Rcpp.h>
#include <vector>
#include <algorithm>

using namespace Rcpp;

double ddelap_C_S(double x, double alpha, double beta, double lambda, bool lg) {
	int k = ceil(x);
	double del_pmf_s = 0.0;
	for(int i = 0; i <= k; ++i) { //using logs to prevent under/overflow
		del_pmf_s += exp(lgamma(alpha + i) + i * log(beta) + (k - i) * log(lambda) -
              lambda - lgamma(alpha) - lgamma(i + 1) - (alpha + i) * log1p(beta) -
              lgamma(k - i + 1));  
	}
  if (lg == TRUE) del_pmf_s = log(del_pmf_s);
	return (del_pmf_s);
}

// [[Rcpp::export]]
NumericVector ddelap_C(NumericVector x, NumericVector alpha, NumericVector beta,
                       NumericVector lambda, bool lg) {
  int n = x.size();
  int a_size = alpha.size();
  int b_size = beta.size();
  int l_size = lambda.size();
  NumericVector del_pmf(n);
  for (int t = 0; t < n; ++t) {
    del_pmf[t] = ddelap_C_S(x[t], alpha[t % a_size], beta[t % b_size], lambda[t % l_size], lg);
  }
  return(del_pmf);
}

// [[Rcpp::export]]
std::vector <double> pdelap_C(std::vector <double> q, std::vector <double> alpha, std::vector <double> beta,
                              std::vector <double> lambda, bool lt, bool lp) {
  std::vector <double>::size_type n = q.size();
  std::vector <double>::size_type a_size = alpha.size();
  std::vector <double>::size_type b_size = beta.size();
  std::vector <double>::size_type l_size = lambda.size();
  std::vector <double> del_cdf(n);
  if (1 == a_size && a_size == b_size && b_size == l_size) {
    /* if parameters are all singletons (not vectors) then the idea is to find the largest value in the vector
     * and build the PDF up to that point. Every other value will be a lookup off of the largest vector.
     * Otherwise each entry will need to build its own value.
     */ 
	  std::vector <double>::iterator MX = std::max_element(q.begin(), q.end());
	  int top = ceil(*MX);
	  std::vector <double> single_cdf_vector(top + 1);
	  single_cdf_vector[0] = exp(-lambda[0]) / pow((1 + beta[0]), alpha[0]);
	  for (int t = 1; t <= top; ++t) {
		  double cdf_at_t = 0.0;
		  for (int i = 0; i <= t; ++i) {
			  cdf_at_t += exp(lgamma(alpha[0] + i) + i * log(beta[0]) + (t - i) * log(lambda[0]) -
                        lambda[0] - lgamma(alpha[0]) - lgamma(i+1) -  (alpha[0] + i) * log1p(beta[0]) -
                        lgamma(t - i + 1));			
			}
		  single_cdf_vector[t] = single_cdf_vector[t - 1] + cdf_at_t;
	  }
	  for (std::vector <double>::size_type t = 0; t < n; ++t) {
		  del_cdf[t] = single_cdf_vector[ceil(q[t])];
	  }
	  if (lt==FALSE) {
		  for (std::vector <double>::size_type t = 0; t < n; ++t) {
			  del_cdf[t] = 1.0 - del_cdf[t];
		  }
	  }
	  if (lp==TRUE) {
		  for (std::vector <double>::size_type t = 0; t < n; ++t) {
			  del_cdf[t] = log(del_cdf[t]);
		  }
	  }
  } else {
    //Here will go code to build each PDF entry individually
  }
	return (del_cdf);
}

bool OUTSIDE01 (double quantile) {return (quantile <= 0.0 || quantile >= 1.0);}
//above needed to trim vector of values <=0 and >= 1

// [[Rcpp::export]]
std::vector <double> qdelap_C(std::vector <double> p, double alpha, double beta, double lambda, bool lt, bool lp) {
// The idea is to find the largest CDF point in the vector, and build counts up to that point.
// Every other value is a lookup off of the largest vector.
	std::vector <double>::size_type n = p.size();
  std::vector <double> pcopy = p;
  if (lp==TRUE) {
    for (std::vector <double>::size_type t = 0; t < n; t++) {
			pcopy[t] = exp(pcopy[t]);
		}
	}
  if (lt==FALSE) {
		for (std::vector <double>::size_type t = 0; t < n; t++) {
			pcopy[t] = 1.0 - pcopy[t];
		}
	}
  std::vector <double>::iterator pcbegin = pcopy.begin();
  std::vector <double>::iterator pcend = pcopy.end();
  pcend = std::remove_if(pcbegin, pcend, OUTSIDE01);
  std::vector <double> pcopy2;
  pcopy2.assign (pcbegin, pcend);
  std::vector <double>::iterator MX = std::max_element(pcopy2.begin(), pcopy2.end());
  pcopy.clear();
	double maxquantile = *MX;
  double cdftop = 0.0;
  std::vector <double> CDFVEC;
  CDFVEC.push_back(exp(-lambda)/pow((1+beta), alpha)); //pre-load 0 value
  int cap = 1; //Will be "max integer"
  while (cdftop < maxquantile) {
    double PMF = 0.0;
  	for(int i = 0; i <= cap; i++) {
			PMF += exp(lgamma(alpha+i)+i*log(beta)+(cap-i)*log(lambda)-lambda-lgamma(alpha)-lgamma(i+1)-(alpha+i)*log(1+beta)-lgamma(cap-i+1));			
		}
    PMF = PMF +  CDFVEC[cap-1]; //add pmf value at cap to previous CDF value
    CDFVEC.push_back(PMF);
    cdftop = CDFVEC[cap];
    ++cap;
  }
  std::vector <double> RETVEC(n);
  std::vector<double>::iterator foundit;
  for (std::vector <double>::size_type t = 0; t < n; t++) {
      if (p[t]< 0) {
        RETVEC[t] = std::numeric_limits<double>::quiet_NaN();
      } else if (p[t]==0) {
        RETVEC[t] = 0;
      } else if (p[t] >= 1) {
        RETVEC[t] = std::numeric_limits<double>::infinity();
      } else {
        foundit = std::upper_bound (CDFVEC.begin(), CDFVEC.end(), p[t]);
        double spot = foundit - CDFVEC.begin();
        RETVEC[t] = spot;
		  }
  }
	return (RETVEC);
}

// [[Rcpp::export]]
std::vector <double> rdelap_C(int p, double alpha, double beta, double lambda) {
// The idea is to find the largest CDF point in the vector, and build counts up to that point.
// Every other value is a lookup off of the largest vector.
  RNGScope scope;
  NumericVector RUNI = runif(p, 0.0, 1.0);
  NumericVector::iterator MX = std::max_element(RUNI.begin(), RUNI.end());
	double maxquantile = *MX;
  double cdftop = 0.0;
  std::vector <double> CDFVEC;
  CDFVEC.push_back(exp(-lambda)/pow((1+beta), alpha)); //pre-load 0 value
  int cap = 1; //Will be "max integer"
  while (cdftop < maxquantile) {
    double CCC = 0.0;
  	for(int i = 0; i <= cap; i++) {
			CCC += exp(lgamma(alpha+i)+i*log(beta)+(cap-i)*log(lambda)-lambda-lgamma(alpha)-lgamma(i+1)-(alpha+i)*log(1+beta)-lgamma(cap-i+1));			
		}
    CCC = CCC +  CDFVEC[cap-1];
    CDFVEC.push_back(CCC);
    cdftop = CDFVEC[cap];
    ++cap;
  }
  std::vector <double> RETVEC(p);
  std::vector <double>::iterator foundit;
  for (int t = 0; t < p; t++) {
      if (RUNI[t] == 0) {
        RETVEC[t] = 0;
      } else if (RUNI[t] == 1) {
        RETVEC[t] = std::numeric_limits<double>::infinity();
      } else {
        foundit = std::upper_bound (CDFVEC.begin(), CDFVEC.end(), RUNI[t]);
        double spot = foundit - CDFVEC.begin();
        RETVEC[t] = spot;
		  }
  }
	return (RETVEC);
}

// [[Rcpp::export]]
NumericVector MoMdelap_C(NumericVector X){
// Using the definitions for the mean, variance and skew of the Delaporte, find the method of moments
// parameter estimates for a vector of data. This is also good starting point for maximum likelihood.
  int n = X.size();
  double nm1 = n - 1.0;
  double P = n * sqrt(nm1) / (n - 2.0);
  double Mu_D = 0;
  double M2 = 0;
	double M3 = 0;
  for (int i = 0; i < n; i++) {
	   double delta = X(i) - Mu_D;
     double delta_i = delta / (i + 1);
      double T1 = delta * delta_i * i;
      Mu_D += delta_i;
      M3 += (T1 * delta_i * (i - 1) - 3 * delta_i * M2);
      M2 += T1;
    }
  double Var_D = M2 / nm1;
  double Skew_D = P * M3 / pow(M2, 1.5);
  double VmM_D = Var_D - Mu_D;
  double beta = 0.5 * (Skew_D * pow(Var_D, 1.5) / VmM_D - 3);
  double alpha = VmM_D / (beta * beta);
  double lambda = Mu_D - alpha * beta;
  return(NumericVector::create(alpha, beta, lambda));
}