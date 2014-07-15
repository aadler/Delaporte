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

double pdelap_C_S(double q, double alpha, double beta, double lambda) {
  int k = ceil(q);
  double del_cdf_s = exp(-lambda) / pow((1 + beta), alpha);
  for (int i = 1; i <= k; ++i) {
    del_cdf_s += ddelap_C_S(i, alpha, beta, lambda, FALSE);
  }
  return(del_cdf_s);
}

// [[Rcpp::export]]
NumericVector pdelap_C(NumericVector q, NumericVector alpha, NumericVector beta,
                              NumericVector lambda, bool lt, bool lp) {
  int n = q.size();
  int a_size = alpha.size();
  int b_size = beta.size();
  int l_size = lambda.size();
  NumericVector del_cdf(n);
  if (1 == a_size && a_size == b_size && b_size == l_size) {
    /* if parameters are all singletons (not vectors) then the idea is to find the largest value in the vector
     * and build the PDF up to that point. Every other value will be a lookup off of the largest vector.
     * Otherwise each entry will need to build its own value.
     */ 
	  NumericVector::iterator MX = std::max_element(q.begin(), q.end());
	  int top = ceil(*MX);
	  NumericVector single_cdf_vector(top + 1);
	  single_cdf_vector[0] = exp(-lambda[0]) / pow((1 + beta[0]), alpha[0]);
	  for (int i = 1; i <= top; ++i) {
		  single_cdf_vector[i] = single_cdf_vector[i - 1] + ddelap_C_S(i, alpha[0], beta[0], lambda[0], FALSE);
	  }
	  for (int i = 0; i < n; ++i) {
		  del_cdf[i] = single_cdf_vector[ceil(q[i])];
	  }
  } else { //Have to build full double summation chain for each entry as the parameters change. Much slower
    for (int i = 0; i < n; ++i) {
    del_cdf[i] = pdelap_C_S(q[i], alpha[i % a_size], beta[i % b_size], lambda[i % l_size]);
    }
  }
	if (lt == FALSE) {
	  for (int i = 0; i < n; ++i) {
		  del_cdf[i] = 1.0 - del_cdf[i];
		}
	}
	if (lp == TRUE) {
	  for (int i = 0; i < n; ++i) {
		  del_cdf[i] = log(del_cdf[i]);
		}
	}
	return (del_cdf);
}

// [[Rcpp::export]]
NumericVector qdelap_C(NumericVector p, NumericVector alpha, NumericVector beta,
                              NumericVector lambda, bool lt, bool lp) {
	int n = p.size();
  int a_size = alpha.size();
  int b_size = beta.size();
  int l_size = lambda.size();
  NumericVector adjusted_p = p;
  NumericVector RETVEC(n);
  if (lp == TRUE) {
    for (int i = 0; i < n; ++i) {
  	  adjusted_p[i] = exp(adjusted_p[i]);
		}
  }
  if (lt == FALSE) {
	  for (int i = 0; i < n; ++i) {
		  adjusted_p[i] = 1.0 - adjusted_p[i];
		}
	}
  for (int i = 0; i < n; ++i) {
    if (adjusted_p[i] < 0) {
      RETVEC[i] = std::numeric_limits<double>::quiet_NaN();
    } else if (adjusted_p[i] == 0) {
      RETVEC[i] = 0;
    } else if (adjusted_p[i] >= 1) {
      RETVEC[i] = std::numeric_limits<double>::infinity();
    } else {
      NumericVector CDFVEC;
      int del_quantile = 0; //Will become "needed integer"
      CDFVEC.push_back(exp(-lambda[i % l_size]) / pow((1 + beta[i % b_size]), alpha[i % a_size])); //pre-load 0 value
      double cdf_curr_top = CDFVEC[0];
      while (cdf_curr_top < adjusted_p[i]) {
        ++del_quantile;
        CDFVEC.push_back(ddelap_C_S(del_quantile, alpha[i % a_size], beta[i % b_size], lambda[i % l_size], FALSE) +  CDFVEC[del_quantile - 1]);
        cdf_curr_top = CDFVEC[del_quantile];
      }
      RETVEC[i] = del_quantile;
  	}
  }
	return (RETVEC);
}

// [[Rcpp::export]]
NumericVector rdelap_C(int n, NumericVector alpha, NumericVector beta,
                              NumericVector lambda) {
  int a_size = alpha.size();
  int b_size = beta.size();
  int l_size = lambda.size();
  NumericVector rand_variates(n);
  RNGScope scope;
  NumericVector RUNI = runif(n, 0.0, 1.0);
  if (1 == a_size && a_size == b_size && b_size == l_size) {
    /* if parameters are all singletons (not vectors). then the idea is to find the largest CDF point in the vector
     * and build counts up to that point. Every other value is a lookup off of the largest vector. Otherwise, each
     * random variate has to be calculated individually which can be much slower.
     */
    NumericVector::iterator MX = std::max_element(RUNI.begin(), RUNI.end());
	  double maxquantile = *MX;
    NumericVector CDFVEC;
    int del_quantile = 0; //Will become "needed integer"
    CDFVEC.push_back(exp(-lambda[0]) / pow((1 + beta[0]), alpha[0])); //pre-load 0 value
    double cdf_curr_top = CDFVEC[0];
    while (cdf_curr_top < maxquantile) {
      ++del_quantile;
      CDFVEC.push_back(ddelap_C_S(del_quantile, alpha[0], beta[0], lambda[0], FALSE) +  CDFVEC[del_quantile - 1]);
      cdf_curr_top = CDFVEC[del_quantile];
    }
    NumericVector::iterator foundit;
    for (int i = 0; i < n; ++i) {
      if (RUNI[i] == 0) {
        rand_variates[i] = 0;
      } else {
        foundit = std::upper_bound (CDFVEC.begin(), CDFVEC.end(), RUNI[i]);
        double spot = foundit - CDFVEC.begin();
        rand_variates[i] = spot;
		  }
    }
  } else {
    rand_variates = qdelap_C(RUNI, alpha, beta, lambda, TRUE, FALSE);
    }
	return (rand_variates);
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