#include <Rcpp.h>
#include <vector>
#include <algorithm>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector ddelap_C(NumericVector x, double alpha, double beta, double lambda, bool lg) {
	int n = x.size();
	NumericVector DDLP(n);
	for (int t = 0; t < n; t++) {
		int k = ceil(x[t]);
		double PV = 0.0;
		for(int i = 0; i <= k; i++) {
			PV += exp(lgamma(alpha+i)+i*log(beta)+(k-i)*log(lambda)-lambda-lgamma(alpha)-lgamma(i+1)-(alpha+i)*log(1+beta)-lgamma(k-i+1));			
		}
		DDLP[t] = PV;
	}
		if (lg==TRUE) {
		for (int t = 0; t < n; t++) {
			DDLP[t] = log(DDLP[t]);
		}
	}
	return (DDLP);
}

// [[Rcpp::export]]
std::vector <double> pdelap_C(std::vector <double> q, double alpha, double beta, double lambda, bool lt, bool lp) {
// The idea is to find the largest value in the vector, and build the PDF up to that point.
// Every other value is a lookup off of the largest vector. This should make the routine
// take O(n) instead of O(n^2) time.	
	std::vector <double>::size_type n = q.size();
	std::vector <double>::iterator MX = std::max_element(q.begin(), q.end());
	int top = ceil(*MX);
	std::vector <double> PDLP(top+1);
	PDLP[0] = exp(-lambda)/pow((1+beta), alpha);
	for (int t = 1; t <= top; t++) {
		double PDFt = 0.0;
		for (int i = 0; i <= t; i++) {
			PDFt += exp(lgamma(alpha+i)+i*log(beta)+(t-i)*log(lambda)-lambda-lgamma(alpha)-lgamma(i+1)-(alpha+i)*log(1+beta)-lgamma(t-i+1));			
			}
		PDLP[t] = PDLP[t-1] + PDFt;
	}
	std::vector <double> PDLPRET(n);
	for (std::vector <double>::size_type t = 0; t < n; t++) {
		PDLPRET[t] = PDLP[ceil(q[t])];
	}
	if (lt==FALSE) {
		for (std::vector <double>::size_type t = 0; t < n; t++) {
			PDLPRET[t] = 1.0 - PDLPRET[t];
		}
	}
	if (lp==TRUE) {
		for (std::vector <double>::size_type t = 0; t < n; t++) {
			PDLPRET[t] = log(PDLPRET[t]);
		}
	}
	return (PDLPRET);
}

bool OUTSIDE01 (double quantile) {return (quantile <= 0.0 || quantile >= 1.0);}

// [[Rcpp::export]]
std::vector <double> qdelap_C(std::vector <double> p, double alpha, double beta, double lambda, bool lt, bool lp) {
// The idea is to find the largest CDF point in the vector, and build counts up to that point.
// Every other value is a lookup off of the largest vector. This should make the routine
// take O(n) instead of O(n^2) time.  
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
  std::vector <double>::iterator MX = std::max_element(pcopy.begin(), pcopy.end());
	double top = *MX;
	std::vector <double> CDFVEC;
	//PDLP[0] = exp(-lambda)/pow((1+beta), alpha);
	//for (int t = 1; t <= top; t++) {
	//	double PDFt = 0.0;
	//	for (int i = 0; i <= t; i++) {
	//		PDFt += exp(lgamma(alpha+i)+i*log(beta)+(t-i)*log(lambda)-lambda-lgamma(alpha)-lgamma(i+1)-(alpha+i)*log(1+beta)-lgamma(t-i+1));			
	//		}
	//	PDLP[t] = PDLP[t-1] + PDFt;
	//}
	//std::vector <double> PDLPRET(n);
	//for (std::vector <double>::size_type t = 0; t < n; t++) {
	//	PDLPRET[t] = PDLP[ceil(q[t])];
	//}

	
	return (pcopy2);
}