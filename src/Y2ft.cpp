#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// [[Rcpp::export]]
List Y2ft(int m, int n, List interval, NumericVector z, NumericMatrix Y, NumericMatrix dN){

	NumericMatrix Y2(m,n);
	NumericMatrix dN2(m,n);


	int i, j, h;

    for(i=0; i<m; i++){
        for(j=0; j<n; j++){
            Y2(i,j)=0;
            dN2(i,j)=0;
        }
    }

    for(i=0; i<n; i++){
        for(h=0; h<m; h++){
            SEXP interval2= interval[h];
            Rcpp::NumericVector interval3(interval2);
            if(interval3(0)<=z(i) && z(i)<interval3(1)){
                for(j=0; j<n; j++){
                    Y2(h,j) = Y2(h,j)+Y(i,j);
                    dN2(h,j)=dN2(h,j)+dN(i,j);
                }
            }
        }
    }

 	List ret;
	ret["Y2"] = Y2; ret["dN2"] = dN2;
	return Rcpp::wrap(ret);
}



