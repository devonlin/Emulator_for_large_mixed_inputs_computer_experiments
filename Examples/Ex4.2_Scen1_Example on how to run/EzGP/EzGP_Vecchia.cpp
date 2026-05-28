//#define BOOST_DISABLE_ASSERTS

#include <RcppArmadillo.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/gamma.hpp>
//#include "onepass.h"
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::depends(BH)]]
#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;
using namespace Rcpp;
using namespace arma;

// give p_e, q_e, m_e as arguments???????????????????????????????????????????????????????????
//'EzGP
// EzGP covariance function
// [[Rcpp::export]]
//add p_e,q_e,m_e
arma::mat EzGP_cov(
    arma::vec covparms, 
    arma::mat locs, 
    int p_e, 
    int q_e, 
    std::vector<int> m_e
) {
    int n = locs.n_rows;
    double tau = 1.490116e-08;  // Constant to add to the diagonal elements
    //int n_cols = locs.n_cols;
   // Function to calculate the sum of subvector
    auto sum_subvector = [](const std::vector<int>& vec, int start, int end) {
        return std::accumulate(vec.begin() + start, vec.begin() + end, 0);
    };

    // Function to calculate the sum of squares
    auto psum = [](const arma::vec& x1, const arma::vec& x2, const arma::vec& par2) {
        int n = x1.n_elem;
        double result = 0.0;
        for (int i = 0; i < n; ++i) {
            result += -par2[i] * pow((x1[i] - x2[i]), 2);
        }
        return result;
    };
    // Calculate the sum of m_e
    int sum_m_e = std::accumulate(m_e.begin(), m_e.end(), 0);

    // Calculate the total number of parameters
    int npar = p_e + q_e + p_e * sum_m_e;
   //Rcpp::Rcout << "npar: " << npar << std::endl;
    // Print n and n_cols
   // Rcpp::Rcout << "n: " << n << ", n_cols: " << n_cols << std::endl;

    // Print the parameters passed
   // Rcpp::Rcout << "p_e: " << p_e << ", q_e: " << q_e << std::endl;
   // Rcpp::Rcout << "m_e: ";
   //for (int val : m_e) {
    //    Rcpp::Rcout << val << " ";
   // }
   //Rcpp::Rcout << std::endl;

    // Calculate covariance matrix
    arma::mat covmat(n, n, arma::fill::zeros);

    for (int i1 = 0; i1 < n; i1++) {
        for (int i2 = 0; i2 <= i1; i2++) {
            arma::vec sigma0 = arma::vec{covparms[0]};
            arma::vec sigmah = covparms.subvec(1, q_e);
            arma::vec theta0 = covparms.subvec(q_e + 1, q_e + p_e);
            arma::vec thetah = covparms.subvec(q_e + p_e + 1, npar);
         // Rcpp::Rcout << "sigma0: " << sigma0 << ", sigmah: " << sigmah << std::endl;
         // Rcpp::Rcout << "theta0: " << theta0 << ", thetah: " << thetah << std::endl;


            // Extract x and z values for rows i1 and i2
            arma::vec x1 = locs.row(i1).subvec(0, p_e - 1).t();
            arma::vec z1 = locs.row(i1).subvec(p_e, p_e + q_e - 1).t();
            arma::vec x2 = locs.row(i2).subvec(0, p_e - 1).t();
            arma::vec z2 = locs.row(i2).subvec(p_e, p_e + q_e - 1).t();
   
            double res1 = sigma0[0] * exp(psum(x1, x2, theta0));
             
            for (int h = 0; h < q_e; ++h) {
                if (z1[h] != z2[h]) {
                    continue;
                }
                int l = z1[h];
                int start_index = sum_subvector(m_e, 0, h) * p_e + (l - 1) * p_e;
                int end_index = start_index + p_e - 1;
                arma::vec theta_segment = thetah.subvec(start_index, end_index);
                double sum_term = psum(x1, x2, theta_segment);
                res1 += sigmah[h] * exp(sum_term);
            }

            covmat(i1, i2) = res1;
            if (i1 != i2) {
                covmat(i2, i1) = res1;
            }
        }
    }
 covmat.diag() += tau;
    return covmat;
}

// give p_e, q_e, m_e as arguments???????????????????????????????????????????????????????????
// Derivatives of EzGP covariance
// [[Rcpp::export]]
arma::cube d_EzGP(
    arma::vec covparms, 
    arma::mat locs, 
    int p_e, 
    int q_e, 
    std::vector<int> m_e
) {
    int n = locs.n_rows;
    //int n_cols = locs.n_cols;
  // Function to calculate the sum of subvector
    auto sum_subvector = [](const std::vector<int>& vec, int start, int end) {
        return std::accumulate(vec.begin() + start, vec.begin() + end, 0);
    };

    // Function to calculate the sum of squares
    auto psum = [](const arma::vec& x1, const arma::vec& x2, const arma::vec& par2) {
        int n = x1.n_elem;
        double result = 0.0;
        for (int i = 0; i < n; ++i) {
            result += -par2[i] * pow((x1[i] - x2[i]), 2);
        }
        return result;
    };

    // Print the parameters passed
   // Rcpp::Rcout << "p_e: " << p_e << ", q_e: " << q_e << std::endl;
   // Rcpp::Rcout << "m_e: ";
   // for (int val : m_e) {
    //    Rcpp::Rcout << val << " ";
    //}
   // Rcpp::Rcout << std::endl;

    // Calculate the sum of m_e
    int sum_m_e = std::accumulate(m_e.begin(), m_e.end(), 0);

    // Calculate the total number of parameters
    int npar = p_e + q_e + p_e * sum_m_e;
   
    // Calculate derivatives
    arma::cube dcovmat(n, n, covparms.n_elem, arma::fill::zeros);

    for (int i1 = 0; i1 < n; i1++) {
        for (int i2 = 0; i2 <= i1; i2++) {
            arma::vec sigma0 = arma::vec{covparms[0]};
            arma::vec sigmah = covparms.subvec(1, q_e);
            arma::vec theta0 = covparms.subvec(q_e + 1, q_e + p_e);
            arma::vec thetah = covparms.subvec(q_e + p_e + 1, npar);

            // Extract x and z values for rows i1 and i2
            arma::vec x1 = locs.row(i1).subvec(0, p_e - 1).t();
            arma::vec z1 = locs.row(i1).subvec(p_e, p_e + q_e - 1).t();
            arma::vec x2 = locs.row(i2).subvec(0, p_e - 1).t();
            arma::vec z2 = locs.row(i2).subvec(p_e, p_e + q_e - 1).t();

            // Derivative w.r.t. sigma_0
            dcovmat(i1, i2, 0) = exp(psum(x1, x2, theta0));

            // Derivative w.r.t. sigma_h
            for (int h = 0; h < q_e; h++) { 
                double gx = 0;
                if (z1[h] == z2[h]) {
                    int l = z1[h];
                    int start_index = sum_subvector(m_e, 0, h) * p_e + (l - 1) * p_e;
                    int end_index = start_index + p_e - 1;
                    arma::vec theta_segment = thetah.subvec(start_index, end_index);
                    gx = exp(psum(x1, x2, theta_segment));
                }
                dcovmat(i1, i2, 1 + h) = gx;
            }

            // Derivative w.r.t. theta_0_s
            for (int s = 0; s < p_e; s++) { 
                dcovmat(i1, i2, 1 + q_e + s) = -sigma0[0] * pow((x1[s] - x2[s]), 2) * exp(psum(x1, x2, theta0));
            }

            // Derivative w.r.t. theta_h
            for (int h = 0; h < q_e; h++) { 
                if (z1[h] != z2[h]) {
                    continue;
                }
                int l = z1[h];
                for (int s = 0; s < p_e; s++) { 
                    int start_index = sum_subvector(m_e, 0, h) * p_e + (l - 1) * p_e;
                    int end_index = start_index + p_e - 1;
                    arma::vec theta_segment = thetah.subvec(start_index, end_index);
                    dcovmat(i1, i2, 1 + q_e + p_e + start_index + s) = -sigmah[h] * pow((x1[s] - x2[s]), 2) * exp(psum(x1, x2, theta_segment));
                }
            }

            if (i1 != i2) {
                for (int p = 0; p < dcovmat.n_slices; ++p) {
                    dcovmat(i2, i1, p) = dcovmat(i1, i2, p);
                }
            }
        }
    }

    return dcovmat;
}




void get_covfun(std::string covfun_name_string, 
                arma::mat (*p_covfun[1])(arma::vec, arma::mat, int, int, std::vector<int>), 
                arma::cube (*p_d_covfun[1])(arma::vec, arma::mat, int, int, std::vector<int>)) 
{
    if (covfun_name_string.compare("EzGP_cov") == 0) { 
        p_covfun[0] = EzGP_cov;
        p_d_covfun[0] = d_EzGP;
    } else { 
        Rcpp::Rcout << "Unrecognized Covariance Function Name \n";
    }
}




arma::vec forward_solve( arma::mat cholmat, arma::vec b ){

    int n = cholmat.n_rows;
    arma::vec x(n);
    x(0) = b(0)/cholmat(0,0);

    for(int i=1; i<n; i++){
        double dd = 0.0;
        for(int j=0; j<i; j++){
            dd += cholmat(i,j)*x(j);
        }
        x(i) = (b(i)-dd)/cholmat(i,i);
    }    
    return x;

} 

arma::mat forward_solve_mat( arma::mat cholmat, arma::mat b ){

    int n = cholmat.n_rows;
    int p = b.n_cols;
    arma::mat x(n,p);
    for(int k=0; k<p; k++){ x(0,k) = b(0,k)/cholmat(0,0); }

    for(int i=1; i<n; i++){
	for(int k=0; k<p; k++){
            double dd = 0.0;
            for(int j=0; j<i; j++){
                dd += cholmat(i,j)*x(j,k);
            }
            x(i,k) = (b(i,k)-dd)/cholmat(i,i);
       	}
    }    
    return x;
} 

arma::vec backward_solve( arma::mat lower, arma::vec b ){

    int n = lower.n_rows;
    arma::vec x(n);
    x(n-1) = b(n-1)/lower(n-1,n-1);

    for(int i=n-2; i>=0; i--){
        double dd = 0.0;
        for(int j=n-1; j>i; j--){
            dd += lower(j,i)*x(j);
        }
        x(i) = (b(i)-dd)/lower(i,i);
    }    
    return x;
} 

arma::mat backward_solve_mat( arma::mat cholmat, arma::mat b ){

    int n = cholmat.n_rows;
    int p = b.n_cols;
    arma::mat x(n,p);
    for(int k=0; k<p; k++){ x(n-1,k) = b(n-1,k)/cholmat(n-1,n-1); }

    for(int i=n-2; i>=0; i--){
	for(int k=0; k<p; k++){
            double dd = 0.0;
            for(int j=n-1; j>i; j--){
                dd += cholmat(j,i)*x(j,k);
            }
            x(i,k) = (b(i,k)-dd)/cholmat(i,i);
       	}
    }    
    return x;
} 


arma::mat mychol( arma::mat A ){

    arma::uword n = A.n_rows;
    arma::mat L(n,n);
    bool pd = true;
    
    // upper-left entry
    if( A(0,0) < 0 ){
	pd = false;
	L(0,0) = 1.0;
    } else {
        L(0,0) = std::sqrt(A(0,0));
    }
    if( n > 1 ){
	// second row
	L(1,0) = A(1,0)/L(0,0);
	double f = A(1,1) - L(1,0)*L(1,0);
	if( f < 0 ){
	    pd = false;
	    L(1,1) = 1.0;
	} else {
	    L(1,1) = std::sqrt( f );
	}
	// rest of the rows
	if( n > 2 ){
            for(uword i=2; i<n; i++){
    	        // leftmost entry in row i
    	        L(i,0) = A(i,0)/L(0,0);
    	        // middle entries in row i 
    	        for(uword j=1; j<i; j++){
    	            double d = A(i,j);
    	            for(uword k=0; k<j; k++){
    	        	d -= L(i,k)*L(j,k);
    	            }
    	            L(i,j) = d/L(j,j);
    	        }
		// diagonal entry in row i
    	        double e = A(i,i);
    	        for(uword k=0; k<i; k++){
    	            e -= L(i,k)*L(i,k);
    	        }
		if( e < 0 ){
		    pd = false;
		    L(i,i) = 1.0;
		} else {
    	            L(i,i) = std::sqrt(e);
		}
	    }
	}
    }
    return L;	
}



// give p_e, q_e, m_e as arguments???????????????????????????????????????????????????????????
void compute_pieces_modif(
    arma::vec covparms, 
    StringVector covfun_name,
    arma::mat locs, 
    arma::mat NNarray,
    arma::vec y, 
    arma::mat X,
    mat* XSX,
    vec* ySX,
    double* ySy,
    double* logdet,
    cube* dXSX,
    mat* dySX,
    vec* dySy,
    vec* dlogdet,
    mat* ainfo,
    int profbeta,
    int grad_info,
    int p_e,  // new argument
    int q_e,  // new argument
    std::vector<int> m_e  // new argument
) {
    // data dimensions
    int n = y.n_elem;
    int m = NNarray.n_cols;
    int p = X.n_cols;
    int nparms = covparms.n_elem;
    int dim = locs.n_cols;
    
    // convert StringVector to std::string to use .compare() below
    std::string covfun_name_string;
    covfun_name_string = covfun_name[0];
    
    // assign covariance fun and derivative based on covfun_name_string
    mat (*p_covfun[1])(arma::vec, arma::mat, int, int, std::vector<int>);
    cube (*p_d_covfun[1])(arma::vec, arma::mat, int, int, std::vector<int>);
    get_covfun(covfun_name_string, p_covfun, p_d_covfun);
  
    // ... other code remains the same ...

#pragma omp parallel 
{   
    arma::mat l_XSX = arma::mat(p, p, fill::zeros);
    arma::vec l_ySX = arma::vec(p, fill::zeros);
    double l_ySy = 0.0;
    double l_logdet = 0.0;
    arma::cube l_dXSX = arma::cube(p,p, nparms, fill::zeros);
    arma::mat l_dySX = arma::mat(p, nparms, fill::zeros);
    arma::vec l_dySy = arma::vec(nparms, fill::zeros);
    arma::vec l_dlogdet = arma::vec(nparms, fill::zeros);
    arma::mat l_ainfo = arma::mat(nparms, nparms, fill::zeros);

    #pragma omp for	    
    for(int i=0; i<n; i++){
    
        int bsize = std::min(i+1,m);

        // first, fill in ysub, locsub, and X0 in reverse order
        arma::mat locsub(bsize, dim);
        arma::vec ysub(bsize);
        arma::mat X0( bsize, p );
        for(int j=bsize-1; j>=0; j--){
            ysub(bsize-1-j) = y( NNarray(i,j)-1 );
            for(int k=0;k<dim;k++){ locsub(bsize-1-j,k) = locs( NNarray(i,j)-1, k ); }
            if(profbeta){ 
                for(int k=0;k<p;k++){ X0(bsize-1-j,k) = X( NNarray(i,j)-1, k ); } 
            }
        }
        
        // compute covariance matrix and derivatives and take cholesky
        arma::mat covmat = p_covfun[0](covparms, locsub, p_e, q_e, m_e);  // updated call

        arma::cube dcovmat;
        if(grad_info){ 
            dcovmat = p_d_covfun[0](covparms, locsub, p_e, q_e, m_e);  // updated call
        }
		
        arma::mat cholmat = eye( size(covmat) );
        chol( cholmat, covmat, "lower" );
        // Print the dimensions and contents of covmat


        // i1 is conditioning set, i2 is response        
        //arma::span i1 = span(0,bsize-2);
        arma::span i2 = span(bsize-1,bsize-1);
       
        // get last row of cholmat
        arma::vec onevec = zeros(bsize);
        onevec(bsize-1) = 1.0;
        arma::vec choli2;
        if(grad_info){
            choli2 = backward_solve( cholmat, onevec );
        }
        
        bool cond = bsize > 1;
        
        // do solves with X and y
        arma::mat LiX0;
        if(profbeta){
            LiX0 = forward_solve_mat( cholmat, X0 );
        }

        arma::vec Liy0 = forward_solve( cholmat, ysub );
        
        // loglik objects
        l_logdet += 2.0*std::log( as_scalar(cholmat(i2,i2)) ); 
        l_ySy +=    pow( as_scalar(Liy0(i2)), 2 );
        if(profbeta){
            l_XSX +=   LiX0.rows(i2).t() * LiX0.rows(i2);
            l_ySX += ( Liy0(i2) * LiX0.rows(i2) ).t();
        }
        
        if(grad_info){
            arma::mat LidSLi2(bsize,nparms);
            
            if(cond){ // if we condition on anything
                
                for(int j=0; j<nparms; j++){
                    arma::vec LidSLi3 = forward_solve( cholmat, dcovmat.slice(j) * choli2 );
                    arma::vec v1 = LiX0.t() * LidSLi3;
                    double s1 = as_scalar( Liy0.t() * LidSLi3 ); 
                    (l_dXSX).slice(j) += v1 * LiX0.rows(i2) + ( v1 * LiX0.rows(i2) ).t() - 
                        as_scalar(LidSLi3(i2)) * ( LiX0.rows(i2).t() * LiX0.rows(i2) );
                    (l_dySy)(j) += as_scalar( 2.0 * s1 * Liy0(i2)  - 
                        LidSLi3(i2) * Liy0(i2) * Liy0(i2) );
                    (l_dySX).col(j) += (  s1 * LiX0.rows(i2) + ( v1 * Liy0(i2) ).t() -  
                        as_scalar( LidSLi3(i2) ) * LiX0.rows(i2) * as_scalar( Liy0(i2))).t();
                    (l_dlogdet)(j) += as_scalar( LidSLi3(i2) );
                    LidSLi2.col(j) = LidSLi3;
                }

                for(int i=0; i<nparms; i++){ for(int j=0; j<i+1; j++){
                    (l_ainfo)(i,j) += 
                        1.0*accu( LidSLi2.col(i) % LidSLi2.col(j) ) - 
                        0.5*accu( LidSLi2.rows(i2).col(j) %
                                  LidSLi2.rows(i2).col(i) );
                }}
                
            } else { // similar calculations, but for when there is no conditioning set
                for(int j=0; j<nparms; j++){
                    arma::mat LidSLi = forward_solve_mat( cholmat, dcovmat.slice(j) );
                    LidSLi = forward_solve_mat( cholmat, LidSLi.t() );
                    (l_dXSX).slice(j) += LiX0.t() *  LidSLi * LiX0; 
                    (l_dySy)(j) += as_scalar( Liy0.t() * LidSLi * Liy0 );
                    (l_dySX).col(j) += ( ( Liy0.t() * LidSLi ) * LiX0 ).t();
                    (l_dlogdet)(j) += trace( LidSLi );
                    LidSLi2.col(j) = LidSLi;
                }
                
                for(int i=0; i<nparms; i++){ for(int j=0; j<i+1; j++){
                    (l_ainfo)(i,j) += 0.5*accu( LidSLi2.col(i) % LidSLi2.col(j) ); 
                }}

            }
        }
 
    }

#pragma omp critical
{
    *XSX += l_XSX;
    *ySX += l_ySX;
    *ySy += l_ySy;
    *logdet += l_logdet;
    *dXSX += l_dXSX;
    *dySX += l_dySX;
    *dySy += l_dySy;
    *dlogdet += l_dlogdet;
    *ainfo += l_ainfo;
}
}
}

    
// give p_e, q_e, m_e as arguments???????????????????????????????????????????????????????????
void synthesize_modif(
    NumericVector covparms, 
    StringVector covfun_name,
    const NumericMatrix locs, 
    NumericMatrix NNarray,
    NumericVector& y, 
    NumericMatrix X,
    NumericVector* ll, 
    NumericVector* betahat,
    NumericVector* grad,
    NumericMatrix* info,
    NumericMatrix* betainfo,
    bool profbeta,
    bool grad_info,
    int p_e,  // new argument
    int q_e,  // new argument
    std::vector<int> m_e  // new argument
) {

    // data dimensions
    int n = y.length();
    int p = X.ncol();
    int nparms = covparms.length();
    
    // likelihood objects
    arma::mat XSX = arma::mat(p, p, fill::zeros);
    arma::vec ySX = arma::vec(p, fill::zeros);
    double ySy = 0.0;
    double logdet = 0.0;
    
    // gradient objects    
    arma::cube dXSX = arma::cube(p,p,nparms,fill::zeros);
    arma::mat dySX = arma::mat(p, nparms, fill::zeros);
    arma::vec dySy = arma::vec(nparms, fill::zeros);
    arma::vec dlogdet = arma::vec(nparms, fill::zeros);
    // fisher information
    arma::mat ainfo = arma::mat(nparms, nparms, fill::zeros);

    // this is where the big computation happens
    // first convert Numeric- to arma
    arma::vec covparms_c = arma::vec(covparms.begin(),covparms.length());
    arma::mat locs_c = arma::mat(locs.begin(),locs.nrow(),locs.ncol());
    arma::mat NNarray_c = arma::mat(NNarray.begin(),NNarray.nrow(),NNarray.ncol());
    arma::vec y_c = arma::vec(y.begin(),y.length());
    arma::mat X_c = arma::mat(X.begin(),X.nrow(),X.ncol());

    // give p_e, q_e, m_e as arguments
    compute_pieces_modif(
        covparms_c, covfun_name, locs_c, NNarray_c, y_c, X_c,
        &XSX, &ySX, &ySy, &logdet, &dXSX, &dySX, &dySy, &dlogdet, &ainfo,
        profbeta, grad_info,
        p_e,  // pass the new argument
        q_e,  // pass the new argument
        m_e   // pass the new argument
    );
        
    // synthesize everything and update loglik, grad, beta, betainfo, info
    
    // betahat and dbeta
    arma::vec abeta = arma::vec(p, fill::zeros);
    if(profbeta){ abeta = solve(XSX, ySX);
    //std::cout << "abeta: " << abeta.t();
    }
    for(int j=0; j<p; j++){ (*betahat)(j) = abeta(j); }

    arma::mat dbeta = arma::mat(p,nparms, fill::zeros);
    if(profbeta && grad_info){
        for(int j=0; j<nparms; j++){
            dbeta.col(j) = solve(XSX, dySX.col(j) - dXSX.slice(j) * abeta);
        }
    }

    // get sigmahatsq
    double sig2 = (ySy - 2.0 * as_scalar(ySX.t() * abeta) + 
                   as_scalar(abeta.t() * XSX * abeta)) / n;
    
    // loglikelihood
    (*ll)(0) = -0.5 * (n * std::log(2.0 * M_PI) + logdet + n * sig2); 
    
    if(profbeta){
        // betainfo
        for(int i=0; i<p; i++){
            for(int j=0; j<i+1; j++){
                (*betainfo)(i,j) = XSX(i,j);
                (*betainfo)(j,i) = XSX(j,i);
            }
        }
    }

    if(grad_info){
        // gradient
        for(int j=0; j<nparms; j++){
            (*grad)(j) = 0.0;
            (*grad)(j) -= 0.5 * dlogdet(j);
            (*grad)(j) += 0.5 * dySy(j);
            (*grad)(j) -= as_scalar(abeta.t() * dySX.col(j));
            (*grad)(j) += as_scalar(ySX.t() * dbeta.col(j));
            (*grad)(j) += 0.5 * as_scalar(abeta.t() * dXSX.slice(j) * abeta);
            (*grad)(j) -= as_scalar(abeta.t() * XSX * dbeta.col(j));
        }

        // fisher information
        for(int i=0; i<nparms; i++){
            for(int j=0; j<i+1; j++){
                (*info)(i,j) = ainfo(i,j);
                (*info)(j,i) = (*info)(i,j);
            }
        }
    }

}
    
 //' Vecchia's loglikelihood, gradient, and Fisher information
//'
//' This function returns Vecchia's (1988) approximation to the Gaussian
//' loglikelihood, profiling out the regression coefficients, and returning
//' the gradient and Fisher information. 
//' Vecchia's approximation modifies the ordered conditional
//' specification of the joint density; rather than each term in the product
//' conditioning on all previous observations, each term conditions on
//' a small subset of previous observations.
//' @inheritParams vecchia_meanzero_loglik
//' @param X Design matrix of covariates. Row \code{i} of \code{X} contains
//' the covariates for the observation at row \code{i} of \code{locs}.
//' @return A list containing 
//' \itemize{
//'     \item \code{loglik}: the loglikelihood
//'     \item \code{grad}: gradient with respect to covariance parameters
//'     \item \code{info}: Fisher information for covariance parameters
//'     \item \code{betahat}: profile likelihood estimate of regression coefs
//'     \item \code{betainfo}: information matrix for \code{betahat}.
//' }
//' The covariance matrix for \code{$betahat} is the inverse of \code{$betainfo}.
//' @examples
//' n1 <- 20
//' n2 <- 20
//' n <- n1*n2
//' locs <- as.matrix( expand.grid( (1:n1)/n1, (1:n2)/n2 ) )
//' X <- cbind(rep(1,n),locs[,2])
//' covparms <- c(2, 0.2, 0.75, 0)
//' y <- X %*% c(1,2) + fast_Gp_sim(covparms, "matern_isotropic", locs, 50 )
//' ord <- order_maxmin(locs)
//' NNarray <- find_ordered_nn(locs,20)
//' #loglik <- vecchia_profbeta_loglik_grad_info( covparms, "matern_isotropic", 
//' #    y, X, locs, NNarray )
//' @export
// [[Rcpp::export]]
// give p_e, q_e, m_e as arguments???????????????????????????????????????????????????????????
List vecchia_profbeta_loglik_grad_info_modif( 
    NumericVector covparms, 
    StringVector covfun_name,
    NumericVector y,
    NumericMatrix X,
    const NumericMatrix locs,
    NumericMatrix NNarray,
    int p_e,  // new argument
    int q_e,  // new argument
    std::vector<int> m_e  // new argument
) {
    NumericVector ll(1);
    NumericVector grad(covparms.length());
    NumericVector betahat(X.ncol());
    NumericMatrix info(covparms.length(), covparms.length());
    NumericMatrix betainfo(X.ncol(), X.ncol());

    // this function calls synthesize_modif
    // then synthesizes the result into loglik, beta, grad, info, betainfo
    synthesize_modif(covparms, covfun_name, locs, NNarray, y, X,
        &ll, &betahat, &grad, &info, &betainfo, true, true,
        p_e,  // pass the new argument
        q_e,  // pass the new argument
        m_e   // pass the new argument
    );
    
    List ret = List::create(Named("loglik") = ll, Named("betahat") = betahat,
                            Named("grad") = grad, Named("info") = info, Named("betainfo") = betainfo);
    return ret;
}