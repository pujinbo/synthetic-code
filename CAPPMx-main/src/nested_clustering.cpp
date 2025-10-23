#include <RcppArmadillo.h>
#include <RcppDist.h>
#include <RcppGSL.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include<omp.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_psi.h>
#include<trunclst.h>
#include <limits>
// [[Rcpp::depends(RcppArmadillo,RcppDist)]]
// [[Rcpp::depends(RcppGSL)]]
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::plugins(openmp)]]

using namespace arma;

#define crossprod(x) symmatu(x.t() * x)
#define tcrossprod(x) symmatu(x * x.t())
#define ssq(x) dot(x,x)

struct eta_list{
  arma::vec eta,sum_of_samples;
};

eta_list gen_eta(int nj,double niw_kap,arma::vec sum_of_samples, arma::vec Meta,
                 const arma::mat &Delta_inv, const arma::mat &chol_Delta,const arma::mat &chol_precmat_eta ){
  unsigned k=sum_of_samples.n_elem;
  
  double niw_kap_post=nj-1+niw_kap;
  // sum_of_samples-=eta;
  arma::vec marg_mu= sum_of_samples/niw_kap_post;
  
  arma::vec z1=arma::vec(k, fill::randn),z2=arma::vec(k, fill::randn);
  
  Meta+=( Delta_inv*  marg_mu  +solve( trimatu(chol_Delta ) ,z2)/sqrt(niw_kap_post) ) ; ///since chol_deta is generated to be lower-triangular
  
  eta_list L;
  L.eta = solve(trimatu( chol_precmat_eta), z1 +solve(trimatl((chol_precmat_eta).t()),  Meta ) ) ;
  
  L.sum_of_samples= sum_of_samples+L.eta;
  return L;
}


inline double log_lomax(double x, double a, double b){
  return log(a/b) -(a+1)* log1p(x/b);
}

inline double surv_fn(double x,unsigned nu, double a, double b){//nu=1 implies failure time
  return nu ?  (log(a/b) -(a+1)* log1p(x/b)) : ( -a* log1p(x/b));
}

double c_lmvgamma (double x, int p) {
  int i;
  double ans = 0;
  if (p < 1)
    Rcpp::stop("p must be greater than or equal to 1.");
  if (x <= 0)
    Rcpp::stop("x must be greater than 0.");
  ans =(p * (p - 1)/4.0) * M_LNPI;
  for (i = 0; i < p; i++){
    ans +=  (lgamma(x  - (i/2.0) ));
  }
  return ans;
}

double log_marg_dens (const arma::mat &chol_marginal_delta,int n,double t1ppp
                        ,double  niw_nu,double  niw_kap,double  diag_psi_iw   ){
  int k=chol_marginal_delta.n_cols;
  double dens;
  dens=-(k/2)*(n* M_LNPI+log1p(n/niw_kap)-niw_nu*log(diag_psi_iw)- (niw_nu+n)*log(t1ppp) ) 
    +c_lmvgamma ( (niw_nu+n)/2 ,  k)- c_lmvgamma ( niw_nu/2 ,  k)- (niw_nu+n)*sum(log(chol_marginal_delta.diag()) );
  return dens;
}

//'  The natural logarithm of the sum of the exponentials of the arguments
//'
//' A numerically stable version of \code{log(sum(exp(x)))}
//' @param x a numeric vector
//'
//' @return The natural logarithm of sum of \code{exp(x)}
//' @export
//' @examples
//' x=c(-1000,-1001)
//' log(sum(exp(x))) ##naive implementation
//' log_sum_exp(x) ##numerically stable implementation
// [[Rcpp::export]]
double log_sum_exp(const arma::vec &x) 
{
  arma::uword max_ind=x.index_max();
  double maxVal= x(max_ind);
  
  double sum_exp=0.0;
  
  for (unsigned  i = 0; i < x.n_elem ; ++i){
    if(i !=max_ind)
      sum_exp += exp(  (x(i) - maxVal));
  }
  
  return log1p(sum_exp)+maxVal ;
}

inline double multiplier_fn(double niw_nu, double niw_kap, int k){
  return (niw_kap +1) / ( niw_kap* (niw_nu-k+1 )  );
}  

inline double log_t_density(double nu,  arma::vec x,arma::vec mu, arma::mat lower_chol){
  int k=x.n_elem;
  double det_sig_half= sum(log(lower_chol.diag() ) );
  arma::vec resid=solve(trimatl ( lower_chol ) , x - mu );
  
  double quad_form=dot(resid,resid);
  double density = lgamma((nu+k)/2 ) -lgamma(nu/2) - 
    (k* log((datum::pi)*nu) + (nu+k)*log1p(quad_form/nu) )/2 -det_sig_half;
  
  return density;
}

inline double log_t_density_empty(double nu, double kap, double diag_psi_iw, arma::vec x){
  int k=x.n_elem;
  double t1ppp=(kap+1)/(kap*nu);
  
  double tmp=t1ppp*diag_psi_iw;
  
  double quad_form=dot(x,x ); 
  double density = lgamma((nu+k)/2 ) -lgamma(nu/2) - 
    (k* log((datum::pi)*nu) + (nu+k)*log1p(quad_form/(tmp*nu)) +k*log(t1ppp*diag_psi_iw) )/2 ;
  // Rcpp::Rcout<<"t1ppp="<<t1ppp<<" det_sig_half"<<det_sig_half<<" quad_form="<<quad_form<<" density="<<density;
  return density;
}

arma::uvec std_setdiff(arma::uvec a, arma::uvec b) {
  
  // std::vector<int> a = arma::conv_to< std::vector<int> >::from(arma::sort(x));
  // std::vector<int> b = arma::conv_to< std::vector<int> >::from(arma::sort(y));
  a=sort(a);b=sort(b);
  
  std::vector<unsigned> out;
  
  std::set_difference(a.begin(), a.end(), b.begin(), b.end(),
                      std::inserter(out, out.end()));
  
  return arma::conv_to<arma::uvec>::from(out);
}

void calculate_prob_for_current_ind(arma::vec &tmp_marginal_mu,arma::vec &tmp_sum_of_samples,
                                    arma::mat &tmp_sum_of_sq,arma::mat &tmp_chol_marginal_delta,arma::mat &tmp_marginal_delta,
                                    const arma::mat &eta,const arma::vec &niw_nu_post,const  arma::vec &niw_kap_post, arma::vec &df_post,
                                    const arma::mat &sum_of_samples, const cube &sum_of_squares,
                                    const arma::uvec &del,const unsigned jj, const arma::uvec &nj_val1,
                                    const double dir_prec, const double diag_psi_iw,
                                    double &dens, double &cluster_prob){
  
  unsigned k=eta.n_cols,j=del(jj);
  arma::mat cross_prod=crossprod(eta.row(jj));
  tmp_sum_of_samples=sum_of_samples .col(j) - (eta.row(jj)).t() ;
  tmp_marginal_mu=(tmp_sum_of_samples ) /(niw_kap_post(j)-1) ;
  
  tmp_sum_of_sq=symmatu( sum_of_squares.slice(j)- cross_prod );
  
  // Rcpp::Rcout<<"eig_sym( tmp_sum_of_sq )=";
  // eig_sym( tmp_sum_of_sq ).print();
  
  // Rcpp::Rcout<<"jj="<<jj<<endl;
  // Rcpp::Rcout<<"Flag 1 current_ind="<<current_ind<<"nj_val-1="<<nj_val(current_ind) -1<<endl;
  
  
  tmp_marginal_delta=symmatu( tmp_sum_of_sq -  tcrossprod(tmp_sum_of_samples)/(niw_kap_post(j)-1 ) );
  // Rcpp::Rcout <<"tmp_marginal_delta.diag()="<< tmp_marginal_delta.diag() << endl;
  tmp_marginal_delta.diag() += diag_psi_iw;  
  
  /////////////////////////////////////////////////////////////////////////
  /*arma::uvec indsss=find(del==j);
   arma::mat tm_ssq=crossprod(eta.rows(indsss)) ;
   arma::vec tm_sum=sum(eta.rows(indsss),0).t();
   
   if( norm(sum_of_samples.col(j)- tm_sum,"inf")>.1){
   Rcpp::Rcout<<"j= "<<j<<" jj = "<<jj<<" nj_val1="<<nj_val1(j)<< endl;
   Rcpp::Rcout<<"sum diff"<<endl;
   Rcpp::stop("");
   }
   
   if(norm(( tm_ssq-sum_of_squares.slice(j) ),"inf")>.1 ){
   Rcpp::Rcout<<"j= "<<j<<" jj = "<<jj<<" nj_val1="<<nj_val1(j)<< endl;
   Rcpp::Rcout<<"ssq diff"<<endl;
   Rcpp::stop("");
   }
   
   arma::uvec indssa(1); indssa(0)=jj;
   indsss=std_setdiff(indsss,indssa);
   // Rcpp::Rcout<<"sort(indsss)"<<sort(indsss).t()<<endl;
   tm_ssq=  crossprod(eta.rows(indsss)) ;
   // tm_ssq.print("tm_ssq before: ");
   tm_ssq-=tmp_sum_of_sq;
   // tm_ssq.print("tm_ssq after: ");
   // tm_ssq=arma::abs(tm_ssq);
   if(norm(( tm_ssq),"inf")>.1){
   Rcpp::Rcout<<"j= "<<j<<" jj = "<<jj<<endl;
   Rcpp::stop("");
   }
   arma::vec tm_sum= sum ( eta.rows(indsss) ,0).t(), tm_mn=mean ( eta.rows(indsss) ,0).t();
   arma::mat tm_marginal_delta =symmatu( tm_ssq - tcrossprod(tm_sum)/(niw_kap_post(j)-1) );
   unsigned tm_n=indsss.n_elem;
   // arma::mat tm_marginal_delta= cov(eta.rows(indsss),1)* tm_n+ ((niw_kap_post(j)-1)*tm_n/(niw_kap_post(j)-1+tm_n)) *tcrossprod(tm_mn);
   tm_marginal_delta.diag() += diag_psi_iw;
   Rcpp::Rcout<<"tm_marginal_delta.is_sympd()= "<<tm_marginal_delta.is_sympd()<<endl;
   tm_marginal_delta.print("tm_marginal_delta before");
   tm_marginal_delta-=tmp_marginal_delta;
   tm_marginal_delta.print("tm_marginal_delta");*/
  /////////////////////////////////////////////////////////////////////////
  
  double t1ppp= multiplier_fn( niw_nu_post(j)-1 , niw_kap_post(j)-1, k);
  // Rcpp::Rcout<<"chol flag calculate"<<endl;
  // tmp_marginal_delta.print("tmp_marginal_delta :");
  // Rcpp::Rcout<<tmp_marginal_delta.is_sympd()<<" t1ppp= "<<t1ppp<< endl;
  tmp_chol_marginal_delta=sqrt(t1ppp)*  chol( tmp_marginal_delta ,"lower") ;
  // Rcpp::Rcout<<"chol flag calculate 2"<<endl;
  
  dens= log_t_density(df_post(j)-1 ,  (eta.row(jj)).t(),tmp_marginal_mu, tmp_chol_marginal_delta );
  cluster_prob=log( nj_val1(j)-1 +dir_prec);
}

void inline update_cluster_params (const arma::vec &tmp_marginal_mu,const arma::vec &tmp_sum_of_samples,
                                   const arma::mat &tmp_sum_of_sq,const arma::mat &tmp_chol_marginal_delta,const arma::mat &tmp_marginal_delta,
                                   const arma::mat &eta, arma::vec &niw_nu_post,arma::vec &niw_kap_post, arma::vec &df_post,
                                   arma::mat &sum_of_samples,arma::mat &marginal_mu, cube &sum_of_squares, cube &tmpcube, cube &marginal_delta,
                                   const arma::uvec &del,const unsigned jj,const unsigned current_ind, arma::uvec &nj_val,
                                   const double dir_prec, const double diag_psi_iw){
  
  //**** update the source cluster parameters ****//
  // arma::uvec indsss;  arma::mat tm_ssq;  arma::vec tm_sum;
  
  unsigned j=current_ind, k=eta.n_cols;
  --nj_val(j);
  // log_nj_val(j ) = log( --nj_val(j) +dir_prec);
  niw_kap_post(j)--; niw_nu_post(j)--; df_post(j)--;
  
  if(0){
    (sum_of_samples.col(j)).zeros();
    (sum_of_squares.slice(j)).zeros();
    marginal_mu.col(j).zeros();
    (marginal_delta.slice(j)).zeros();
    (tmpcube.slice(j)).zeros();
  }
  if(nj_val(j)>=1){
    // Rcpp::Rcout<<"update flag 1"<<endl;
    sum_of_samples.col(j)=tmp_sum_of_samples;
    // (tmp_sum_of_samples.t()).print("tmp_sum_of_samples");
    // Rcpp::Rcout<<size(sum_of_samples)<<"\t"<<size(tmp_sum_of_samples)<<endl;
    // Rcpp::Rcout<<size(sum_of_squares)<<"\t"<<size(tmp_sum_of_sq)<<endl;
    sum_of_squares.slice(j)=tmp_sum_of_sq;
    // tmp_sum_of_sq.print("tmp_sum_of_sq");
    
    marginal_mu.col(j)=tmp_marginal_mu;
    // (tmp_marginal_mu.t()).print("tmp_marginal_mu");
    
    // Rcpp::Rcout<<"update flag 1 1"<<endl;
    marginal_delta.slice(j)=tmp_marginal_delta;
    // tmp_marginal_delta.print("tmp_marginal_delta");
    // Rcpp::Rcout<<"update flag 1 2"<<endl;
    // chol_marginal_delta.slice(j)=tmp_chol_marginal_delta;
    tmpcube.slice(j)=tmp_chol_marginal_delta;
    // tmp_chol_marginal_delta.print("tmp_chol_marginal_delta");
    // Rcpp::Rcout<<"update flag 2"<<endl;
    
    /*indsss=find(del==j);
     tm_ssq=crossprod(eta.rows(indsss)) ;
     tm_sum=sum(eta.rows(indsss),0).t();
     
     if( norm(sum_of_samples.col(j)- tm_sum,"inf")>.1){
     Rcpp::Rcout<<"j= "<<j<<" jj = "<<jj<<" nj_val1="<<nj_val(j)<< endl;
     Rcpp::Rcout<<"sum diff source"<<endl;
     Rcpp::stop("");
     }
     
     if(norm(( tm_ssq-sum_of_squares.slice(j) ),"inf")>.1 ){
     Rcpp::Rcout<<"j= "<<j<<" jj = "<<jj<<" nj_val1="<<nj_val(j)<< endl;
     Rcpp::Rcout<<"ssq diff source"<<endl;
     Rcpp::stop("");
     }*/
  }
  
  
  ///////////////////////////////////////////////////////
  
  //**** update the destination cluster parameters ****//
  j=del(jj);
  ++nj_val(j);// log_nj_val(j ) = log( ++nj_val(j) +dir_prec);
  niw_kap_post(j)++; niw_nu_post(j)++; df_post(j)++;
  
  arma::mat cross_prod=crossprod(eta.row(jj));
  if(nj_val(j)==1){
    // Rcpp::Rcout<<"update flag 3"<<endl;
    sum_of_samples.col(j)=  (eta.row(jj)).t() ;
    sum_of_squares.slice(j) = cross_prod ;
    // Rcpp::Rcout<<"update flag 4"<<endl;
  }
  else{
    // Rcpp::Rcout<<"update flag 5"<<endl;
    sum_of_samples.col(j)+=  (eta.row(jj)).t() ;
    sum_of_squares.slice(j) += cross_prod ;
    // Rcpp::Rcout<<"update flag 6"<<endl;
  }
  
  marginal_mu.col(j) =(sum_of_samples.col(j) ) /(niw_kap_post(j)) ;
  // Rcpp::Rcout<<"update flag 7"<<endl;
  marginal_delta.slice(j) =symmatu( sum_of_squares.slice(j) - tcrossprod(sum_of_samples .col(j))/niw_kap_post(j) );
  (marginal_delta.slice(j)).diag() += diag_psi_iw;
  double t1ppp= multiplier_fn( niw_nu_post(j) , niw_kap_post(j), k);
  
  tmpcube.slice(j) =  chol(t1ppp*marginal_delta.slice(j) ,"lower") ;
  // Rcpp::Rcout<<"update flag 8"<<endl;
  // tmpcube.slice(j)=chol_marginal_delta.slice(j);
  
  /*indsss=find(del==j);
   tm_ssq=crossprod(eta.rows(indsss)) ;
   tm_sum=sum(eta.rows(indsss),0).t();
   
   if( norm(sum_of_samples.col(j)- tm_sum,"inf")>.1){
   Rcpp::Rcout<<"j= "<<j<<" jj = "<<jj<<" nj_val1="<<nj_val(j)<< endl;
   Rcpp::Rcout<<"sum diff dest"<<endl;
   Rcpp::stop("");
   }
   
   if(norm(( tm_ssq-sum_of_squares.slice(j) ),"inf")>.1 ){
   Rcpp::Rcout<<"j= "<<j<<" jj = "<<jj<<" nj_val1="<<nj_val(j)<< endl;
   Rcpp::Rcout<<tm_ssq-sum_of_squares.slice(j)<<endl;
   Rcpp::Rcout<<cross_prod-crossprod(eta.row(jj))<<endl;
   Rcpp::Rcout<<"ssq diff dest"<<endl;
   Rcpp::stop("");
   }*/
  ///////////////////////////////////////////////////////
}  


/*arma::umat sanity_check(arma::uvec &inds,  field<arma::uvec> &non_na_obs, arma::umat &eta, arma::uvec del){
 for(auto j:inds){
 // inds_eq_j[j] =find(del==j);
 inds_eq_j1[j] =find(del.head(n1)==j); nj_val1(j)=(inds_eq_j1[j]) .n_elem; 
 inds_eq_j2[j] =find(del.tail(n2)==j)+n1; nj_val2(j)=(inds_eq_j2[j]) .n_elem; 
 inds_eq_j[j] =join_cols(inds_eq_j1[j],inds_eq_j2[j]);
 nj_val(j)=nj_val1(j)+nj_val2(j);
 
 if(nj_val(j)){
 Rcpp::Rcout<<"n_inds at j="<<j<<" is "<<(inds_eq_j[j]) .n_elem<<endl;
 
 for(auto it1:inds_eq_j[j])//it1 iterates through cluster membership indicators
 for(auto it2:non_na_obs(it1)){//it2 iterates through non-missing variables for each observation
 // Rcpp::Rcout<<"it1= "<<it1<<" it2 ="<<it2<<" eta(it1,it2)= "<<eta(it1,it2)<<endl;
 ++noccu(j,it2)(eta(it1,it2));
 }
 }
 }
}*/

// [[Rcpp::export]]
void sanity(arma::uvec x){
  unsigned z=sum(x);
  Rcpp::Rcout<<"z= "<<z<<endl;
}

// [[Rcpp::export]]
Rcpp::List mixture_cat(const double alpha, const unsigned nmix, arma::uvec ncat,
                       int nrun, int burn, int thin, 
                       arma::umat eta1, arma::umat eta2,
                       Rcpp::List non_na_obs1,// Rcpp::List non_na_obs1, Rcpp::List non_na_obs2,
                       arma::uvec del1,arma::uvec del2){
  int acceptance=0,mh_count=0;
  const unsigned n1=eta1.n_rows,n2=eta2.n_rows, k=eta1.n_cols;
  const unsigned n=n1+n2;
  double t1ppp;
  
  arma::umat eta=join_cols( eta1, eta2 );
  arma::uvec del=join_cols(del1,del2);
  // Rcpp::Rcout<<"size(del)"<<size(del)<<endl;
  field<arma::uvec> non_na_obs(n);
  for(unsigned i=0;i<n;++i)
    non_na_obs(i)=Rcpp::as<arma::uvec>(non_na_obs1[i]);
  
  unsigned j;
  
  arma::umat alloc_var_mat(std::floor((nrun+burn)/thin), n); 
  arma::mat weights(std::floor((nrun+burn)/thin), n2), weights2(std::floor((nrun+burn)/thin), n2);
  arma::mat pimat1(std::floor((nrun+burn)/thin), nmix), pimat2(std::floor((nrun+burn)/thin), nmix);
  // --- initialise loop objects --- //
  ivec d(nmix); ////multinomial indicator for each sample
  arma::vec probs(nmix),log_probs(nmix), alpha_vec(k,fill::ones); ////assignment probability for each sample
  
  field<arma::uvec> noccu(nmix,k);
  for(unsigned i=0;i<nmix;++i){
    for(unsigned j=0;j<k;++j){
      if(ncat(j)<2){
        Rcpp::Rcout<<"j= "<<j<<" ncat(j)<2"<<endl;
        Rcpp::stop("");
      }
      noccu(i,j).set_size(ncat(j));
      noccu(i,j).zeros();
    }
  }
  
  bool thincheck, printcheck;
  /* arma::uvec *inds_eq_j, *inds_eq_j1,*inds_eq_j2;
   inds_eq_j=new arma::uvec[nmix ]; inds_eq_j1=new arma::uvec[nmix ]; inds_eq_j2=new arma::uvec[nmix ];*/
  field<arma::uvec> inds_eq_j(nmix), inds_eq_j1(nmix), inds_eq_j2(nmix);
  
  ///// Setting-up GSL random number generator for sampling from dirichlet
  const gsl_rng_type * T;
  gsl_rng * r;
  /* create a generator chosen by the
   environment variable GSL_RNG_TYPE */
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  gsl_rng_set(r,500);
  ///////////////////////////////////////////////////////////////////
  double dir_prec=alpha/nmix, prob_empty=log(dir_prec);
  ///////initialize parameters corresponding to clusters//////
  arma::uvec nj_val1(nmix),nj_val2(nmix), nj_val(nmix);  //occupancy number corresponding to each cluster
  // arma::vec log_nj_val1(nmix), log_nj_val2(nmix);
  
  Rcpp::Rcout<<"cluster occupancy finding starts!!"<<endl;
  for(unsigned j=0;j<nmix;++j){
    // inds_eq_j[j] =find(del==j);
    inds_eq_j1[j] =find(del.head(n1)==j); nj_val1(j)=(inds_eq_j1[j]) .n_elem; 
    inds_eq_j2[j] =find(del.tail(n2)==j)+n1; nj_val2(j)=(inds_eq_j2[j]) .n_elem; 
    inds_eq_j[j] =join_cols(inds_eq_j1[j],inds_eq_j2[j]);
    nj_val(j)=nj_val1(j)+nj_val2(j);
    
    if(nj_val(j)){
      Rcpp::Rcout<<"n_inds at j="<<j<<" is "<<(inds_eq_j[j]) .n_elem<<endl;
      /*if(nj_val2(j))
       log_nj_val2(j)=log(nj_val2(j)+dir_prec);*/
      
      for(auto it1:inds_eq_j[j])//it1 iterates through cluster membership indicators
        for(auto it2:non_na_obs(it1)){//it2 iterates through non-missing variables for each observation
          // Rcpp::Rcout<<"it1= "<<it1<<" it2 ="<<it2<<" eta(it1,it2)= "<<eta(it1,it2)<<endl;
          ++noccu(j,it2)(eta(it1,it2));
        }
    }
  }
  
  ///find the number of non-NA observations for each variable
  arma::umat nobs(nmix,k,fill::zeros);
  for(unsigned l=0;l<nmix;++l)
    for(unsigned j=0;j<k;++j){
      nobs(l,j)= sum(noccu(l,j));
      Rcpp::Rcout<<"nobs(j,l)= "<<nobs(l,j)<< " nj_val(j)= "<<nj_val(l)<<endl;
    }
    
    //////////////////////////
    //////////////////////////////////////////////////////////
    
    //----------------------------------------------//
    Rcpp::Rcout<<"loop starts"<<endl;
  // --- loop --- //
  arma::mat post_mat;
  // for(int i=0; i<nrun; i++, start++){
  arma::uvec non_empty_clusters1=(find(nj_val1)), non_empty_clusters2=(find(nj_val2));
  unsigned nmix1=non_empty_clusters2.n_elem;
  double dir_prec1=alpha/nmix1;
  
  for(int i=0; i<nrun+burn; i++){
    /////////////Sanity check /////////////
    /*Rcpp::Rcout<<"flag sanity"<<endl;
     arma::uvec nobs_check(k,fill::zeros);
     for(j=0;j<k;++j){
     for(unsigned l=0;l<nmix;++l)
     nobs_check(j)+= sum(noccu(l,j));
     if(nobs_check(j) != nobs(j)){
     Rcpp::Rcout<<"j = "<<j<<"nobs_check(j) != nobs(j)"<<endl;
     Rcpp::Rcout<<"nobs_check(j)= "<<nobs_check(j) <<"  nobs(j)= "<<nobs(j)<<endl;
     Rcpp::stop("");
     }
     }*/
    //////////////////////////
    
    // --- UPDATE ALLOCATION VARIABLES --- //
    arma::vec normal_exp(nmix),probs1;//exp2(nmix),exp1(nmix);
    double  log_probs_max;//,normal_exp; 
    long double log_DEN;
    // arma::mat resid(k,nmix);
    for(unsigned jj=0;jj<n;++jj){
      double dens,cluster_prob, pdf_empty;
      log_probs.fill(datum::log_min); //resetting log(probs) variable
      // Rcpp::Rcout<<"log_pdf_empty="<<pdf_empty<<endl;
      // Rcpp::stop("");
      // Rcpp::Rcout<<"jj="<<jj<<"\t"<<"del(jj)="<<del(jj)<<endl;
      unsigned current_ind=del(jj); //find the current index of data point jj
      if(jj<n1){//this is for G_1
        // Rcpp::Rcout<<" G1"<<endl;
        for(j=0;j<nmix;++j){
          // Rcpp::Rcout<<"j= "<<j<<" nj_val2(j)="<<nj_val2(j)<<endl;
          if(nj_val2(j)){
            if(j!=current_ind){
              dens=0.0;
              for(auto it:non_na_obs(jj))
                dens+=( log( noccu(j, it  )(eta(jj,it) )+  alpha_vec(it))-log(nobs(j,it) + ncat(it)*alpha_vec(it)  ) );
              
              // Rcpp::Rcout<<"at j="<<j<<" log_density="<<dens<<"\t clust prob="<<cluster_prob<<endl;
              cluster_prob=log(nj_val1(j)+dir_prec1);  //log_nj_val2(j);
              
              log_probs(j)=  dens +cluster_prob;
              // Rcpp::Rcout<<"log_probs "<<j<< "\t"<<log_probs(j)<<endl;
            }
          }
        }
        //calculating allocation probabilities
        j=current_ind;
        if(nj_val(j)==0 || nj_val2(j)==0)   Rcpp::stop("nj_val(current_ind)=0 || nj_val2(current_ind)=0 jj in G_1");
        // Rcpp::Rcout<<"flag 1"<<endl;
        dens=0.0;
        for(auto it:non_na_obs(jj))
          dens+=( log( noccu(j, it  )(eta(jj,it) )-1 +  alpha_vec(it))-log( ncat(it)*alpha_vec(it) + nobs(j,it) -1) );
        
        // Rcpp::Rcout<<"at j="<<j<<" log_density="<<dens<<"\t clust prob="<<cluster_prob<<endl;
        cluster_prob=log(nj_val1(j)-1+dir_prec1);  //log_nj_val2(j);
        // Rcpp::Rcout<<"flag 2"<<endl;
        log_probs(j)=  dens +cluster_prob;
        // (log_probs.t()).print("log_probs:");
        
        // non_empty_clusters2=  find(nj_val2);
        arma::vec log_probs1=log_probs(non_empty_clusters2);
        log_DEN=log_sum_exp(log_probs1);
        probs1=normalise(exp(log_probs1-log_DEN) ,1);
        // Rcpp::Rcout<<"flag 3"<<endl;
      } else{//this is for G_2
        // Rcpp::Rcout<<" G2"<<endl;
        j=current_ind;
        if(nj_val1(j) && nj_val2(j)==1){//the degenerate case
          log_probs.fill(datum::log_min);
          log_probs(j)=datum::log_max;
        } else if(nj_val1(j) && !nj_val2(j)  ){//the infeasible case
          Rcpp::stop("nj_val(current_ind)=0 || nj_val2(current_ind)=0 jj in G_2");
        } else {//other possible cases
          double pdf_empty= 0.0;
          pdf_empty=-sum(log(ncat(non_na_obs(jj))) );
          
          for(j=0;j<nmix;++j){
            if(j!=current_ind){
              if(!nj_val(j) ){
                dens=pdf_empty;
                cluster_prob=prob_empty;
              } else{
                dens=0.0;
                for(auto it:non_na_obs(jj))
                  dens+=( log( noccu(j, it  )(eta(jj,it) ) +  alpha_vec(it))-log( ncat(it)*alpha_vec(it) + nobs(j,it) ) );
                // Rcpp::Rcout<<"at j="<<j<<" log_density="<<dens<<"\t clust prob="<<cluster_prob<<endl;
                cluster_prob=log(nj_val2(j)+dir_prec);  //log_nj_val2(j);
              }
              log_probs(j)=  dens +cluster_prob;
            }
          }
          
          if(nj_val(current_ind)==1){
            dens= pdf_empty;
            cluster_prob=prob_empty;
          } else{ 
            j=current_ind;
            // Rcpp::Rcout<<"jj= "<<jj<<" del(jj)= "<<del(jj)<<" current_ind="<<current_ind<<endl;
            dens=0.0;
            for(auto it:non_na_obs(jj))
              dens+=( log( noccu(j, it  )(eta(jj,it) )-1 +  alpha_vec(it))-log( ncat(it)*alpha_vec(it) + nobs(j,it) -1) );
            cluster_prob=log(nj_val2(j)-1+dir_prec);
          }
          
          log_probs(current_ind)=  dens +cluster_prob;
        }
      }
      
      if(jj<n1){
        // Rcpp::Rcout<<"flag 3 0"<<endl;
        probs=zeros(nmix);
        probs(non_empty_clusters2)=probs1;
        // Rcpp::Rcout<<"flag 3 1"<<endl;
        // (probs1.t()).print("probs1 :");
        // (probs.t()).print("probs :");
      } else{
        log_DEN=log_sum_exp(log_probs);
        probs= normalise(exp(log_probs-log_DEN) ,1);
      }
      
      // log_probs.print("log_probs: ");
      // probs.print("probs: ");
      
      /*if(   gsl_fcmp(sum(probs),1.0,1e-5) ){
       // Rcpp::Rcout<< "At jj="<<jj<<"sum_prob is 0"<<endl;
       log_probs_max= max(log_probs);
       log_probs-=log_probs_max;
       
       probs=normalise(exp(log_probs) ,1);
      }*/
      
      ////check if sum of the allocation probbilities is zero
      if(sum(probs)==0) Rcpp::stop("sum(probs)=0");
      ///////////////////////////////////////////////////////
      // Rcpp::Rcout<<"flag 3 2"<<endl;
      // Rcpp::Rcout<<size(probs)<<"\t"<<size(d)<<"\t"<<size(del)<< "\t nmix= "<<nmix<<endl;
      // Rcpp::Rcout<<"jj= "<<jj<<endl;
      // gsl_ran_multinomial (r, nmix, 1, probs.begin(), d.begin());
      R::rmultinom(1, probs.begin(), nmix, d.begin());
      // (d.t()).print("d :");
      arma::uvec dd=find(d==1,1,"first");
      // (dd.t()).print("dd :");
      del(jj)=dd(  0);
      // Rcpp::Rcout<<"flag jj 4"<<endl;
      
      //updating cluster occupancies
      if(del(jj)!=current_ind){
        if(jj<n1){
          --nj_val1(current_ind);
          ++nj_val1(del(jj));
        } else{
          --nj_val2(current_ind);
          ++nj_val2(del(jj));
        }
        --nj_val(current_ind);
        ++nj_val(del(jj));
        // Rcpp::Rcout<<"flag 3 3"<<endl;
        
        for(auto it:non_na_obs(jj)){
          --noccu(current_ind, it  )(eta(jj,it) );
          ++noccu(del(jj), it  )(eta(jj,it) );
          
          --nobs(current_ind, it);
          ++nobs(del(jj), it);
          
          /*if(nobs(del(jj), it)>nj_val(del(jj)))
           Rcpp::stop("nobs(del(jj), it)>nj_val(del(jj))");*/
        }
        
        if(nj_val(del(jj)) != (nj_val2(del(jj)) + nj_val1(del(jj))) )
          Rcpp::stop("nj_val(del(jj)) != (nj_val2(del(jj)) + nj_val1(del(jj)))");
      }
    }
    // Rcpp::stop("1 MC iteration done!");
    non_empty_clusters1=(find(nj_val1));    non_empty_clusters2=(find(nj_val2));
    unsigned nmix1=non_empty_clusters2.n_elem;
    double dir_prec1=alpha/nmix1;
    Rcpp::Rcout<<"allocation variables updated"<<endl;
    
    // --- UPDATE MIXTURE PROBABILITY --- //
    unsigned counttttt=0,count1=0,count2=0;
    arma::vec pi1(nmix,fill::zeros), pi2(nmix),pi1_tmp(nmix1);
    for(j=0;j<nmix;++j){
      if(nj_val(j))
        ++counttttt;
      if(nj_val1(j)){
        ++count1;
        inds_eq_j1[j] =find(del.head(n1)==j);
      } else inds_eq_j1[j].reset();
      if(nj_val2(j)){
        ++count2;
        inds_eq_j2[j] =find(del.tail(n2)==j)+n1;
      } else inds_eq_j2[j].reset();
      inds_eq_j[j] =join_cols(inds_eq_j1[j],inds_eq_j2[j]);
      
      if((inds_eq_j1[j]).n_elem != nj_val1(j) || (inds_eq_j2[j]).n_elem != nj_val2(j) ){
        Rcpp::Rcout<<"inds_eq_j1[j]).n_elem="<<(inds_eq_j1[j]).n_elem<<"\t"<<"nj_val1("<<j<<")="<<nj_val1(j)<<endl;
        Rcpp::Rcout<<"inds_eq_j2[j]).n_elem="<<(inds_eq_j2[j]).n_elem<<"\t"<<"nj_val2("<<j<<")="<<nj_val2(j)<<endl;
        Rcpp::stop("Occupancy mismatch!!!");
      }
    }
    
    arma::vec dir_alpha=conv_to<arma::vec>::from( nj_val2)+ dir_prec;
    // gsl_ran_dirichlet(r, nmix, dir_alpha.begin(), pi2.begin());
    pi2=normalise(dir_alpha,1);
    
    arma::vec dir_alpha1=conv_to<arma::vec>::from( nj_val1(non_empty_clusters1))+ dir_prec1;
    /*arma::vec dir_alpha1=conv_to<arma::vec>::from( nj_val1(non_empty_clusters2))+ dir_prec1;
     gsl_ran_dirichlet(r, nmix1, dir_alpha1.begin(), pi1_tmp.begin());
     pi1(non_empty_clusters2)=pi1_tmp;
     pi1.zeros(); pi1(non_empty_clusters1)=normalise(dir_alpha1,1);
     */
    pi1.zeros(); pi1(non_empty_clusters1)=normalise(dir_alpha1,1); //assigning weights to only non-empty clusters in the posterior
    
    arma::vec wght(nmix,fill::zeros),wght2(nmix,fill::zeros); 
    ///assigning weights only to the non-empty clusters in group 1
    wght.zeros(); wght2.zeros();
    wght(non_empty_clusters1)=pi1(non_empty_clusters1)*(n2-1+alpha)/(dir_alpha(non_empty_clusters1)-1) ;   ///ratio weight
    wght2(non_empty_clusters1)=pi1(non_empty_clusters1)/nj_val2(non_empty_clusters1); //\pi_{1h}/n_{2h}
    
    
    arma::vec wght_x2=wght(del.tail(n2)), wght_x2_new=wght2(del.tail(n2));
    
    Rcpp::Rcout<< "# total clusters= "<< counttttt<<endl;
    Rcpp::Rcout<< "# total clusters in clinical arm= "<< count1<<endl;
    Rcpp::Rcout<< "# total clusters in synthetic arm= "<< count2<<endl;
    
    thincheck = i - std::floor(i/thin) * thin; // % operator stolen by arma
    if(!thincheck ) {
      Rcpp::Rcout<<"Iteration: "<<i<<endl;
    }
    
    
    int remainder= (i+1 );
    int quotient= (int) std::floor(remainder/thin);
    remainder-= (quotient*thin) ;
    
    if(remainder==0){
      alloc_var_mat.row(quotient-1)=del.t();
      pimat1.row(quotient-1)=pi1.t();
      pimat2.row(quotient-1)=pi2.t();
      weights.row(quotient-1)=wght_x2.t(); //Ratio of cluster probabilities
      weights2.row(quotient-1)=wght_x2_new.t(); //probability of the cluster
    }
  }
  
  // delete[] inds_eq_j;   delete[] inds_eq_j1;   delete[] inds_eq_j2;
  gsl_rng_free (r);
  
  
  return Rcpp::List::create(Rcpp::Named("pimat1") =pimat1,
                            Rcpp::Named("pimat2") =pimat2,
                            Rcpp::Named("Allocation variables") = alloc_var_mat,
                            Rcpp::Named("Weights")=weights,
                            Rcpp::Named("Weights2")=weights2);
}


/*******functions required for HMC on Dirichlet mixture parameters**********/
double ll_alpha(const double l_alpha, const unsigned K, const arma::vec &nj_vec,
                const double mu_alp, const double sig_alp){
  const double alpha=exp(l_alpha);
  const unsigned Kn=nj_vec.n_elem;
  if(K<Kn)
    Rcpp::stop("Error in ll_alpha: # Mixtures < # Occupied clusters! ");
  
  const double al_by_k= alpha/K;
  
  double ret= lgamma(alpha) - lgamma(alpha+ accu(nj_vec)) + accu(lgamma(nj_vec+al_by_k))
    - Kn * lgamma(al_by_k) + log_normpdf(l_alpha , mu_alp, sig_alp );
  
  return ret;
}

double del_ll_alpha(const double l_alpha, const unsigned K, const arma::vec &nj_vec,
                    const double mu_alp, const double sig_alp){
  const double alpha=exp(l_alpha);
  const unsigned Kn=nj_vec.n_elem;
  if(K<Kn)
    Rcpp::stop("Error in ll_alpha: # Mixtures < # Occupied clusters! ");
  
  const double al_by_k= alpha/((double) K);
  
  vec tmp=nj_vec+ al_by_k;
  tmp.transform( [](double val) { return ( gsl_sf_psi(val) ); } );
  
  // Rcpp::Rcout<<"alpha="<<alpha<<" al_by_k="<<al_by_k<<" alpha+ accu(nj_vec)="<<alpha+ accu(nj_vec) <<endl;
  double ret= (gsl_sf_psi(alpha) - gsl_sf_psi(alpha+ accu(nj_vec)) + (accu(tmp) -  Kn*gsl_sf_psi(al_by_k) )/ ((double) K) )*alpha 
  - (l_alpha - mu_alp)/gsl_pow_2(sig_alp) ;
  
  return ret;
}

void leapfrog_dir_alpha(const unsigned nstep,const double delta, double &v_old,  double &p_lam, double &l_alpha, 
                        const unsigned K, const arma::vec &nj_vec, const double mu_alp, const double sig_alp){
  for(unsigned i=0;i<nstep;++i){
    l_alpha+=(delta)*(p_lam -(delta/2)*v_old);
    l_alpha = std::clamp(l_alpha, -1e1, 1e1);
    
    double v_new=-del_ll_alpha(l_alpha, K, nj_vec, mu_alp, sig_alp);
    p_lam-=(delta/2)* ( v_old+v_new);
    v_old=v_new;
  }
}

void update_alpha(const unsigned K, const arma::vec &nj_vec, const double mu_alp, const double sig_alp, 
                  const arma::vec &del_range_alp, const unsigned nleapfrog_alp, double &l_alpha, unsigned &acceptance){
  /*double mu_alp, sig_alp;
   sig_alp=log1p(ps_hyper(1) /gsl_pow_2(ps_hyper(0)));mu_alp=log(ps_hyper(0))-sig_alp/2;sig_alp=sqrt(sig_alp);*/
  
  double ll_old=-ll_alpha(l_alpha, K, nj_vec, mu_alp, sig_alp), v_old=-del_ll_alpha(l_alpha, K, nj_vec, mu_alp, sig_alp);
  
  double p_lam=randn();
  double kin_energy= gsl_pow_2(p_lam);
  
  double l_alpha_new=l_alpha, v_new=v_old;
  
  unsigned pois_draw=(unsigned) R::rpois(nleapfrog_alp); //randomly generating no. of leapfrog steps
  unsigned nstep=GSL_MAX_INT(1,pois_draw); 
  // nstep=GSL_MIN_INT(nstep,leapmax);
  double delta= R::runif(del_range_alp(0),del_range_alp(1) );//randomly generating \delta t 
  // Rcpp::Rcout<<"Flag X params -.5!!"<<endl;
  leapfrog_dir_alpha(nstep, delta, v_new,  p_lam, l_alpha_new, K, nj_vec, mu_alp, sig_alp);
  // params_new.print("params_new");
  // Rcpp::Rcout<<"Flag X params 0!!"<<endl;
  double ll_new=- ll_alpha(l_alpha_new, K, nj_vec, mu_alp, sig_alp); //log-likelihood at the new proposed value
  
  //get H_ll_new
  double H_new= ll_new+  gsl_pow_2(p_lam) /2, H_old=ll_old+kin_energy/2;
  
  if(log(randu())< -(H_new-H_old) ){
    l_alpha=l_alpha_new;
    /*ll_old=ll_new;
     v_old=v_new;*/
    ++acceptance;
  }
  /*Rcpp::Rcout<<"i= "<<i<< "Acceptance rate="<<(((double)acceptance)/ ((double) i))<<" alpha= "<<exp(l_alpha)<<endl;
   Rcpp::Rcout<<"Acceptance rate="<<(acceptance/n_mc)<<endl;
   return exp(alpha_vec);*/
}
/***************************************************/


/*******functions required for MALA on exponential parameters**********/
inline double log_gamma_dens(double x, double a, double b){
  //exp(x)~ Ga(a,b) with mean a/b and variance a/b^2
  return a*x -b*exp(x);
}

inline double del_log_gamma_dens(double x, double a, double b){
  //exp(x)~ Ga(a,b) with mean a/b and variance a/b^2
  return a -b*exp(x);
}


double logpi(const arma::uvec &nj_val1, const arma::uvec &nj_val2, 
             const arma::uvec &isfailure1, const arma::uvec &isfailure2,
             const arma::vec &survtime1, const arma::vec &survtime2,
             const arma::vec &params,
             double aa,double ba, double al, double bl){
  double a0=exp(params(0)),  l0=exp(params(1));
  double sum= log_gamma_dens( params(0), aa, ba) + log_gamma_dens( params(1), al, bl);
  
  // gsl_sf_lnpoch(a,b)=log (gamma(a+b)/gamma(a)  )
  double a_times_log_l= a0* params(1);
  for(unsigned j=0;j< nj_val2.n_elem;++j){
    if(nj_val2(j)){
      sum+= (gsl_sf_lnpoch (a0, isfailure2(j) )  - (a0+isfailure2(j))*log(l0+survtime2(j) ) +a_times_log_l );
      if(nj_val1(j))
        sum+= (gsl_sf_lnpoch (a0, isfailure1(j) )  - (a0+isfailure1(j))*log(l0+survtime1(j) ) + a_times_log_l);
    }
  }
  return sum;
}

arma::vec delpi(const arma::uvec &nj_val1, const arma::uvec &nj_val2, 
                const arma::uvec &isfailure1, const arma::uvec &isfailure2,
                const arma::vec &survtime1, const arma::vec &survtime2,
                const arma::vec &params,
                double aa,double ba, double al, double bl){
  double  a0=exp(params(0)),  l0=exp(params(1)), dela= 0.0, dell=0.0;
  
  // if(failcheck){
  //   (params.t()).print("params: ");
  //   Rcpp::Rcout<<"a0= "<<a0<<" l0= "<<l0<<endl;
  //   Rcpp::Rcout<<"a0/l0= "<<a0/l0<<" , "<< exp( params(0)- params(1)) <<endl;
  // } 
  for(unsigned j=0;j< nj_val2.n_elem;++j){
    // gsl_sf_psi is digamma function =d/dx (log(\Gamma(x)))
    if(nj_val2(j)){
      dela+=( gsl_sf_psi(a0 + isfailure2(j) ) - gsl_sf_psi(a0) + params(1) - log(l0+survtime2(j) ) );
      dell+=( - ((a0 + isfailure2(j))/(l0+survtime2(j) )) *l0 +a0);
      // if(failcheck)       Rcpp::Rcout<< ( -(a0 + isfailure2(j))/(l0+survtime2(j) ))<<endl;
      
      if(nj_val1(j)){
        dela+= (gsl_sf_psi(a0 + isfailure1(j) ) - gsl_sf_psi(a0) + params(1) - log(l0+survtime1(j) ));
        dell+=( -((a0 + isfailure1(j))/(l0+survtime1(j) )) *l0 +a0 );
        // if(failcheck)       Rcpp::Rcout<< ( -(a0 + isfailure1(j))/(l0+survtime1(j) ))<<endl;
      }
    }
  }
  
  /*if(failcheck)
   Rcpp::Rcout<<"del_a= "<<dela*a0<< " del_l= "<<dell<<endl;*/
  
  arma::vec del_pi(2);
  del_pi(0)=dela* a0 + del_log_gamma_dens(params(0), aa,ba);
  del_pi(1)=dell+ del_log_gamma_dens(params(1), al,bl);
  return del_pi;
}
/****************************************************************/

double post_t_dens(const double x, const  double ss_survtime,const  double survtime,
                   const double df0,const double a0,const double mu0,const double beta0,const unsigned n){
  ///ss_survtime is ssq and survtime is just the sum
  double df_post=df0+ n , alpha_post=a0+ n/2.0;
  double tmp=(   n *df0 )/df_post, ss_j, survtime_j, mean_j, mu_post, beta_post, sigma_post;
  
  if(n){
    survtime_j = survtime ;
    mean_j= survtime_j/ n  ;
    ss_j= ss_survtime -n *gsl_pow_2(mean_j);
  } else ss_j=  survtime_j = mean_j= 0.0;
  
  
  mu_post = (df0*mu0 + survtime_j) / df_post;
  
  beta_post=beta0+  (ss_j   + tmp* gsl_pow_2(mean_j - mu0) )/2.0; //ss_j differently defined from the rest of the code
  sigma_post =sqrt( beta_post * (df_post+1)/(df_post *alpha_post) );
  
  
  // gsl_sf_lnpoch(a,b)=log (gamma(a+b)/gamma(a)  )
  
  /*double df_final=2.0*alpha_post;
   double denss=gsl_sf_lnpoch(alpha_post,.5) - (M_LNPI/2 + log( sqrt(df_final)* sigma_post ) +(alpha_post+.5)* log1p( gsl_pow_2( (x-mu_post)/sigma_post )/df_final  ) );
   Rcpp::Rcout<<"mu_post="<<mu_post<<" ss_j="<<ss_j <<  " manual dens="<<denss<<" RcppDIST dens="<<d_lst( x,  2.0*alpha_post, mu_post, sigma_post, 1)<<endl;*/
  
  return d_lst( x,  2.0*alpha_post, mu_post, sigma_post, 1);
}

inline double surv_fn_lognorm(const double x, const unsigned nu,const  double ss_survtime,const  double survtime,
                              const double df0,const double a0,const double mu0,const double beta0,const unsigned n){//nu=1 implies failure time
  ///don't need x_aug argument here!!
  double df_post=df0+ n , alpha_post=a0+ n/2.0;
  double tmp=(   n *df0 )/df_post, ss_j, survtime_j, mean_j, mu_post, beta_post, sigma_post;
  
  if(n){
    survtime_j = survtime ;
    mean_j= survtime_j/ n  ;
    ss_j= ss_survtime -n *gsl_pow_2(mean_j);
  } else ss_j=  survtime_j = mean_j= 0.0;
  
  
  mu_post = (df0*mu0 + survtime_j) / df_post;
  beta_post=beta0+  (ss_j   + tmp* gsl_pow_2(mean_j - mu0) )/2.0; //ss_j differently defined from the rest of the code
  sigma_post =sqrt( beta_post * (df_post+1)/(df_post *alpha_post) );
  
  /*Rcpp::Rcout<<"mu0= "<<mu0 << " beta0= "<<beta0<<endl;
   Rcpp::Rcout<<"mean_j="<<mean_j<<" survtime_j="<<survtime_j;
   Rcpp::Rcout<< " ss_j=" << ss_j<< " gsl_pow_2(mean_j - mu0)=" <<gsl_pow_2(mean_j - mu0)<<" n="<<n << " tmp="<<tmp<< " df_post="<<df_post<< " alpha_post="<<alpha_post<< 
   " beta_post="<<beta_post<<" sigma_post="<<sigma_post<<endl;
   Rcpp::Rcout<<"survtime_j= "<<survtime_j<<" ss_j= "<<ss_j<<  " df = "<< 2*df_post<< " mu_post = "<< mu_post <<" sigma_post = "<<sigma_post<< endl;
   Rcpp::Rcout<<"x= "<<x<<" x_aug= "<< x_aug<<endl;*/
  
  return nu ?  ( d_lst( x,  2.0*alpha_post, mu_post, sigma_post, 1) ) : ( p_lst(x, 2.0*alpha_post, mu_post, sigma_post, 0, 1 ) );
}

arma::vec sim_lognorm_params( const  double ss_j,const  double survtime_j,
                              const double df0,const double a0,const double mu0,const double beta0,const unsigned n){
  double df_post=df0+ n, alpha_post=a0+ n/2.0;
  double tmp=( n*df0 )/df_post;
  double mean_j=  n ? (survtime_j/ n) : 0.0;
  double mu_post = (df0*mu0 + survtime_j) / df_post;
  double beta_post=beta0+  (ss_j -  ((double) n)  * gsl_pow_2(mean_j )  + tmp* gsl_pow_2(mean_j - mu0) )/2.0;
  double sigma_post =sqrt( beta_post  / (alpha_post * df_post ) );
  
  arma::vec ret(2);
  // Rcpp::Rcout<<"Flag lognormal_param 0 alpha_post= "<<alpha_post<<" beta_post= "<< beta_post<<" data beta= "<<(ss_j -  ((double) n)  * gsl_pow_2(mean_j )  + tmp* gsl_pow_2(mean_j - mu0) )/2.0<<endl;
  ret(0)=r_lst(2*alpha_post,  mu_post, sigma_post );
  ret(1)= 1/ randg( distr_param(alpha_post, 1/beta_post) );
  return ret;
}

double logpi_lognorm(const arma::vec &nj_val1, const arma::vec &nj_val2, 
                     const arma::vec &survtime1, const arma::vec &ss_survtime1,
                     const arma::vec &survtime2, const arma::vec &ss_survtime2,
                     const arma::vec &params,
                     const double a0, const double df0, 
                     const double mu_m, const double mu_v,
                     const double b_m, const double b_v){
  //mu_v and b_v ae prior variances, not sds!!!
  double mu0=params(0),  beta0=exp(params(1)), log_b0=params(1);
  double sum= - (gsl_pow_2(mu0-mu_m)/mu_v + gsl_pow_2( log_b0 - b_m) / b_v) /2; //normal prior
  //log_normpdf(mu0,mu_m, sqrt(mu_v)) + log_normpdf(log_b0 , b_m, sqrt(b_v));
  
  
  // Rcpp::Rcout<<"sum= "<<sum<<endl;
  for(unsigned j=0; j< nj_val2.n_elem; ++j){//Iterate over the clusters
    ///notations mostly follow Wikipedia conjugate prior NIG
    if(nj_val2(j)){
      double df_post=df0+ nj_val2(j), alpha_post=a0+ nj_val2(j)/2.0; 
      double tmp=( nj_val2(j)*df0 )/df_post;
      double ss_j= ss_survtime2(j) , survtime_j = survtime2(j) ;
      double mean_j= survtime_j/nj_val2(j) ;
      double beta_post=beta0+  (ss_j -  nj_val2(j) * gsl_pow_2(mean_j )  + tmp* gsl_pow_2(mean_j - mu0) )/2.0;
      sum+= (a0 *log_b0 -alpha_post *  log( beta_post));
      /*Rcpp::Rcout<<"df_post= "<<df_post<<" tmp="<<tmp<<" beta_post= "<<beta_post<<" ss_j= "<<ss_j<<" nj_val2(j) * gsl_pow_2(mean_j )="<<nj_val2(j) * gsl_pow_2(mean_j )<<  endl;
       Rcpp::Rcout<<"density at j= "<< j<<"is "<<(mu0 *log_b0 -alpha_post *  log( beta_post)) <<endl;*/
      
      if(nj_val1(j)){
        df_post=df0+ nj_val1(j); alpha_post=a0+ nj_val1(j) /2.0;
        tmp=( nj_val1(j) *df0 )/df_post;
        ss_j= ss_survtime1(j), survtime_j = survtime1(j);
        mean_j= survtime_j/nj_val1(j);
        beta_post=beta0+  (ss_j -  nj_val1(j) * gsl_pow_2(mean_j )  + tmp* gsl_pow_2(mean_j - mu0) )/2.0;
        sum+= (a0 *log_b0 -alpha_post *  log( beta_post));
      }
    }
  }
  return sum;
}

arma::vec delpi_lognorm(const arma::vec &nj_val1, const arma::vec &nj_val2, 
                        const arma::vec &survtime1, const arma::vec &ss_survtime1,
                        const arma::vec &survtime2, const arma::vec &ss_survtime2,
                        const arma::vec &params,
                        double a0, double df0, 
                        double mu_m, double mu_v,
                        double b_m, double b_v){
  //mu_v and b_v are prior variances, not sds!!!
  double mu0=params(0),  beta0=exp(params(1));
  double del_m=-(mu0-mu_m)/mu_v, del_b= -( params(1) - b_m) / b_v;//derivative of normal prior
  
  for(unsigned j=0; j<nj_val2.n_elem; ++j){//Iterate over the clusters
    ///notations mostly follow Wikipedia conjugate prior NIG
    ////UPDATE for the experimental arm
    if(nj_val2(j)){
      double df_post=df0+ nj_val2(j), alpha_post=a0+ nj_val2(j)/2.0;
      double tmp=( nj_val2(j)*df0 )/df_post;
      double ss_j= ss_survtime2(j) , survtime_j = survtime2(j) ;
      double mean_j= survtime_j/nj_val2(j) ;
      double beta_post=beta0+  (ss_j -   nj_val2(j) * gsl_pow_2(mean_j )  + tmp* gsl_pow_2(mean_j - mu0) )/2.0;
      double common_part=alpha_post / beta_post ;
      del_m+= (common_part * tmp*(mean_j-mu0));
      del_b+= (a0 -common_part*beta0) ;
      
      if(nj_val1(j)){
        df_post=df0+ nj_val1(j); alpha_post=a0+ nj_val1(j) /2.0;
        tmp=( nj_val1(j) *df0 )/df_post;
        ss_j= ss_survtime1(j), survtime_j = survtime1(j);
        mean_j= survtime_j/nj_val1(j);
        beta_post=beta0+  (ss_j -  nj_val1(j) * gsl_pow_2(mean_j )  + tmp* gsl_pow_2(mean_j - mu0) )/2.0;
        common_part=alpha_post / beta_post ;
        del_m+= (common_part * tmp*(mean_j-mu0));
        del_b+= (a0 -common_part*beta0) ;
      }
    }
  }
  arma::vec ret(2);
  ret(0)=del_m; ret(1)= del_b;
  return ret;
}


void leapfrog_lognorm_hyper(const unsigned nstep,const double delta, arma::vec &v_old,  arma::vec &p_lam, arma::vec &params, 
                            const arma::vec &nj_val1, const arma::vec &nj_val2, 
                            const arma::vec &survtime1, const arma::vec &ss_survtime1,
                            const arma::vec &survtime2, const arma::vec &ss_survtime2,
                            const double a0,const  double df0, 
                            const double mu_m, const  double mu_v,
                            const double b_m, const double b_v){
  for(unsigned i=0;i<nstep;++i){
    // Rcpp::Rcout<<"flag -1 leap i="<<i<<endl;
    params+=(delta)*(p_lam -(delta/2)*v_old);
    // Rcpp::Rcout<<"flag 0 leap i="<<i<<endl;
    params(0)=std::clamp(params(0), -1e3, 1e5);
    params(1)=std::clamp(params(1), -1e1, 1e1);
    
    arma::vec v_new=-delpi_lognorm(nj_val1, nj_val2, survtime1, ss_survtime1, survtime2, ss_survtime2,
                                   params, a0, df0, mu_m, mu_v, b_m, b_v);
    // Rcpp::Rcout<<"flag 1.5 leap i="<<i<<endl;
    p_lam-=(delta/2)* ( v_old+v_new);
    v_old=v_new;
  }
}

void update_lognorm_hyper(const arma::vec &nj_val1, const arma::vec &nj_val2, 
                          const arma::vec &survtime1, const arma::vec &ss_survtime1,
                          const arma::vec &survtime2, const arma::vec &ss_survtime2,
                          const double a0,const  double df0, 
                          const double mu_m, const  double mu_v,
                          const double b_m, const double b_v,
                          const arma::vec &del_range, const unsigned nleapfrog, arma::vec &params, unsigned &acceptance){
  double ll_old=-logpi_lognorm(nj_val1, nj_val2, survtime1, ss_survtime1, survtime2, ss_survtime2,
                               params, a0, df0, mu_m, mu_v, b_m, b_v);
  arma::vec v_old=-delpi_lognorm(nj_val1, nj_val2, survtime1, ss_survtime1, survtime2, ss_survtime2,
                                 params, a0, df0, mu_m, mu_v, b_m, b_v);
  
  arma::vec p_lam(2, fill::randn);
  double kin_energy= dot(p_lam, p_lam);
  
  arma::vec params_new=params, v_new=v_old;
  
  unsigned pois_draw=(unsigned) R::rpois(nleapfrog);
  unsigned nstep=GSL_MAX_INT(1,pois_draw);
  double delta= R::runif(del_range(0),del_range(1) );
  // Rcpp::Rcout<<"Flag X params -.5!!"<<endl;
  leapfrog_lognorm_hyper(nstep, delta, v_new, p_lam, params_new, 
                         nj_val1, nj_val2, survtime1, ss_survtime1, survtime2, ss_survtime2,
                         a0, df0, mu_m, mu_v, b_m, b_v);
  // params_new.print("params_new");
  // Rcpp::Rcout<<"Flag X params 0!!"<<endl;
  double ll_new=- logpi_lognorm(nj_val1, nj_val2, survtime1, ss_survtime1, survtime2, ss_survtime2,
                                params_new, a0, df0, mu_m, mu_v, b_m, b_v); //log-likelihood
  
  //get H_ll_new
  double H_new= ll_new+  ssq(p_lam) /2, H_old=ll_old+kin_energy/2;
  
  if(log(randu())< -(H_new-H_old) ){
    params=params_new;
    /*ll_old=ll_new;
     v_old=v_new;*/
    ++acceptance;
  }
  /*Rcpp::Rcout<<"i= "<<i<< "Acceptance rate="<<(((double)acceptance)/ ((double) i))<<" alpha= "<<exp(l_alpha)<<endl;
   Rcpp::Rcout<<"Acceptance rate="<<(acceptance/n_mc)<<endl;
   return exp(alpha_vec);*/
}
/*****************************************************************/

// [[Rcpp::export]]
Rcpp::List common_atoms_cat(const double alpha, const unsigned nmix, arma::uvec ncat,
                            const double aa,const double ba,const double al,const double bl, //gamma hyperparameters for response
                            int nrun, int burn, int thin, 
                            arma::umat eta1, arma::umat eta2,
                            arma::vec st1, arma::uvec nu1, //responses for the treatment
                            arma::vec st2, arma::uvec nu2, //responses for the synthetic control
                            Rcpp::List non_na_obs1,// Rcpp::List non_na_obs1, Rcpp::List non_na_obs2,
                            arma::uvec del1,arma::uvec del2, double tau){
  unsigned acceptance=0;
  const unsigned n1=eta1.n_rows,n2=eta2.n_rows, k=eta1.n_cols;
  const unsigned n=n1+n2;
  double t1ppp;
  
  double a0=2.0, l0=120.0;
  arma::vec current_params(2);  current_params(0)=log(a0); current_params(1)=log(l0);
  
  arma::umat eta=join_cols( eta1, eta2 );
  arma::uvec del=join_cols(del1,del2), nu=join_cols(nu1,nu2) ;
  arma::vec st=join_cols(st1,st2);
  // Rcpp::Rcout<<"size(del)"<<size(del)<<endl;
  field<arma::uvec> non_na_obs(n);
  for(unsigned i=0;i<n;++i)
    non_na_obs(i)=Rcpp::as<arma::uvec>(non_na_obs1[i]);
  
  unsigned j;
  
  
  /////Define MCMC storage matrices
  unsigned n_mc=std::floor((nrun+burn)/thin);
  arma::umat alloc_var_mat(n_mc, n); 
  arma::mat weights2(n_mc, n2);//, weights(n_mc, n2);
  arma::mat pimat1(n_mc, nmix), pimat2(n_mc, nmix);
  arma::mat thetamat1(n_mc, nmix), thetamat2(n_mc, nmix);//, HRmat(n_mc, nmix);
  arma::mat hyperparams(n_mc, 2);
  arma::mat unifmat(n_mc,n);
  ///////////////////////
  
  
  // --- initialise loop objects --- //
  ivec d(nmix); ////multinomial indicator for each sample
  arma::vec probs(nmix),log_probs(nmix), alpha_vec(k,fill::ones); ////assignment probability for each sample
  
  field<arma::uvec> noccu(nmix,k);
  for(unsigned i=0;i<nmix;++i){
    for(unsigned j=0;j<k;++j){
      if(ncat(j)<2){
        Rcpp::Rcout<<"j= "<<j<<" ncat(j)<2"<<endl;
        Rcpp::stop("");
      }
      noccu(i,j).set_size(ncat(j));
      noccu(i,j).zeros();
    }
  }
  
  bool thincheck, printcheck;
  /* arma::uvec *inds_eq_j, *inds_eq_j1,*inds_eq_j2;
   inds_eq_j=new arma::uvec[nmix ]; inds_eq_j1=new arma::uvec[nmix ]; inds_eq_j2=new arma::uvec[nmix ];*/
  field<arma::uvec> inds_eq_j(nmix), inds_eq_j1(nmix), inds_eq_j2(nmix);
  
  ///// Setting-up GSL random number generator for sampling from dirichlet
  const gsl_rng_type * T;
  gsl_rng * r;
  /* create a generator chosen by the
   environment variable GSL_RNG_TYPE */
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  gsl_rng_set(r,500);
  ///////////////////////////////////////////////////////////////////
  double dir_prec=alpha/nmix, prob_empty=log(dir_prec);
  ///////initialize parameters corresponding to clusters//////
  arma::uvec nj_val1(nmix),nj_val2(nmix), nj_val(nmix),//occupancy number corresponding to each cluster
  isfailure1(nmix,fill::zeros), isfailure2(nmix,fill::zeros);  //sum of failure time indicators in each cluster
  arma::vec survtime1(nmix,fill::zeros), survtime2(nmix,fill::zeros) ;
  // arma::vec log_nj_val1(nmix), log_nj_val2(nmix);
  
  Rcpp::Rcout<<"cluster occupancy finding starts!!"<<endl;
  for(unsigned j=0;j<nmix;++j){
    // inds_eq_j[j] =find(del==j);
    inds_eq_j1[j] =find(del.head(n1)==j); nj_val1(j)=(inds_eq_j1[j]) .n_elem; 
    inds_eq_j2[j] =find(del.tail(n2)==j)+n1; nj_val2(j)=(inds_eq_j2[j]) .n_elem; 
    inds_eq_j[j] =join_cols(inds_eq_j1[j],inds_eq_j2[j]);
    nj_val(j)=nj_val1(j)+nj_val2(j);
    
    survtime1(j)= sum(st(inds_eq_j1[j]) ); isfailure1(j)=sum(nu(inds_eq_j1[j] ));
    survtime2(j)= sum(st(inds_eq_j2[j]) ); isfailure2(j)=sum(nu(inds_eq_j2[j] ));
    
    if(nj_val(j)){
      Rcpp::Rcout<<"n_inds at j="<<j<<" is "<<(inds_eq_j[j]) .n_elem<<endl;
      /*if(nj_val2(j))
       log_nj_val2(j)=log(nj_val2(j)+dir_prec);*/
      
      for(auto it1:inds_eq_j[j])//it1 iterates through cluster membership indicators
        for(auto it2:non_na_obs(it1)){//it2 iterates through non-missing variables for each observation
          // Rcpp::Rcout<<"it1= "<<it1<<" it2 ="<<it2<<" eta(it1,it2)= "<<eta(it1,it2)<<endl;
          ++noccu(j,it2)(eta(it1,it2));
        }
    }
  }
  const arma::uvec censored_indices=find(nu==0);
  
  ///find the number of non-NA observations for each variable
  arma::umat nobs(nmix,k,fill::zeros);
  for(unsigned l=0;l<nmix;++l)
    for(unsigned j=0;j<k;++j){
      nobs(l,j)= sum(noccu(l,j));
      Rcpp::Rcout<<"nobs(j,l)= "<<nobs(l,j)<< " nj_val(j)= "<<nj_val(l)<<endl;
    }
    //////////////////////////////////////////////////////////
    
    //----------------------------------------------//
    Rcpp::Rcout<<"loop starts"<<endl;
  // --- loop --- //
  arma::uvec non_empty_clusters1=(find(nj_val1)), non_empty_clusters2=(find(nj_val2));
  unsigned nmix1=non_empty_clusters2.n_elem;
  double dir_prec1=alpha/nmix1;
  arma::vec unif1, unif2, unif;
  
  for(int i=0; i<nrun+burn; i++){
    /////////////Sanity check /////////////
    /*Rcpp::Rcout<<"flag sanity"<<endl;
     arma::uvec nobs_check(k,fill::zeros);
     for(j=0;j<k;++j){
     for(unsigned l=0;l<nmix;++l)
     nobs_check(j)+= sum(noccu(l,j));
     if(nobs_check(j) != nobs(j)){
     Rcpp::Rcout<<"j = "<<j<<"nobs_check(j) != nobs(j)"<<endl;
     Rcpp::Rcout<<"nobs_check(j)= "<<nobs_check(j) <<"  nobs(j)= "<<nobs(j)<<endl;
     Rcpp::stop("");
     }
     }*/
    //////////////////////////
    
    // --- UPDATE ALLOCATION VARIABLES --- //
    arma::vec normal_exp(nmix),probs1;//exp2(nmix),exp1(nmix);
    double  log_probs_max;//,normal_exp; 
    long double log_DEN;
    for(unsigned jj=0;jj<n;++jj){
      double dens,cluster_prob, pdf_empty;
      log_probs.fill(datum::log_min);
      // Rcpp::Rcout<<"log_pdf_empty="<<pdf_empty<<endl;
      // Rcpp::stop("");
      // Rcpp::Rcout<<"jj="<<jj<<"\t"<<"del(jj)="<<del(jj)<<endl;
      unsigned current_ind=del(jj); //find the current index of data point jj
      if(jj<n1){//this is for G_1
        // Rcpp::Rcout<<" G1"<<endl;
        for(j=0;j<nmix;++j){
          // Rcpp::Rcout<<"j= "<<j<<" nj_val2(j)="<<nj_val2(j)<<endl;
          if(nj_val2(j)){
            if(j!=current_ind){
              dens= surv_fn(st(jj),nu(jj),a0+isfailure1(j), l0+survtime1(j) ) ;
              //lgamma(a0+isfailure1(j)+ nu(jj) ) -(a0+isfailure1(j)+ nu(jj) )  *log(l0+survtime1(j)+ st(jj));
              for(auto it:non_na_obs(jj))
                dens+=( log( noccu(j, it  )(eta(jj,it) )+  alpha_vec(it))-log(nobs(j,it) + ncat(it)*alpha_vec(it)  ) );
              
              // Rcpp::Rcout<<"at j="<<j<<" log_density="<<dens<<"\t clust prob="<<cluster_prob<<endl;
              cluster_prob=log(nj_val1(j)+dir_prec1);  //log_nj_val2(j);
              
              log_probs(j)=  dens +cluster_prob;
              // Rcpp::Rcout<<"log_probs "<<j<< "\t"<<log_probs(j)<<endl;
            }
          }
        }
        //calculating allocation probabilities
        j=current_ind;
        if(nj_val(j)==0 || nj_val2(j)==0)   Rcpp::stop("nj_val(current_ind)=0 || nj_val2(current_ind)=0 jj in G_1");
        // Rcpp::Rcout<<"flag 1"<<endl;
        dens=surv_fn(st(jj),nu(jj),a0+isfailure1(j)-nu(jj), l0+survtime1(j) -st(jj)) ;
        //lgamma(a0+isfailure1(j) ) -(a0+isfailure1(j) )  *log(l0+survtime1(j) );
        for(auto it:non_na_obs(jj))
          dens+=( log( noccu(j, it  )(eta(jj,it) )-1 +  alpha_vec(it))-log( ncat(it)*alpha_vec(it) + nobs(j,it) -1) );
        
        // Rcpp::Rcout<<"at j="<<j<<" log_density="<<dens<<"\t clust prob="<<cluster_prob<<endl;
        cluster_prob=log(nj_val1(j)-1+dir_prec1);  //log_nj_val2(j);
        // Rcpp::Rcout<<"flag 2"<<endl;
        log_probs(j)=  dens +cluster_prob;
        // (log_probs.t()).print("log_probs:");
        
        // non_empty_clusters2=  find(nj_val2);
        arma::vec log_probs1=log_probs(non_empty_clusters2);
        log_DEN=log_sum_exp(log_probs1);
        probs1=normalise(exp(log_probs1-log_DEN) ,1);
        // Rcpp::Rcout<<"flag 3"<<endl;
      } else{//this is for G_2
        // Rcpp::Rcout<<" G2"<<endl;
        j=current_ind;
        if(nj_val1(j) && nj_val2(j)==1){//the degenerate case
          log_probs.fill(datum::log_min);
          log_probs(j)=datum::log_max;
        } else if(nj_val1(j) && !nj_val2(j)  ){//the infeasible case
          Rcpp::stop("nj_val(current_ind)=0 || nj_val2(current_ind)=0 jj in G_2");
        } else {//other possible cases
          double pdf_empty = surv_fn(st(jj),nu(jj),a0, l0  ) ;
          //-sum(log(ncat(non_na_obs(jj))) )+ (lgamma(a0+ nu(jj) ) -(a0+ nu(jj) )  *log(l0+ st(jj)) ) ;
          
          for(j=0;j<nmix;++j){
            if(j!=current_ind){
              if(!nj_val(j) ){
                dens=pdf_empty;
                cluster_prob=prob_empty;
              } else{
                dens= surv_fn(st(jj),nu(jj),a0+isfailure2(j), l0+survtime2(j) ) ;
                //lgamma(a0+isfailure2(j)+ nu(jj) ) -(a0+isfailure2(j)+ nu(jj) )  *log(l0+survtime2(j)+ st(jj));
                for(auto it:non_na_obs(jj))
                  dens+=( log( noccu(j, it  )(eta(jj,it) ) +  alpha_vec(it))-log( ncat(it)*alpha_vec(it) + nobs(j,it) ) );
                // Rcpp::Rcout<<"at j="<<j<<" log_density="<<dens<<"\t clust prob="<<cluster_prob<<endl;
                cluster_prob=log(nj_val2(j)+dir_prec);  //log_nj_val2(j);
              }
              log_probs(j)=  dens +cluster_prob;
            }
          }
          
          if(nj_val(current_ind)==1){
            dens= pdf_empty;
            cluster_prob=prob_empty;
          } else{ 
            j=current_ind;
            // Rcpp::Rcout<<"jj= "<<jj<<" del(jj)= "<<del(jj)<<" current_ind="<<current_ind<<endl;
            dens= surv_fn(st(jj),nu(jj),a0+isfailure2(j)-nu(jj), l0+survtime2(j)-st(jj) ) ;
            //lgamma(a0+isfailure2(j) ) -(a0+isfailure2(j) )  *log(l0+survtime2(j) );
            for(auto it:non_na_obs(jj))
              dens+=( log( noccu(j, it  )(eta(jj,it) )-1 +  alpha_vec(it))-log( ncat(it)*alpha_vec(it) + nobs(j,it) -1) );
            cluster_prob=log(nj_val2(j)-1+dir_prec);
          }
          
          log_probs(current_ind)=  dens +cluster_prob;
        }
      }
      
      if(jj<n1){
        // Rcpp::Rcout<<"flag 3 0"<<endl;
        probs.zeros(nmix);
        probs(non_empty_clusters2)=probs1;
        // Rcpp::Rcout<<"flag 3 1"<<endl;
        // (probs1.t()).print("probs1 :");
        // (probs.t()).print("probs :");
      } else{
        log_DEN=log_sum_exp(log_probs);
        probs= normalise(exp(log_probs-log_DEN) ,1);
      }
      
      // log_probs.print("log_probs: ");
      // probs.print("probs: ");
      
      /*if(   gsl_fcmp(sum(probs),1.0,1e-5) ){
       // Rcpp::Rcout<< "At jj="<<jj<<"sum_prob is 0"<<endl;
       log_probs_max= max(log_probs);
       log_probs-=log_probs_max;
       
       probs=normalise(exp(log_probs) ,1);
      }*/
      
      ////check if sum of the allocation probbilities is zero
      if(sum(probs)==0) Rcpp::stop("sum(probs)=0");
      ///////////////////////////////////////////////////////
      // Rcpp::Rcout<<"flag 3 2"<<endl;
      // Rcpp::Rcout<<size(probs)<<"\t"<<size(d)<<"\t"<<size(del)<< "\t nmix= "<<nmix<<endl;
      // Rcpp::Rcout<<"jj= "<<jj<<endl;
      R::rmultinom(1, probs.begin(), nmix, d.begin());
      // (d.t()).print("d :");
      arma::uvec dd=find(d==1,1,"first");
      // (dd.t()).print("dd :");
      del(jj)=dd(  0);
      // Rcpp::Rcout<<"flag jj 4"<<endl;
      
      //updating cluster occupancies
      if(del(jj)!=current_ind){
        if(jj<n1){
          --nj_val1(current_ind);
          ++nj_val1(del(jj));
          
          survtime1(current_ind)-= st(jj);
          isfailure1(current_ind)-= nu(jj);
          
          survtime1(del(jj))+= st(jj);
          isfailure1(del(jj))+= nu(jj);
        } else{
          --nj_val2(current_ind);
          ++nj_val2(del(jj));
          
          survtime2(current_ind)-= st(jj);
          isfailure2(current_ind)-= nu(jj);
          
          survtime2(del(jj))+= st(jj);
          isfailure2(del(jj))+= nu(jj);
        }
        --nj_val(current_ind);
        ++nj_val(del(jj));
        // Rcpp::Rcout<<"flag 3 3"<<endl;
        
        
        
        for(auto it:non_na_obs(jj)){
          --noccu(current_ind, it  )(eta(jj,it) );
          ++noccu(del(jj), it  )(eta(jj,it) );
          
          --nobs(current_ind, it);
          ++nobs(del(jj), it);
          
          /*if(nobs(del(jj), it)>nj_val(del(jj)))
           Rcpp::stop("nobs(del(jj), it)>nj_val(del(jj))");*/
        }
        
        if(nj_val(del(jj)) != (nj_val2(del(jj)) + nj_val1(del(jj))) )
          Rcpp::stop("nj_val(del(jj)) != (nj_val2(del(jj)) + nj_val1(del(jj)))");
      }
    }
    // Rcpp::stop("1 MC iteration done!");
    non_empty_clusters1=(find(nj_val1)); non_empty_clusters2=(find(nj_val2));
    
    unsigned nmix1=non_empty_clusters2.n_elem;
    double dir_prec1=alpha/nmix1;
    // Rcpp::Rcout<<"allocation variables updated"<<endl;
    
    // --- UPDATE MIXTURE PROBABILITY --- //
    unsigned counttttt=0,count1=0,count2=0;
    arma::vec pi1(nmix,fill::zeros), pi2(nmix),pi1_tmp(nmix1);
    arma::vec lam1(nmix,fill::zeros), lam2(nmix,fill::zeros);
    for(j=0;j<nmix;++j){
      if(nj_val(j))
        ++counttttt;
      if(nj_val1(j)){
        ++count1;
        inds_eq_j1[j] =find(del.head(n1)==j);
      } else inds_eq_j1[j].reset();
      if(nj_val2(j)){
        ++count2;
        inds_eq_j2[j] =find(del.tail(n2)==j)+n1;
        
      } else inds_eq_j2[j].reset();
      inds_eq_j[j] =join_cols(inds_eq_j1[j],inds_eq_j2[j]);
      
      if((inds_eq_j1[j]).n_elem != nj_val1(j) || (inds_eq_j2[j]).n_elem != nj_val2(j) ){
        Rcpp::Rcout<<"inds_eq_j1[j]).n_elem="<<(inds_eq_j1[j]).n_elem<<"\t"<<"nj_val1("<<j<<")="<<nj_val1(j)<<endl;
        Rcpp::Rcout<<"inds_eq_j2[j]).n_elem="<<(inds_eq_j2[j]).n_elem<<"\t"<<"nj_val2("<<j<<")="<<nj_val2(j)<<endl;
        Rcpp::stop("Occupancy mismatch!!!");
      }
      
      ///simulate parameters corresponding to the response variable
      /* we have to do this for all non-empty clusters of X2, that's why simulating
       * lam1(j)'s also. For nj_val1(j)=0, lam1(j) will get simulated from the prior distribution.       */
      lam1(j)=randg(distr_param ( a0+ isfailure1(j), 1/(l0+survtime1(j)) ));
      lam2(j)=randg(distr_param ( a0+ isfailure2(j), 1/(l0+survtime2(j)) ));
      /////////////////////////////////////
    }
    
    
    arma::vec dir_alpha=conv_to<arma::vec>::from( nj_val2)+ dir_prec; //posterior parameters for pi2
    gsl_ran_dirichlet(r, nmix, dir_alpha.begin(), pi2.begin()); //simulate pi2 from Dirichlet
    // pi2=normalise(dir_alpha,1);
    
    
    // arma::vec dir_alpha1=conv_to<arma::vec>::from( nj_val1(non_empty_clusters1))+ dir_prec1; //only clusters with observations (cheating!)
    arma::vec dir_alpha1=conv_to<arma::vec>::from( nj_val1(non_empty_clusters2))+ dir_prec1; //posterior parameters for pi1
    gsl_ran_dirichlet(r, nmix1, dir_alpha1.begin(), pi1_tmp.begin()); //simulate pi1 from Dirichlet (restricted to the non-empty clusters of X2)
    pi1.zeros(); pi1(non_empty_clusters2)=pi1_tmp; //for programming convenience, we set length(pi1)=nmix and set the indices corresponding to the empty clusters of X2=0
    // pi1.zeros(); pi1(non_empty_clusters)=normalise(dir_alpha1,1);
    // pi1.zeros(); pi1(non_empty_clusters1)=normalise(dir_alpha1,1); //assigning weights to only non-empty clusters in the posterior
    
    arma::vec  wght2(nmix,fill::zeros), mean_pi1(nmix,fill::zeros);
    mean_pi1(non_empty_clusters2)=normalise(dir_alpha1,1);
    // Rcpp::Rcout<<"flag 101"<<endl;
    wght2(non_empty_clusters2)=mean_pi1(non_empty_clusters2)  /nj_val2(non_empty_clusters2); //E(\pi_{1h}| c_{1:n}) /n_{2h}
    // Rcpp::Rcout<<"flag 102"<<endl;
    ////////////////////////////////////////////////
    
    ///assigning prior weights to the empty clusters in group 1
    /*wght(non_empty_clusters2)= ( (n2-1+alpha)/(n1+alpha) )* (dir_alpha1/(dir_alpha(non_empty_clusters2)-1)) ;
     wght2(non_empty_clusters2)=   (dir_alpha1/(n1+alpha)) ; //weight of clusters in the clinical arm*/
    //////////////////////
    
    arma::vec  wght_x2_new=wght2(del.tail(n2));
    
    
    /////******* average treatment effects *******/////
    // Rcpp::Rcout<<"isfailure/survtime"<<endl;
    /*arma::vec theta1 = (a0 + conv_to<arma::vec>::from(isfailure1))/(l0+survtime1);
     arma::vec theta2 = (a0 + conv_to<arma::vec>::from(isfailure2))/(l0+survtime2);*/
    //////////////////////
    
    
    
    //for exponential model HR= theta_clinic/theta_synthetic . 
    //Now E(\lambda)= (a0+isfailure)/(l0+survtime). E(1/\lambda)= (l0+survtime)/(a0+isfailure-1).
    /*arma::vec HR(nmix,fill::zeros); 
     HR.zeros();
     HR(non_empty_clusters1)= (theta1(non_empty_clusters1)%(l0+survtime2(non_empty_clusters1)))/ ((a0-1) + conv_to<arma::vec>::from(isfailure2(non_empty_clusters1))) ;*/
    
    // Rcpp::Rcout<<"isfailure/survtime 2"<<endl;
    /*************************************************/
    
    
    /*Rcpp::Rcout<< "# total clusters= "<< counttttt<<endl;
     Rcpp::Rcout<< "# total clusters in clinical arm= "<< count1<<endl;
     Rcpp::Rcout<< "# total clusters in synthetic arm= "<< count2<<endl;*/
    
    ///////Generate Uniforms for model validation
    unif1= 1- exp( - st.head(n1)% lam1(del.head(n1)) );
    unif2= 1- exp( - st.tail(n2)% lam2(del.tail(n2)) );
    unif=join_cols(unif1, unif2);
    unif(censored_indices)+= randu(censored_indices.n_elem)% (1- unif(censored_indices)  ); //adjustment from cao et al (2010) Biometrics paper
    ///////////////////////////////////
    
    /****************MALA UPDATE OF a0,l0**********/
    arma::vec  proposed_params(2), del_pi_current, del_pi_proposed;
    // current_params(0)=loga0; current_params(1)=l0;
    
    del_pi_current= delpi(nj_val1, nj_val2, 
                          isfailure1, isfailure2,
                          survtime1,  survtime2,
                          current_params, aa, ba,  al,  bl);
    proposed_params = current_params +tau* del_pi_current +sqrt(2*tau) * arma::vec(2, fill::randn);
    
    del_pi_proposed= delpi(nj_val1, nj_val2, 
                           isfailure1, isfailure2,
                           survtime1,  survtime2,
                           proposed_params, aa, ba,  al,  bl);
    arma::vec del_tmp= proposed_params - current_params -tau*del_pi_current;
    double log_q01= -dot(del_tmp, del_tmp )/(4*tau);
    
    
    del_tmp= current_params- proposed_params -tau*del_pi_proposed;
    // del_tmp.print("del_tmp for q10");
    double log_q10= -dot(del_tmp, del_tmp )/(4*tau);
    
    double log_pi_current=logpi(nj_val1, nj_val2, isfailure1, isfailure2,
                                survtime1,  survtime2, current_params, aa, ba,  al,  bl), 
                                log_pi_proposed=logpi(nj_val1, nj_val2, isfailure1, isfailure2,
                                                      survtime1,  survtime2, proposed_params, aa, ba,  al,  bl);
    
    
    double log_prob= log_pi_proposed + log_q10 -( log_pi_current +log_q01);
    /*current_params.print("current_params");
     proposed_params.print("proposed_params");
     Rcpp::Rcout<<"log_pi_proposed = "<< log_pi_proposed<<" log_pi_current= "<<log_pi_current<< endl;
     Rcpp::Rcout<<"log_q10 = "<< log_q10<<" log_q01= "<<log_q01<< endl;*/
    
    
    if(log(randu()) <log_prob){
      current_params=proposed_params;
      a0=exp(proposed_params(0)); l0=exp(proposed_params(1));
      ++acceptance;
    }
    //////////////////////////////////////////////// 
    
    
    thincheck = i - std::floor(i/thin) * thin; // % operator stolen by arma
    if(!thincheck ) {
      Rcpp::Rcout<<"Iteration: "<<i<<" Acceptance prob= "<< ((double)acceptance)/i<<" log_prob="<<log_prob<< endl;
    }
    
    
    int remainder= (i+1 );
    int quotient= (int) std::floor(remainder/thin);
    remainder-= (quotient*thin) ;
    
    if(remainder==0){
      alloc_var_mat.row(quotient-1)=del.t();
      pimat1.row(quotient-1)=pi1.t();
      pimat2.row(quotient-1)=pi2.t();
      // weights.row(quotient-1)=wght_x2.t(); //Ratio of cluster probabilities
      weights2.row(quotient-1)=wght_x2_new.t(); //probability of the cluster
      thetamat1.row(quotient-1)=lam1.t(); // theta1.t();
      thetamat2.row(quotient-1)=lam2.t(); //theta2.t();
      hyperparams.row(quotient-1)= current_params.t();
      // HRmat.row(quotient-1)=HR.t();
      unifmat.row(quotient-1)=unif.t();
    }
  }
  
  // delete[] inds_eq_j;   delete[] inds_eq_j1;   delete[] inds_eq_j2;
  gsl_rng_free (r);
  
  
  return Rcpp::List::create(Rcpp::Named("pimat1") =pimat1,
                            Rcpp::Named("pimat2") =pimat2,
                            Rcpp::Named("Allocation variables") = alloc_var_mat,
                            // Rcpp::Named("Weights")=weights,
                            Rcpp::Named("Weights2")=weights2,
                            Rcpp::Named("Avg_response_clinical")=thetamat1,
                            Rcpp::Named("Avg_response_control")=thetamat2,
                            Rcpp::Named("Exponential_hyperparams")=exp(hyperparams),
                            Rcpp::Named("Unifs")=unifmat
  );
}

// [[Rcpp::export]]
Rcpp::List common_atoms_cat_lognormal(const unsigned nmix, arma::uvec ncat,
                                      const double a0, const double df0,  const double mu_m,const double mu_v, const double b_m, const double b_v,//nomral and lognormal hyperparameters for the response
                                      const int nrun,const  int burn, const int thin, 
                                      arma::umat eta1, arma::umat eta2,
                                      arma::mat eta_cont1, arma::mat eta_cont2,
                                      arma::vec st1, arma::uvec nu1, //responses for the treatment
                                      arma::vec st2, arma::uvec nu2, //responses for the synthetic control
                                      Rcpp::List non_na_obs1,// NA observations in the categorical variable
                                      Rcpp::List non_na_obs1_cont, // NA observations in the cont variable
                                      arma::uvec del1, arma::uvec del2, 
                                      const arma::vec &del_range_lognorm, const unsigned nleapfrog_lognorm,
                                      const arma::vec &alpha_hyper, const arma::vec &del_range_alp1, const unsigned nleapfrog_alp1,
                                      const arma::vec &del_range_alp2, const unsigned nleapfrog_alp2){
  unsigned acceptance=0, acceptance_alph1=0, acceptance_alph2=0;
  double mu_alp, sig_alp;
  sig_alp=log1p(alpha_hyper(1) /gsl_pow_2(alpha_hyper(0))); mu_alp=log(alpha_hyper(0))-sig_alp/2;sig_alp=sqrt(sig_alp);
  
  const unsigned n1=st1.n_elem,n2=st2.n_elem, k=eta1.n_cols, k_cont=eta_cont1.n_cols;
  const unsigned n=n1+n2;
  double t1ppp;
  
  // double mu0=sqrt(mu_v)* randn()+mu_m, beta0=exp(sqrt(b_v)* randn()+b_m);
  double mu0=mu_m, beta0=exp( b_v/2 +b_m);
  arma::vec current_params(2);  current_params(0)=mu0; current_params(1)=log(beta0);
  
  arma::umat eta=join_cols( eta1, eta2 ); //categorical covariates
  arma::mat eta_cont=join_cols( eta_cont1, eta_cont2 ); //cont covariates
  arma::mat eta_cont_sq=square(eta_cont);
  arma::uvec del=join_cols(del1,del2), nu=join_cols(nu1,nu2) ;
  arma::vec st=join_cols(st1,st2);
  const double max_st=5*max(st);
  // Rcpp::Rcout<<"size(del)"<<size(del)<<endl;
  field<arma::uvec> non_na_obs(n), non_na_obs_cont(n);
  unsigned cat_na_count=0;
  for(unsigned i=0;i<n;++i){ /**changing the data-type to field from input list**/ 
    non_na_obs(i)=Rcpp::as<arma::uvec>(non_na_obs1[i]);
    cat_na_count+= (k-(non_na_obs(i)).n_elem);
    
    non_na_obs_cont(i)=Rcpp::as<arma::uvec>(non_na_obs1_cont[i]);
  }
  bool cat_na= (cat_na_count==(eta.n_rows* eta.n_cols)); //TRUE if cat covariate is NULL
  
  unsigned j;
  
  /////Define MCMC storage matrices
  unsigned n_mc=std::floor((nrun+burn)/thin);
  arma::umat alloc_var_mat(n_mc, n); 
  arma::mat weights2(n_mc, n2), weights(n_mc, nmix);
  arma::mat pimat1(n_mc, nmix), pimat2(n_mc, nmix), dir_alpha_mat(n_mc, 2);
  arma::mat lognormal_mu1(n_mc, nmix), lognormal_sig1(n_mc, nmix), lognormal_mu2(n_mc, nmix), lognormal_sig2(n_mc, nmix);
  arma::mat hyperparams(n_mc, 2);
  arma::mat unifmat(n_mc, n);
  
  ///////////////////////
  
  
  // --- initialize loop objects --- //
  uvec d(nmix); ////multinomial indicator for each sample
  arma::vec probs(nmix),log_probs(nmix), alpha_vec(k,fill::value(1.0)); ////assignment probability for each sample
  
  ////set atoms for the categorical covariate
  
  unsigned count_cat=0;
  eta.for_each( [ &count_cat](umat::elem_type val) {count_cat+= (!std::isfinite(val) )  ; } );
  
  field<arma::uvec> noccu(nmix,k);
  for(unsigned i=0;i<nmix;++i){
    for(unsigned j=0;j<k;++j){
      if(!cat_na)
        if(ncat(j)<2){
          Rcpp::Rcout<<"j= "<<j<<" ncat(j)<2"<<endl;
          Rcpp::stop("");
        }
      
      noccu(i,j).set_size(ncat(j));
      noccu(i,j).zeros();
    }
  }
  //////
  
  ////set atoms for the cont covariates
  arma::umat nj_x(nmix,k_cont,fill::zeros); //counts the number of non-missing observations for the cont covariates in each atom 
  // arma::mat df_post_x(nmix,k_cont), alpha_post_x(nmix,k_cont), beta_post_x(nmix,k_cont), 
  arma::mat sum_j_x(nmix,k_cont, fill::zeros), ss_j_x(nmix,k_cont, fill::zeros);
  const double df_x=1, alpha_x=k_cont+ 30.0, beta_x=1.0, mu_x=0.0; //This is for the simulations
  // const double df_x=3, alpha_x=k_cont+ 1.0, beta_x=.50, mu_x=0.0; //This is for the diagnostic
  // df_post_x.fill(df_x); alpha_post_x.fill(alpha_x); beta_post_x.fill(beta_x);
  ///////////////////////////////////////
  
  bool thincheck, printcheck;
  field<arma::uvec> inds_eq_j(nmix), inds_eq_j1(nmix), inds_eq_j2(nmix);
  
  ///// Setting-up GSL random number generator for sampling from dirichlet
  const gsl_rng_type * T;
  gsl_rng * r;
  /** create a generator chosen by the
   environment variable GSL_RNG_TYPE **/
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  gsl_rng_set(r,500);
  /*************************************/
  
  /********initialize parameters corresponding to clusters*******/
  arma::vec nj_val1(nmix),nj_val2(nmix), nj_val(nmix);//occupancy number corresponding to each cluster
  arma::vec survtime1(nmix,fill::zeros), ss_survtime1(nmix,fill::zeros), survtime2(nmix,fill::zeros), ss_survtime2(nmix,fill::zeros) ;
  
  Rcpp::Rcout<<"cluster occupancy finding starts!!"<<endl;
  for(unsigned j=0;j<nmix;++j){
    inds_eq_j1[j] =find(del.head(n1)==j); nj_val1(j)=(inds_eq_j1[j]) .n_elem; 
    inds_eq_j2[j] =find(del.tail(n2)==j)+n1; nj_val2(j)=(inds_eq_j2[j]) .n_elem; 
    inds_eq_j[j] =join_cols(inds_eq_j1[j],inds_eq_j2[j]);
    nj_val(j)=nj_val1(j)+nj_val2(j);
    
    survtime1(j)= sum(st(inds_eq_j1[j]) ); //isfailure1(j)=sum(nu(inds_eq_j1[j] ));
    survtime2(j)= sum(st(inds_eq_j2[j]) ); //isfailure2(j)=sum(nu(inds_eq_j2[j] ));
    ss_survtime1(j)= ssq(st(inds_eq_j1[j]) ); 
    ss_survtime2(j)= ssq(st(inds_eq_j2[j]) ); 
    
    if(nj_val(j)){
      Rcpp::Rcout<<"n_inds at j="<<j<<" is "<<(inds_eq_j[j]) .n_elem<<endl;
      /*if(nj_val2(j))
       log_nj_val2(j)=log(nj_val2(j)+dir_prec);*/
      
      for(auto it1:inds_eq_j[j]){//it1 iterates through cluster membership indicators
        for(auto it2:non_na_obs(it1)){//it2 iterates through non-missing variables of the categorical covariate for each observation
          // Rcpp::Rcout<<"it1= "<<it1<<" it2 ="<<it2<<" eta(it1,it2)= "<<eta(it1,it2)<<endl;
          ++noccu(j,it2)(eta(it1,it2));
        }
        for(auto it2:non_na_obs_cont(it1)){//it2 iterates through non-missing variables of the continuous covariates for each observation
          // Rcpp::Rcout<<"it1= "<<it1<<" it2 ="<<it2<<" eta_cont(it1,it2)= "<<eta_cont(it1,it2)<<endl;
          ++nj_x(j,it2);
          sum_j_x(j,it2)+= eta_cont(it1,it2);
          ss_j_x(j,it2)+= eta_cont_sq(it1,it2);// gsl_pow_2( eta_cont(it1,it2));
        }
      } 
    }
  }
  if(sum_j_x.has_nan()|| ss_j_x.has_nan())
    Rcpp::stop("sum_j_x.has_nan()|| ss_j_x.has_nan() ");
  const arma::vec st_original=st;
  const arma::uvec censored_indices= find(nu==0);
  
  //----------------------------------------------//
  
  ///find the number of non-NA observations for each variable
  arma::umat nobs(nmix,k,fill::zeros);
  for(unsigned l=0;l<nmix;++l)
    for(unsigned j=0;j<k;++j){
      nobs(l,j)= sum(noccu(l,j));
      // Rcpp::Rcout<<"nobs(j,l)= "<<nobs(l,j)<< " nj_val(j)= "<<nj_val(l)<<endl;
    }
    //////////////////////////////////////////////////////////  
    
    Rcpp::Rcout<<"loop starts"<<endl;
  // --- loop --- //
  arma::uvec non_empty_clusters1=(find(nj_val1)), non_empty_clusters2=(find(nj_val2));
  unsigned nmix1=non_empty_clusters2.n_elem;
  double alpha1=1, alpha2=1;
  double dir_prec1=alpha1/nmix1, dir_prec2=alpha2/nmix;
  double prob_empty1= log(dir_prec1), prob_empty2= log(dir_prec2);
  double df_post, alpha_post, tmp, ss_j, survtime_j, mean_j, mu_post, beta_post, sigma_post;
  arma::vec unif1, unif2, unif;
  
  for(int i=0; i<nrun+burn; ++i){
    /////////////Sanity check /////////////
    /*Rcpp::Rcout<<"flag sanity"<<endl;
     arma::uvec nobs_check(k,fill::zeros);
     for(j=0;j<k;++j){
     for(unsigned l=0;l<nmix;++l)
     nobs_check(j)+= sum(noccu(l,j));
     if(nobs_check(j) != nobs(j)){
     Rcpp::Rcout<<"j = "<<j<<"nobs_check(j) != nobs(j)"<<endl;
     Rcpp::Rcout<<"nobs_check(j)= "<<nobs_check(j) <<"  nobs(j)= "<<nobs(j)<<endl;
     Rcpp::stop("");
     }
     }*/
    //////////////////////////
    
    // --- Augment censored observations ---//
    for(j=0; j<nmix; ++j){//Iterate over the clusters
      ///notations mostly follow Wikipedia conjugate prior NIG
      
      ////UPDATE for the experimental arm
      if(nj_val1(j)){
        // Rcpp::Rcout<<"df0= "<<df0<<" a0= "<<a0<< " mu0= "<<mu0 <<" beta0= "<<beta0<< endl;
        df_post=df0+ ( (double) nj_val1(j))-1; alpha_post=a0+ (( (double) nj_val1(j))-1) /2.0;
        tmp=( ((double) (nj_val1(j)-1))*df0 )/df_post;
        for(auto jj:inds_eq_j1[j] ){//Iterate over the observations
          if(!nu(jj)){//censored observations
            if(nj_val1(j)>1){
              ss_j= ss_survtime1(j) - gsl_pow_2(st(jj) ), survtime_j = survtime1(j)- st(jj);
              mean_j= survtime_j/(nj_val1(j) - 1 );
            } else ss_j=  survtime_j =  mean_j= 0.0;
            mu_post = (df0*mu0 + survtime_j) / df_post;
            beta_post=beta0+  (ss_j -  ((double) nj_val1(j)-1) * gsl_pow_2(mean_j )  + tmp* gsl_pow_2(mean_j - mu0) )/2.0;
            sigma_post =sqrt( beta_post * (df_post+1)/(df_post *alpha_post) );
            st(jj)=r_trunclst(2*alpha_post,  mu_post, sigma_post,  st_original(jj), std::numeric_limits<double>::max());
            st(jj)= GSL_MIN_DBL(st(jj), max_st);
            /*if(st(jj)>log(100*52))
             Rcpp::Rcout<<"jj="<<jj << " st_original(jj)="<<st_original(jj)<< " st(jj)="<<st(jj)<<endl;*/
            ss_survtime1(j)=ss_j + gsl_pow_2(st(jj)); survtime1(j)= survtime_j +st(jj);
          }
        }
      }
      if(nj_val2(j)){
        ////UPDATE for the synthetic control arm
        df_post=df0+ nj_val2(j)-1, alpha_post=a0+ ((double) nj_val2(j)-1)/2.0;
        tmp=( (nj_val2(j)-1)*df0 )/df_post;
        for(auto jj:inds_eq_j2[j] ){//Iterate over the observations
          if(!nu(jj)){//censored observations
            if(nj_val2(j)>1){
              ss_j= ss_survtime2(j) - gsl_pow_2(st(jj) ), survtime_j = survtime2(j)- st(jj);
              mean_j= survtime_j/(nj_val2(j) - 1 );
            } else ss_j=  survtime_j =  mean_j= 0.0;
            mu_post = (df0*mu0 + survtime_j) / df_post;
            beta_post=beta0+  (ss_j -  ((double) nj_val2(j)-1) * gsl_pow_2(mean_j )  + tmp* gsl_pow_2(mean_j - mu0) )/2.0;
            sigma_post =sqrt( beta_post * (df_post+1)/(df_post *alpha_post) );
            st(jj)=r_trunclst(2*alpha_post,  mu_post, sigma_post,  st_original(jj), std::numeric_limits<double>::max()); 
            //df_post is not the degrees of freedom, its the precision parameter
            st(jj)= GSL_MIN_DBL(st(jj), max_st );
            /*if(st(jj)>log(100*52))
             Rcpp::Rcout<<"jj="<<jj << " st_original(jj)="<<st_original(jj)<< " st(jj)="<<st(jj)<<endl;*/
            ss_survtime2(j)=ss_j + gsl_pow_2(st(jj)); survtime2(j)= survtime_j +st(jj);
          }
        }
      }
    }
    // Rcpp::Rcout<<"data augmentation done!"<<endl;
    
    
    /***** UPDATE ALLOCATION VARIABLES *****/
    arma::vec normal_exp(nmix),probs1;//exp2(nmix),exp1(nmix);
    double  log_probs_max;//,normal_exp;
    long double log_DEN;
    for(unsigned jj=0;jj<n;++jj){
      double dens,cluster_prob, pdf_empty,  st_sq=gsl_pow_2(st(jj));
      log_probs.fill(datum::log_min);
      // Rcpp::Rcout<<"log_pdf_empty="<<pdf_empty<<endl;
      // Rcpp::stop("");
      // Rcpp::Rcout<<"jj="<<jj<<"\t"<<"del(jj)="<<del(jj)<<endl;
      unsigned current_ind=del(jj); //find the current index of data point jj
      if(jj<n1){//this is for G_1
        // Rcpp::Rcout<<" G1"<<endl;
        for(j=0;j<nmix;++j){
          // Rcpp::Rcout<<"j= "<<j<<" nj_val1(j)="<<nj_val1(j)<<endl;
          if(nj_val2(j)){
            if(j!=current_ind){
              dens= surv_fn_lognorm( st_original(jj),  nu(jj), ss_survtime1(j), survtime1(j),
                                     df0, a0, mu0, beta0, nj_val1(j));
              
              // Rcpp::Rcout<<"jj= "<<jj<< " nu(jj)= "<<nu(jj)<< " st(jj)= "<<st(jj)<< "st_original(jj)= "<<st_original(jj) << " j= "<<j<<" nj_val1(j)= "<<nj_val1(j) <<" dens= "<<dens<<endl;
              
              for(auto it:non_na_obs(jj))//density part from categorical covs
                dens+=( log( noccu(j, it  )(eta(jj,it) )+  alpha_vec(it))-log(nobs(j,it) + ncat(it)*alpha_vec(it)  ) );
              
              for(auto it:non_na_obs_cont(jj))///density part from cont covs
                dens+=( post_t_dens(eta_cont(jj,it) , ss_j_x(j,it), sum_j_x(j,it), 
                                    df_x, alpha_x, mu_x, beta_x, nj_x(j,it) ) );
              // Rcpp::Rcout<<"at j="<<j<<" log_density="<<dens<<"\t clust prob="<<cluster_prob<<endl;
              cluster_prob=log(nj_val1(j)+dir_prec1);  //log_nj_val2(j);
              
              log_probs(j)=  dens +cluster_prob;
              // Rcpp::Rcout<<"log_probs "<<j<< "\t"<<log_probs(j)<<endl;
            }
          }
        }
        //calculating allocation probabilities
        j=current_ind;
        if(nj_val(j)==0 || nj_val2(j)==0)   Rcpp::stop("nj_val(current_ind)=0 || nj_val2(current_ind)=0 jj in G_1");
        // Rcpp::Rcout<<"flag 1"<<endl;
        dens= surv_fn_lognorm( st_original(jj),  nu(jj), ss_survtime1(j)-st_sq, survtime1(j) -st(jj),
                               df0, a0, mu0, beta0, nj_val1(j)-1);
        // Rcpp::Rcout<<"jj= "<<jj<< " nu(jj)= "<<nu(jj)<< " st(jj)= "<<st(jj)<< "st_original(jj)= "<<st_original(jj) << " j= "<<j<<" nj_val1(current_ind)= "<<nj_val1(j) <<" dens= "<<dens<<endl;
        for(auto it:non_na_obs(jj)){
          /*double tmp_dens=( log( noccu(j, it  )(eta(jj,it) )-1 +  alpha_vec(it))-log( ncat(it)*alpha_vec(it) + nobs(j,it) -1) );
           Rcpp::Rcout<<"tmp_dens cat="<<tmp_dens<<endl;
           dens+=tmp_dens;*/
          dens+=( log( noccu(j, it  )(eta(jj,it) )-1 +  alpha_vec(it))-log( ncat(it)*alpha_vec(it) + nobs(j,it) -1) );
        }
        for(auto it:non_na_obs_cont(jj)){///density part from cont covs
          /*double tmp_dens=( post_t_dens(eta_cont(jj,it) , ss_j_x(j,it)-eta_cont_sq(jj,it), sum_j_x(j,it)- eta_cont(jj,it), 
           df_x, alpha_x, mu_x, beta_x, nj_x(j,it)-1 ) );
           Rcpp::Rcout<<"eta_cont(jj,it)="<<eta_cont(jj,it)<<" ss_j_x(j,it)="<<ss_j_x(j,it)<<" sum_j_x(j,it)="<<sum_j_x(j,it)<<" nj_x(j,it)="<<nj_x(j,it)<<" tmp_dens cont="<<tmp_dens<<endl;
           dens+=tmp_dens;*/
          dens+=( post_t_dens(eta_cont(jj,it) , ss_j_x(j,it)-eta_cont_sq(jj,it), sum_j_x(j,it)- eta_cont(jj,it), 
                              df_x, alpha_x, mu_x, beta_x, nj_x(j,it)-1 ) );
        }
        
        // Rcpp::Rcout<<"at j="<<j<<" log_density="<<dens<<"\t clust prob="<<cluster_prob<<endl;
        cluster_prob=log(nj_val1(j)-1+dir_prec1);  //log_nj_val2(j);
        // Rcpp::Rcout<<"flag 2"<<endl;
        log_probs(j)=  dens +cluster_prob;
        // Rcpp::Rcout<<"at current_ind dens="<<dens<<" cluster_prob="<<cluster_prob<<" log_probs(current_ind)"<<log_probs(j)<<endl;
        
        // (log_probs.t()).print("log_probs:");
        
        // non_empty_clusters2=  find(nj_val2);
        arma::vec log_probs1=log_probs(non_empty_clusters2);
        log_DEN=log_sum_exp(log_probs1);
        probs1=exp(log_probs1-log_DEN); //normalise(exp(log_probs1-log_DEN) ,1);
        // Rcpp::Rcout<<"flag 3"<<endl;
      } else{//this is for G_2
        // Rcpp::Rcout<<" G2"<<endl;
        j=current_ind;
        if(nj_val1(j) && nj_val2(j)==1){//the degenerate case
          log_probs.fill(datum::log_min);
          log_probs(j)=datum::log_max;
        } else if(nj_val1(j) && !nj_val2(j)  ){//the infeasible case
          Rcpp::stop("nj_val(current_ind)=0 || nj_val2(current_ind)=0 jj in G_2");
        } else {//other possible cases
          double pdf_empty = surv_fn_lognorm( st_original(jj),  nu(jj), 0.0, 0.0, df0, a0, mu0, beta0, 0);
          //log(dens) for the categorical covariates for an empty cluster is 0. Verify!
          for(auto it:non_na_obs_cont(jj))///density part from cont covs
            pdf_empty+=( post_t_dens(eta_cont(jj,it) , 0.0, 0.0, df_x, alpha_x, mu_x, beta_x, 0 ) );
          
          for(j=0;j<nmix;++j){
            if(j!=current_ind){
              if(!nj_val(j) ){
                dens=pdf_empty;
                cluster_prob=prob_empty2;
              } else{
                dens= surv_fn_lognorm( st_original(jj),  nu(jj), ss_survtime2(j), survtime2(j),
                                       df0, a0, mu0, beta0, nj_val2(j));
                for(auto it:non_na_obs(jj))
                  dens+=( log( noccu(j, it  )(eta(jj,it) ) +  alpha_vec(it))-log( ncat(it)*alpha_vec(it) + nobs(j,it) ) );
                
                for(auto it:non_na_obs_cont(jj))///density part from cont covs
                  dens+=( post_t_dens(eta_cont(jj,it) , ss_j_x(j,it), sum_j_x(j,it), 
                                      df_x, alpha_x, mu_x, beta_x, nj_x(j,it) ) );
                // Rcpp::Rcout<<"at j="<<j<<" log_density="<<dens<<"\t clust prob="<<cluster_prob<<endl;
                cluster_prob=log(nj_val2(j)+dir_prec2);  //log_nj_val2(j);
              }
              log_probs(j)=  dens +cluster_prob;
            }
          }
          
          if(nj_val(current_ind)==1){
            dens= pdf_empty;
            cluster_prob=prob_empty2;
          } else{
            j=current_ind;
            // Rcpp::Rcout<<"jj= "<<jj<<" del(jj)= "<<del(jj)<<" current_ind="<<current_ind<<endl;
            dens= surv_fn_lognorm( st_original(jj),  nu(jj), ss_survtime2(j)-st_sq, survtime2(j) -st(jj),
                                   df0, a0, mu0, beta0, nj_val2(j)-1);
            for(auto it:non_na_obs(jj))
              dens+=( log( noccu(j, it  )(eta(jj,it) )-1 +  alpha_vec(it))-log( ncat(it)*alpha_vec(it) + nobs(j,it) -1) );
            for(auto it:non_na_obs_cont(jj))///density part from cont covs
              dens+=( post_t_dens(eta_cont(jj,it) , ss_j_x(j,it)-eta_cont_sq(jj,it), sum_j_x(j,it)- eta_cont(jj,it), 
                                  df_x, alpha_x, mu_x, beta_x, nj_x(j,it)-1 ) );
            cluster_prob=log(nj_val2(j)-1+dir_prec2);
          }
          
          log_probs(current_ind)=  dens +cluster_prob;
        }
      }
      
      if(jj<n1){
        // Rcpp::Rcout<<"flag 3 0"<<endl;
        probs.zeros(nmix);
        probs(non_empty_clusters2)=probs1;
        // Rcpp::Rcout<<"flag 3 1"<<endl;
        // (probs1.t()).print("probs1 :");
        // (probs.t()).print("probs :");
      } else{
        log_DEN=log_sum_exp(log_probs);
        probs= exp(log_probs-log_DEN);//normalise(exp(log_probs-log_DEN) ,1);
      }
      
      // log_probs.print("log_probs: ");
      // probs.print("probs: ");
      
      ////check if sum of the allocation probabilities is zero
      if(sum(probs)==0) Rcpp::stop("sum(probs)=0");
      ///////////////////////////////////////////////////////
      gsl_ran_multinomial (r, nmix, 1, probs.begin(), d.begin());
      // R::rmultinom(1, probs.begin(), nmix, d.begin());
      // (d.t()).print("d :");
      arma::uvec dd=find(d==1,1,"first");
      // (dd.t()).print("dd :");
      del(jj)=dd(  0);
      // Rcpp::Rcout<<"flag jj 4"<<endl;
      
      //updating cluster occupancies
      if(del(jj)!=current_ind){
        if(jj<n1){
          --nj_val1(current_ind);
          ++nj_val1(del(jj));
          
          survtime1(current_ind)-= st(jj);
          ss_survtime1(current_ind)-= st_sq;
          
          survtime1(del(jj))+= st(jj);
          ss_survtime1(del(jj))+= st_sq;
        } else{
          --nj_val2(current_ind);
          ++nj_val2(del(jj));
          
          survtime2(current_ind)-= st(jj);
          ss_survtime2(current_ind)-= st_sq;
          
          survtime2(del(jj))+= st(jj);
          ss_survtime2(del(jj))+= st_sq;
        }
        --nj_val(current_ind);
        ++nj_val(del(jj));
        // Rcpp::Rcout<<"flag 3 3"<<endl;
        
        //update atoms of the categorical covs
        for(auto it:non_na_obs(jj)){ 
          --noccu(current_ind, it  )(eta(jj,it) );
          ++noccu(del(jj), it  )(eta(jj,it) );
          
          --nobs(current_ind, it);
          ++nobs(del(jj), it);
          
          /*if(nobs(del(jj), it)>nj_val(del(jj)))
           Rcpp::stop("nobs(del(jj), it)>nj_val(del(jj))");*/
        }
        /////////////
        
        //update atoms of the cont covs
        arma::uvec tmp_current_ind={current_ind}, tmp_del_jj={del(jj)}, tmp_jj={jj};;
        
        nj_x(tmp_current_ind,non_na_obs_cont(jj))-=1;
        nj_x(tmp_del_jj,non_na_obs_cont(jj))+=1;
        
        ss_j_x(tmp_current_ind,non_na_obs_cont(jj))-= eta_cont_sq(tmp_jj,non_na_obs_cont(jj));
        ss_j_x(tmp_del_jj,non_na_obs_cont(jj))+= eta_cont_sq(tmp_jj,non_na_obs_cont(jj));
        
        sum_j_x(tmp_current_ind,non_na_obs_cont(jj))-= eta_cont(tmp_jj,non_na_obs_cont(jj));
        sum_j_x(tmp_del_jj,non_na_obs_cont(jj))+= eta_cont(tmp_jj,non_na_obs_cont(jj));
        //////////////////////////////
      }
    }
    // Rcpp::stop("1 MC iteration done!");
    non_empty_clusters1=(find(nj_val1)); non_empty_clusters2=(find(nj_val2));
    
    unsigned nmix1=non_empty_clusters2.n_elem;
    double dir_prec1=alpha1/nmix1;
    // Rcpp::Rcout<<"allocation variables updated"<<endl;
    
    /*** UPDATE MIXTURE PROBABILITY ***/
    unsigned counttttt=0,count1=0,count2=0;
    arma::vec pi1(nmix,fill::zeros), pi2(nmix),pi1_tmp(nmix1);
    arma::vec mu1(nmix,fill::zeros), sig1(nmix,fill::zeros), mu2(nmix,fill::zeros), sig2(nmix,fill::zeros);
    
    for(j=0; j<nmix; ++j){
      if(nj_val(j))
        ++counttttt;
      if(nj_val1(j)){
        ++count1;
        inds_eq_j1[j] =find(del.head(n1)==j);
      } else inds_eq_j1[j].reset();
      if(nj_val2(j)){
        ++count2;
        inds_eq_j2[j] =find(del.tail(n2)==j)+n1;
      } else inds_eq_j2[j].reset();
      inds_eq_j[j] =join_cols(inds_eq_j1[j],inds_eq_j2[j]);
      
      if((inds_eq_j1[j]).n_elem != nj_val1(j) || (inds_eq_j2[j]).n_elem != nj_val2(j) ){
        Rcpp::Rcout<<"inds_eq_j1[j]).n_elem="<<(inds_eq_j1[j]).n_elem<<"\t"<<"nj_val1("<<j<<")="<<nj_val1(j)<<endl;
        Rcpp::Rcout<<"inds_eq_j2[j]).n_elem="<<(inds_eq_j2[j]).n_elem<<"\t"<<"nj_val2("<<j<<")="<<nj_val2(j)<<endl;
        Rcpp::stop("Occupancy mismatch!!!");
      }
      
      ///simulate parameters corresponding to the response variable
      /* we have to do this for all non-empty clusters of X2, that's why simulating
       * lam1(j)'s also. For nj_val1(j)=0, lam1(j) will get simulated from the prior distribution.       */
      
      // Rcpp::Rcout<<"Treatment"<<endl;
      arma::vec tmpvec=sim_lognorm_params(ss_survtime1(j) , survtime1(j), df0, a0, mu0, beta0, nj_val1(j));
      mu1(j)=tmpvec(0); sig1(j)=tmpvec(1);
      
      // Rcpp::Rcout<<"Synthetic j= "<<j<<" ss_survtime2(j)= "<<ss_survtime2(j)<<" survtime2(j)= "<< survtime2(j)  <<endl;
      tmpvec=sim_lognorm_params(ss_survtime2(j) , survtime2(j), df0, a0, mu0, beta0, nj_val2(j) );
      mu2(j)=tmpvec(0); sig2(j)=tmpvec(1);
      // if(sig2(j)>100)
      //   Rcpp::Rcout<<"j= "<<j<<" sig= "<<sig2(j) <<" nj_val2(j)= "<<nj_val2(j)<< "ss_survtime2(j)= "<< ss_survtime2(j)<<endl;
      /////////////////////////////////////
    }
    
    
    arma::vec dir_alpha= nj_val2+ dir_prec2; //posterior parameters for pi2
    gsl_ran_dirichlet(r, nmix, dir_alpha.begin(), pi2.begin()); //simulate pi2 from Dirichlet
    // pi2=normalise(dir_alpha,1);
    
    
    // arma::vec dir_alpha1=conv_to<arma::vec>::from( nj_val1(non_empty_clusters1))+ dir_prec1; //only clusters with observations (cheating!)
    arma::vec dir_alpha1= nj_val1(non_empty_clusters2)+ dir_prec1; //posterior parameters for pi1
    gsl_ran_dirichlet(r, nmix1, dir_alpha1.begin(), pi1_tmp.begin()); //simulate pi1 from Dirichlet (restricted to the non-empty clusters of X2)
    pi1.zeros(); pi1(non_empty_clusters2)=pi1_tmp; //for programming convenience, we set length(pi1)=nmix and set the indices corresponding to the empty clusters of X2=0
    // pi1.zeros(); pi1(non_empty_clusters)=normalise(dir_alpha1,1);
    // pi1.zeros(); pi1(non_empty_clusters1)=normalise(dir_alpha1,1); //assigning weights to only non-empty clusters in the posterior
    
    arma::vec  wght2(nmix,fill::zeros), mean_pi1(nmix,fill::zeros);
    mean_pi1(non_empty_clusters2)=normalise(dir_alpha1-dir_prec1,1);
    // Rcpp::Rcout<<"flag 101"<<endl;
    wght2(non_empty_clusters2)=mean_pi1(non_empty_clusters2)  /nj_val2(non_empty_clusters2); //E(\pi_{1h}| c_{1:n}) /n_{2h}
    // wght2.t().print("wght2:");
    // wght2(non_empty_clusters2)=( (n2-1+alpha)/(n1+alpha) )* (dir_alpha1/(dir_alpha(non_empty_clusters2)-1)); //E(\pi_{1h}/\pi_{2h}| c_{1:n}) 
    // Rcpp::Rcout<<"flag 102"<<endl;
    ////////////////////////////////////////////////
    
    ///assigning prior weights to the empty clusters in group 1
    /*wght(non_empty_clusters2)= ( (n2-1+alpha)/(n1+alpha) )* (dir_alpha1/(dir_alpha(non_empty_clusters2)-1)) ;
     wght2(non_empty_clusters2)=   (dir_alpha1/(n1+alpha)) ; //weight of clusters in the clinical arm*/
    //////////////////////
    
    arma::vec  wght_x2_new=wght2(del.tail(n2));
    
    
    /*****Generate Uniforms for model validation *****/
    unif1= normcdf(st_original.head(n1), mu1(del.head(n1)) , sqrt (sig1(del.head(n1))) ) ;
    unif2= normcdf(st_original.tail(n2), mu2(del.tail(n2)) , sqrt(sig2(del.tail(n2))) );
    unif=join_cols(unif1, unif2);
    unif(censored_indices)+= randu(censored_indices.n_elem)% (1- unif(censored_indices)  ); //adjustment from cao et al (2010) Biometrics paper
    /*************************************************/
    
    
    /****************HMC UPDATE OF mu0, beta0**********/
    update_lognorm_hyper(nj_val1, nj_val2, survtime1, ss_survtime1,
                         survtime2, ss_survtime2,
                         a0, df0, mu_m, mu_v, b_m, b_v,
                         del_range_lognorm, nleapfrog_lognorm, current_params, acceptance);
    mu0=current_params(0), beta0=exp(current_params(1));
    /*******************************************************/
    
    
    /****** HMC UPDATE of Dirichlet mixture parameters**********/
    //update alpha_1 (for the nested mixture)
    double l_alpha=log(alpha1); 
    update_alpha(nmix1,   nj_val1(non_empty_clusters1), mu_alp, sig_alp, del_range_alp1, nleapfrog_alp1, l_alpha, acceptance_alph1);
    alpha1=exp(l_alpha);    dir_prec1=alpha1/nmix1; 
    //////////////////////////
    
    //update alpha_2 (for the rwd mixture)
    l_alpha=log(alpha2); 
    update_alpha(nmix,   nj_val2(non_empty_clusters2), mu_alp, sig_alp, del_range_alp2, nleapfrog_alp2, l_alpha, acceptance_alph2);
    alpha2=exp(l_alpha);    dir_prec2=alpha2/nmix; prob_empty2= log(dir_prec2);
    //////////////////////////
    /***********************************************************/
    
    
    thincheck = i - std::floor(i/thin) * thin; // % operator stolen by arma
    if(!thincheck )
      Rcpp::Rcout<<"Iteration: "<<i<<" Acceptance lognorm= "<< ((double)acceptance)/i<< 
        " Acceptance alpha1= "<< ((double)acceptance_alph1)/i<< " Acceptance alpha2= "<< ((double)acceptance_alph2)/i<< endl;
    
    // Rcpp::Rcout<<"Iteration: "<<i<<" # Acceptance= "<<acceptance<<" Acceptance prob= "<< ((double)acceptance)/i<<" log_prob="<<log_prob<< endl;
    
    int remainder= (i+1 );
    int quotient= (int) std::floor(remainder/thin);
    remainder-= (quotient*thin) ;
    
    if(remainder==0){
      alloc_var_mat.row(quotient-1)=del.t();
      pimat1.row(quotient-1)=pi1.t();
      pimat2.row(quotient-1)=pi2.t();
      weights2.row(quotient-1)=wght_x2_new.t(); //probability of the cluster
      // weights.row(quotient-1)=wght2.t();
      
      dir_alpha_mat.row(quotient-1)={alpha1,alpha2};
      
      lognormal_mu1.row(quotient-1)=mu1.t(); lognormal_sig1.row(quotient-1)= sig1.t();
      lognormal_mu2.row(quotient-1)=mu2.t(); lognormal_sig2.row(quotient-1)= sig2.t();
      
      hyperparams.row(quotient-1)= current_params.t();
      unifmat.row(quotient-1)=unif.t();
    }
  }
  
  gsl_rng_free (r);
  
  field<arma::mat> lognormal_params1(2), lognormal_params2(2);
  lognormal_params1(0)= lognormal_mu1; lognormal_params1(1)= lognormal_sig1;
  lognormal_params2(0)= lognormal_mu2; lognormal_params2(1)= lognormal_sig2;
  
  return Rcpp::List::create(Rcpp::Named("pimat1") =pimat1,
                            Rcpp::Named("pimat2") =pimat2,
                            // Rcpp::Named("Allocation variables") = alloc_var_mat,
                            Rcpp::Named("Weights2")=weights2,
                            Rcpp::Named("Lognormal_params1")=lognormal_params1,
                            Rcpp::Named("Lognormal_params2")=lognormal_params2,
                            Rcpp::Named("Unifs")=unifmat,
                            Rcpp::Named("Lognormal_hyperparams")=hyperparams,
                            Rcpp::Named("Dirichlet_params")=dir_alpha_mat,
                            Rcpp::Named("Acceptance rates")= 
                              Rcpp::NumericVector::create((double)acceptance_alph1, (double)acceptance_alph2, (double)acceptance )/((double)(nrun+burn) )
  );
}

