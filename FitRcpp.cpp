// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;

// -----------------------------------------------------------------
// ----------------- helper functions ----------------------------------
// [[Rcpp::export]]
Rcpp::NumericMatrix make_mat(Rcpp::List input_list){
  
  unsigned int n = input_list.length();
  
  if(n == 0) {
    Rcpp::stop("Must supply a list with more than 1 element.");
  }
  
  Rcpp::NumericVector testvals = input_list[0];
  unsigned int elems = testvals.length();
  
  Rcpp::NumericMatrix result_mat = Rcpp::no_init(n, elems);
  
  // fill by row
  for(unsigned int i = 0; i < n; i++) {
    NumericVector row_val = input_list[i];
    
    if(elems != row_val.length()) {
      Rcpp::stop("Length of row does not match matrix requirements");
    }
    
    result_mat(i, Rcpp::_) = row_val;
    
  }
  
  return result_mat;
}

// [[Rcpp::export]]
arma::colvec rm_indx(arma::colvec vec, int ind) {
  
  if(ind == vec.n_elem){
    vec = vec(arma::span(0, ind - 2) );
  } else if(ind == 1){
    vec = vec(arma::span(1, vec.n_elem - 1) );
  } else {
    arma::colvec vec1 = vec(arma::span(0, ind - 2) );
    arma::colvec vec2 = vec(arma::span(ind , vec.n_elem - 1 ) );
    
    vec = arma::join_cols(vec1, vec2);
  }
  
  return vec;
}

// [[Rcpp::export]]
List get_allParameters(List dpobj) {
  
  List clustPars = dpobj["clusterParameters"];
  IntegerVector clustLabels = dpobj["clusterLabels"];
  
  int p = clustPars.length();
  int n = dpobj["n"];
  
  List out(p);
  
  
  for(int i = 0; i < p; ++i){
    arma::cube ccCube = clustPars[i];
    arma::cube ncCube =  arma::zeros<arma::cube>(ccCube.n_rows, ccCube.n_cols, n);
    for(int j = 0; j < n; ++j){
      ncCube.slice(j) = ccCube.slice(clustLabels(j) - 1);
    }
    out[i] = ncCube;
  }
  
  return out;
}

// [[Rcpp::export]]
NumericVector combine(const List& list)
{
  std::size_t n = list.size();
  
  // Figure out the length of the output vector
  std::size_t total_length = 0;
  for (std::size_t i = 0; i < n; ++i)
    total_length += Rf_length(list[i]);
  
  // Allocate the vector
  NumericVector output = no_init(total_length);
  
  // Loop and fill
  std::size_t index = 0;
  for (std::size_t i = 0; i < n; ++i)
  {
    NumericVector el = list[i];
    std::copy(el.begin(), el.end(), output.begin() + index);
    
    // Update the index
    index += el.size();
  }
  
  return output;
  
}


// [[Rcpp::export]]
arma::colvec seq_id(int first, int last){
  arma::colvec y(abs(last - first) + 1);
  if(first < last) {
    std::iota(y.begin(), y.end(), first);
  } else {
    std::iota(y.begin(), y.end(), last);
    std::reverse(y.begin(), y.end());
  }
  return y;
}

// [[Rcpp::export]]
int nsample_label(int n, int size, bool replace, arma::colvec probs){
  
  arma::colvec x = seq_id(1, n);
  
  int newLabel = (int)(RcppArmadillo::sample(x, size, replace, probs)[0] );
  
  return newLabel;
}


// [[Rcpp::export]]
double rcpp_sum(arma::vec x){
  double s = std::accumulate(x.begin(), x.end(), 0.0);
  return s;
}

// [[Rcpp::export]]
void rcpp_rprintf(arma::colvec v){
  // printing values of all the elements of Rcpp vector
  for(int i=0; i < v.n_elem; ++i){
    Rprintf("the value of v[%i] : %f \n", i, v[i]);
  }
}

// [[Rcpp::export]]
NumericVector concatenate(NumericVector x, NumericVector y) {
  int nx = x.length();
  int ny = y.length();
  NumericVector out(nx + ny);
  int i;
  
  for(i = 0; i < nx; ++i){
    out[i] = x[i];
  }
  for(i = 0; i < ny; ++i){
    out[i+nx] = y[i];
  }
  return out;
}
/**
* Extended division reminder to vectors
*
* @param a Dividend
* @param n Divisor
*/
// [[Rcpp::export]]
int mod(arma::uvec a, int n)
{
  return (int)(a - arma::floor(a/n)*n)[0];
}

// [[Rcpp::export]]
arma::colvec ewm(arma::colvec x, arma::colvec y){
  return x % y;
}

arma::colvec conc(const arma::colvec & x, const arma::colvec & y) {
  return arma::join_cols(x, y);
}

// [[Rcpp::export]]
arma::colvec erasing(arma::colvec x, int ind){
  
  arma::uvec ids = arma::find(x == ind);
  unsigned int id = ids(0);
  arma::mat X = arma::zeros<arma::mat>(x.n_elem, 1);
  X.col(0) = x;
  X.shed_row(id);
  x = arma::vectorise(X);
  
  return x;
}

// ------------------------------------------------------------
// ------------- Likelihood functions -------------------------

// [[Rcpp::export]]
NumericVector NormalLikelihood(NumericVector x, NumericVector means, NumericVector sds){
  int n = x.size() ;
  NumericVector res(n);
  for(int i = 0; i < n; ++i) res[i] = R::dnorm( x[i], means[i], sds[i], 0);
  return res;
}

//' @title Generate Random Wishart Distribution
//' @description Creates a random wishart distribution when given degrees of freedom and a sigma matrix. 
//' @param df An \code{int}, which gives the degrees of freedom of the Wishart.  (> 0)
//' @param S A \code{matrix} with dimensions m x m that provides Sigma, the covariance matrix. 
//' @return A \code{matrix} that is a Wishart distribution, aka the sample covariance matrix of a Multivariate Normal Distribution
//' @seealso \code{\link{riwishart}} 
//' @author James J Balamuta
//' @examples 
//' #Call with the following data:
//' rwishart(3, diag(2))
//' 
//' # Validation
//' set.seed(1337)
//' S = toeplitz((10:1)/10)
//' n = 10000
//' o = array(dim = c(10,10,n))
//' for(i in 1:n){
//' o[,,i] = rwishart(20, S)
//' }
//' mR = apply(o, 1:2, mean)
//' Va = 20*(S^2 + tcrossprod(diag(S)))
//' vR = apply(o, 1:2, var)
//' stopifnot(all.equal(vR, Va, tolerance = 1/16))
//' 
// [[Rcpp::export]]
arma::mat c_rwishart(unsigned int df, const arma::mat& S){
  // Dimension of returned wishart
  unsigned int m = S.n_rows;
  
  // Z composition:
  // sqrt chisqs on diagonal
  // random normals below diagonal
  // misc above diagonal
  arma::mat Z(m,m);
  
  // Fill the diagonal
  for(unsigned int i = 0; i < m; i++){
    Z(i,i) = sqrt(R::rchisq(df-i));
  }
  
  // Fill the lower matrix with random guesses
  for(unsigned int j = 0; j < m; j++){  
    for(unsigned int i = j+1; i < m; i++){    
      Z(i,j) = R::rnorm(0,1);
    }
  }
  
  // Lower triangle * chol decomp
  arma::mat C = arma::trimatl(Z).t() * arma::chol(S);
  
  // Return random wishart
  return C.t()*C;
}

// [[Rcpp::export]]
arma::colvec rcpp_NormalLikelihood(arma::mat datax, List thetaList){
  
  // get the  parameters for the model
  arma::vec muvec = thetaList[0];
  arma::vec sdvec = thetaList[1];
  // count the current number of clusters
  int nclust = muvec.n_elem;
  
  // create output
  arma::colvec likelihoodValues = arma::zeros<arma::colvec>(nclust);
  
  for(int p = 0; p < nclust; ++p){
    likelihoodValues(p) = (double)(R::dnorm( datax(0,0), muvec(p), sdvec(p), false));
  }
  
  return likelihoodValues;
}

// [[Rcpp::export]]
arma::colvec rcpp_DPNormalLikelihood(List dpobj){
  
  List thetaList = get_allParameters(dpobj);
  arma::mat datax = dpobj["data"];
  
  arma::vec muvec = thetaList[0];
  arma::vec sdvec = thetaList[1];
  int n = dpobj["n"];
  
  arma::colvec likelihoodValues = arma::zeros<arma::colvec>(n);
  
  for(unsigned int idx = 0; idx < n; ++idx){
    likelihoodValues(idx) = (double)(R::dnorm( datax(idx, 0), muvec(idx),  sdvec(idx), false));
  }
  
  return likelihoodValues;
}

// [[Rcpp::export]]
arma::mat rcpp_PosteriorParametersNormal(List mdobj, arma::mat X){
  
  arma::vec priorParameters = mdobj["priorParameters"];
  
  // if Gaussian Mixture Model
  
  arma::vec x = X.col(0);
  int nx = x.n_elem;
  double ybar = rcpp_sum(x)/nx;
  
  double mu0 = priorParameters(0);
  double kappa0 = priorParameters(1);
  double alpha0 = priorParameters(2);
  double beta0 = priorParameters(3);
  
  double mun, kappan, alphan, betan;
  
  mun = (kappa0 * mu0 + nx * ybar)/(kappa0 + nx);
  kappan = kappa0 + nx;
  alphan = alpha0 + nx/2;
  betan = beta0 + 0.5 * rcpp_sum(pow(x - ybar, 2)) + kappa0 * nx * pow(ybar - mu0, 2)/(2 * (kappa0 + nx));
  
  arma::rowvec PPvec = {mun, kappan, alphan, betan};
  arma::mat PPmat = arma::zeros<arma::mat>(1, 4);
  
  PPmat.row(0) = PPvec;
  return PPmat;
}

// [[Rcpp::export]]
List rcpp_PosteriorDrawNormal(List mdobj, arma::mat X, int n = 1) {
  // create PosteriorPars matrix
  arma::mat PosteriorParameters_calc = rcpp_PosteriorParametersNormal(mdobj, X);
  
  // if Gaussian Mixture Model
  
  arma::vec lambda = rgamma(n, PosteriorParameters_calc(0, 2), 1/PosteriorParameters_calc(0, 3));
  arma::vec mu = arma::zeros<arma::vec>(n);
  for(int i = 0; i < n; ++i){
    mu[i] = (double)rnorm(1, PosteriorParameters_calc(0, 0), 1/sqrt(PosteriorParameters_calc(0, 1) * lambda[i])  )[0];
  }
  
  arma::cube muArray = arma::zeros<arma::cube>(1, 1, n);
  arma::cube sdArray = arma::zeros<arma::cube>(1, 1, n);
  
  // update elements
  muArray.slices(0, n - 1) = mu;
  sdArray.slices(0, n - 1) = sqrt(1/lambda);
  
  return  List::create(muArray,
                       sdArray);
}

const double log2pi = std::log(2.0 * M_PI);

//' @title Multivariate Normal Density
// [[Rcpp::export]]
arma::vec dmvnrm_arma(arma::mat x,  
                      arma::rowvec mean,  
                      arma::mat sigma, 
                      bool logd = false) { 
  int n = x.n_rows;
  int xdim = x.n_cols;
  arma::vec out(n);
  arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
  double rootisum = arma::sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0) * log2pi;
  
  for (int i=0; i < n; i++) {
    arma::vec z = rooti * arma::trans( x.row(i) - mean) ;    
    out(i)      = constants - 0.5 * arma::sum(z%z) + rootisum;     
  }  
  
  if (logd == false) {
    out = exp(out);
  }
  return(out);
}

//' @title Generating Multivariate Gaussian
// [[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

// [[Rcpp::export]]
arma::mat rcpp_DPMVNormalLikelihood(List dpobj){
  
  List thetaList = get_allParameters(dpobj);
  arma::mat datax = dpobj["data"];
  
  arma::cube Mu = thetaList[0];
  arma::cube Sigma = thetaList[1];
  int n = dpobj["n"];
  
  arma::mat likelihoodValues = arma::zeros<arma::mat>(n, n);
  
  for(unsigned int id = 0; id < n; ++id){
    likelihoodValues.col(id) = dmvnrm_arma(datax, Mu.slice(id), Sigma.slice(id));
  }
  
  return likelihoodValues;
}

// [[Rcpp::export]]
arma::colvec rcpp_MVNormalLikelihood(arma::mat datax, List thetaList){
  arma::cube Mu = thetaList[0];
  arma::cube Sigma = thetaList[1];
  // count the current number of clusters
  int nclust = Mu.n_slices;
  // create output
  arma::colvec likelihoodValues = arma::zeros<arma::colvec>(nclust);
  
  for(int nc = 0; nc < nclust; ++nc){
    likelihoodValues(nc) = dmvnrm_arma(datax, Mu.slice(nc), Sigma.slice(nc))[0];
  }
  return(likelihoodValues);
}


// [[Rcpp::export]]
List rcpp_PosteriorParametersMVNormal(List mdobj, arma::mat X){
  
  int nx = X.n_rows;
  int np = X.n_cols;
  
  // get Prior parameters
  List priorParameters = mdobj["priorParameters"];
  arma::rowvec mu0 = priorParameters[0];
  arma::mat Lambda = priorParameters[1];
  double kappa0 = priorParameters[2];
  double nu = priorParameters[3];
  
  double kappa_n = kappa0 + nx;
  double nu_n = nu + nx;
  
  // calculate mu_n
  arma::mat mu_n = (kappa0 * mu0 + nx * arma::mean(X)) / (nx + kappa0);
  arma::mat sum_squares = arma::zeros<arma::mat>(np, np);
  // calculate sum of squares
  if(nx > 1){
  sum_squares = (nx - 1) * cov( X ); 
  // replace non_finite values
  sum_squares.replace(arma::datum::nan, 0);
  }
  
  // calculate T_n
  arma::mat t_n = Lambda + sum_squares + ((kappa0 * nx) / (kappa0 + nx) ) * vectorise((mu0 - arma::mean(X))) * vectorise((mu0 - arma::mean(X))).t();
  
  return List::create( _["mu_n"] = mu_n,
                       _["t_n"] = t_n,
                       _["kappa_n"] = kappa_n,
                       _["nu_n"] = nu_n);
}


// [[Rcpp::export]]
List rcpp_PosteriorDrawMVNormal(List mdobj, arma::mat X, int n = 1) {
  
  // create PosteriorPars matrix
  List PosteriorParameters_calc = rcpp_PosteriorParametersMVNormal(mdobj, X);
  
  // declare input
  double nu_n =  PosteriorParameters_calc["nu_n"];
  arma::mat t_n = PosteriorParameters_calc["t_n"]; 
  unsigned int p = t_n.n_rows;
  arma::vec mu_n = PosteriorParameters_calc["mu_n"];
  double kappa_n = PosteriorParameters_calc["kappa_n"];
  
  
  // declare output
  arma::cube sig(p, p, n);
  arma::cube mu(1, p, n);
  //arma::mat Sigma(p, p);
  
  for(int i = 0; i < n; ++i ){
    sig.slice(i) =  arma::wishrnd(t_n, nu_n);
    mu.slice(i) = arma::trans(arma::mvnrnd(mu_n, arma::inv(kappa_n * sig.slice(i) ), 1 ) );
  }
  
  sig = sig / pow( kappa_n, 2);
  
  return List::create(_["mu"] = mu,
                      _["sig"] = sig);
}


// ------------------------------------------------------------
// ------------------- Cluster Label Change -------------------
// [[Rcpp::export]]
List rcpp_ClusterLabelChange(List dpObj, int i, int newLabel, int currentLabel) {

  arma::mat Y = dpObj["data"];
  arma::mat X = Y.row(i);
  int numLabs = dpObj["numberClusters"];
  List mdObj = dpObj["mixingDistribution"];
  
  arma::colvec ppclust = dpObj["pointsPerCluster"];
  IntegerVector clustLabs = dpObj["clusterLabels"];
  List clustPars = dpObj["clusterParameters"];
  int n = dpObj["n"];
  StringVector dist = mdObj["distribution"];

  if (newLabel <= numLabs) {
    //Rcout << "taking this route 1 \n";

    ++ppclust[newLabel - 1];
    clustLabs[i] = newLabel;

    if (ppclust[currentLabel - 1] == 0){
      //Rcout << "way 1 and 2\n";
      // Remove Empty Cluster
      numLabs = numLabs - 1;

      ppclust = rm_indx(ppclust, currentLabel);

      // Number of elements in the List input
      for (int j = 0; j < clustPars.length(); ++j) {
        arma::cube cubeArray = clustPars[j];
        cubeArray.shed_slice(currentLabel - 1);
        clustPars[j] = cubeArray;
      }

      LogicalVector inds = clustLabs > currentLabel;
      for(int j = 0; j < n; ++j){
        if(inds[j]){
          clustLabs[j] = clustLabs[j] - 1;
        }
      }
    }
  } else {

    if(ppclust[currentLabel - 1] == 0) {

      //Rcout << "way 3\n";
      List post_draw;
      if(dist[0] == "normal"){
        post_draw = rcpp_PosteriorDrawNormal(mdObj, X);
      } else if(dist[0] == "mvnormal") {
        post_draw = rcpp_PosteriorDrawMVNormal(mdObj, X);
        } else {
        // printing error message
        Rcerr << "The distribution is not suported\n";
      }
      for(int j = 0; j < clustPars.length(); ++j){
        arma::cube cubeArray = clustPars[j];
        arma::cube par_post_draw = post_draw[j];
        cubeArray.slice(currentLabel - 1) = par_post_draw;
        clustPars[j] = cubeArray;
      }

      ppclust[currentLabel - 1] = ppclust[currentLabel - 1] + 1;

    } else {

      //Rcout << "way 4\n";

      clustLabs[i] = newLabel;
      numLabs = numLabs + 1;
      ppclust = arma::join_cols(ppclust, arma::ones<arma::colvec>(1));
      
      List post_draw;

      if(dist[0] == "normal"){
        post_draw = rcpp_PosteriorDrawNormal(mdObj, X);
      } else if(dist[0] == "mvnormal") {
        post_draw = rcpp_PosteriorDrawMVNormal(mdObj, X);
      } else {
        // printing error message
        Rcerr << "The distribution is not suported\n";
      }

      int np = clustPars.length();

      for(int j = 0; j < np; ++j){
        NumericVector parVec = clustPars[j];
        NumericVector postVec = post_draw[j];
        arma::cube parCube = clustPars[j];
        arma::cube postCube = post_draw[j];

        NumericVector vecArray = concatenate(parVec, postVec);

        arma::cube cubeArray(vecArray.begin(), postCube.n_rows,
                             postCube.n_cols, parCube.n_slices + 1);
        clustPars[j] = cubeArray;
      }
    }
  }

  dpObj["pointsPerCluster"] = ppclust;
  dpObj["clusterLabels"] = clustLabs;
  dpObj["clusterParameters"] = clustPars;
  dpObj["numberClusters"] = numLabs;
  return dpObj;
}


// [[Rcpp::export]]
List rcpp_clusterComponenetUpdate(List dpobj) {

  arma::mat Y = dpobj["data"];
  int n = dpobj["n"];
  double alpha = dpobj["alpha"];
  IntegerVector clustLabs = dpobj["clusterLabels"];
  int newLabel, currentLabs, nclust;
  List clustPars = dpobj["clusterParameters"];
  int numLabs = dpobj["numberClusters"];
  List mdObj = dpobj["mixingDistribution"];
  StringVector dist = mdObj["distribution"];
  arma::colvec ppclust = dpobj["pointsPerCluster"];
  arma::colvec predArray = dpobj["predictiveArray"];

  for(int i = 0; i < n; ++i){

    currentLabs = clustLabs[i];

    --ppclust[currentLabs - 1];
    nclust = ppclust.n_elem;

    arma::colvec probs = arma::zeros<arma::colvec>(nclust + 1);

    // element-wise multiplication
    arma::colvec pp = arma::zeros<arma::colvec>(1);
    pp(0) = alpha * predArray(i);

    
    if(dist[0] == "normal"){ //normal distribution
      probs = arma::join_cols( ppclust % rcpp_NormalLikelihood(Y.row(i), clustPars) , pp);
    } else if(dist[0] == "mvnormal") { 
      probs = arma::join_cols( ppclust % rcpp_MVNormalLikelihood(Y.row(i), clustPars) , pp);
    } else {
      // printing error message
      Rcerr << "The distribution is not suported\n";
    }
    //Rcout << "The value of probs : " << probs << "\n";
    
    probs.elem(arma::find_nonfinite(probs) ).zeros();

    if(bool status1 = arma::all(probs == 0)){
      probs.ones();
    }

    newLabel = nsample_label(numLabs + 1, 1, false, probs);

    dpobj["pointsPerCluster"] = ppclust;

    List dpobj_new = rcpp_ClusterLabelChange(dpobj, i, newLabel, currentLabs);

    arma::colvec ppclust_new = dpobj_new["pointsPerCluster"];
    IntegerVector clustLabs_new = dpobj_new["clusterLabels"];
    List clustPars_new = dpobj_new["clusterParameters"];
    int numLabs_new = dpobj_new["numberClusters"];

    ppclust = ppclust_new;
    clustLabs = clustLabs_new;
    clustPars = clustPars_new;
    numLabs = numLabs_new;

  }

  dpobj["pointsPerCluster"] = ppclust;
  dpobj["clusterLabels"] = clustLabs;
  dpobj["clusterParameters"] = clustPars;
  dpobj["numberClusters"] = numLabs;

  return dpobj;
}


// ------ Cluster Parameter Update ---------------
// -----------------------------------------------
// [[Rcpp::export]]
List rcpp_clusterParameterUpdate(List dpObj){
  List dpobj = clone(dpObj);

  arma::mat Y = dpobj["data"];
  int numLabs = dpobj["numberClusters"];

  arma::uvec clustLabs = dpobj["clusterLabels"];
  //IntegerVector ppclust = dpobj["pointsPerCluster"];
  //int newLabel, currentLabs, nclust;
  List clustPars = dpobj["clusterParameters"];
  int npars = clustPars.length();
  List mdobj = dpobj["mixingDistribution"];
  StringVector dist = mdobj["distribution"];

  for(int i = 0; i < numLabs; ++i){
    arma::uvec ids = arma::find( clustLabs == (i + 1) );
    arma::mat pts = Y.rows( ids );
    List post_draw;
    if(dist[0] == "normal"){
      post_draw = rcpp_PosteriorDrawNormal(mdobj, pts);
    } else if(dist[0] == "mvnormal") {
      post_draw = rcpp_PosteriorDrawMVNormal(mdobj, pts);
    } else {
      // printing error message
      Rcerr << "The distribution is not suported\n";
    }

    for(int j = 0; j < npars; ++ j) {
      arma::cube cubeArray = clustPars[j];
      arma::cube cubePostdraw = post_draw[j];
      cubeArray.slice(i) = cubePostdraw;
      clustPars[j] = cubeArray;
    }

  }

  dpobj["clusterParameters"] = clustPars;

  return(dpobj);
}

// ---------------------------------
// ---------- Update alpha ---------
// [[Rcpp::export]]
double update_concentration(double oldParam, int n, int nParams, arma::colvec priorParameters){

  double newalpha = 0;

  double x = (double)Rcpp::rbeta(1, oldParam + 1, n)[0];

  double pi1 = priorParameters[0] + nParams - 1;

  double pi2 = n * (priorParameters[1] - log(x));

  pi1 = pi1/(pi1 + pi2);

  bool dice = (double)Rcpp::runif(1)[0] < pi1;

  if(dice) {
    newalpha =  (double)rgamma(1, priorParameters[0] + nParams, 1/ (priorParameters[1] - log(x)) )[0];
  } else {
    newalpha =  (double)rgamma(1, priorParameters[0] + nParams - 1, 1/ (priorParameters[1] - log(x)) )[0];
  }
  return newalpha;
}

// [[Rcpp::export]]
List update_alpha(List dpobj){

  double alpha = (double)dpobj["alpha"];
  int n = (int)dpobj["n"];
  int nClust = (int)dpobj["numberClusters"];
  arma::colvec priorParams = dpobj["alphaPriorParameters"];

  double newAlpha = update_concentration(alpha , n, nClust, priorParams );
  dpobj["alpha"] = newAlpha;

  return dpobj;

}


// -------- Fit function -------------------------
// -----------------------------------------------

// [[Rcpp::export]]
List rcpp_Fit(List dpobj, unsigned int its, bool updatePrior = false, bool progressBar = false) {
  List dpObj = clone(dpobj);
  arma::colvec alphaChain(its);
  List weightsChain = Rcpp::List(its);
  List clusterParametersChain = Rcpp::List(its);
  List priorParametersChain = Rcpp::List(its);
  List labelsChain = Rcpp::List(its);
  arma::colvec likelihoodChain(its);
  arma::mat datax = dpObj["data"];

  for(int i = 0; i < its; ++i){

    alphaChain[i] = (double)dpObj["alpha"];
    arma::colvec ppclust = dpObj["pointsPerCluster"];
    weightsChain[i] = ppclust / double(dpObj["n"]);
    List clustParmas = dpObj["clusterParameters"];
    clusterParametersChain[i] = clustParmas;
    List mdObj = dpObj["mixingDistribution"];
    StringVector dist = mdObj["distribution"];
    List priorParams = mdObj["priorParameters"];
    priorParametersChain[i] = priorParams;
    arma::colvec clustLabels = dpObj["clusterLabels"];
    labelsChain[i] = clustLabels;
    

    if(dist[0] == "normal"){
      likelihoodChain[i] =  rcpp_sum(log(rcpp_DPNormalLikelihood(dpObj) ) );
    } else if(dist[0] == "mvnormal") {
      likelihoodChain[i] =  rcpp_sum(log(arma::vectorise(rcpp_DPMVNormalLikelihood(dpObj) ) ) );
    } else {
      // printing error message
      Rcerr << "The distribution is not suported\n";
    }
    

    dpObj = rcpp_clusterComponenetUpdate(dpObj);
    dpObj = rcpp_clusterParameterUpdate(dpObj);
    dpObj = update_alpha(dpObj);
    

    if(updatePrior){
      //mdObj = prior_parameters_update(mdObj, arma::colvec dpObj["clusterLabels"])
    }
    if(progressBar){
      //set_txt_progress_bar(pb, i);
    }

  }
  arma::colvec ppclust = dpObj["pointsPerCluster"];
  dpObj["weight"] = ppclust / double(dpObj["n"]);
  dpObj["alphaChain"] = alphaChain;
  dpObj["likelihoodChain"] = likelihoodChain;
  dpObj["weightChain"] = weightsChain;
  dpObj["clusterParametersChain"] = clusterParametersChain;
  dpObj["priorParametersChain"] = priorParametersChain;
  dpObj["labelsChain"] = labelsChain;

  if(progressBar){
    //close(pb);
  }

  return(dpObj);
}


/*** R
dp = dirichletprocess::DirichletProcessMvnormal(scale(faithful))
#dp = dirichletprocess::DirichletProcessGaussian(rt(200, 3) + 3)
dp2 = rcpp_Fit(dp, 20)
#test_out <- rcpparma_outerproduct(c(1:4))
*/