///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// Henry Scharf
//
// This c++ code fits the model to data.
//
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

#include "RcppArmadillo.h"
#include <RcppArmadilloExtensions/sample.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

///////////////////
// BEGIN HELPERS //
///////////////////
double CalcNormDensity(const double x,
                       const double mu = 0,
                       const double sigma_sq = 1,
                       const int ln = 1) {
  double out = -0.5*log(2*datum::pi*sigma_sq) - 0.5*pow(x - mu, 2) / sigma_sq;
  if (ln == 0)
    out = exp(out);
  return out;
}

arma::cube ReadCube(const NumericVector myArray){
  NumericVector vecArray(myArray);
  IntegerVector arrayDims = vecArray.attr("dim");
  if (arrayDims.size() !=3)
    throw std::domain_error("array not 3D");
  arma::cube cubeArray(vecArray.begin(), arrayDims(0), arrayDims(1), arrayDims(2));
  return cubeArray;
}

field<cube> ReadField(const NumericVector myArray){
  NumericVector vecArray(myArray);
  IntegerVector arrayDims = vecArray.attr("dim");
  if (arrayDims.size() != 4)
    throw std::domain_error("4D only, please");
  field<cube> fieldArray(arrayDims(3));
  for (int d = 0; d<arrayDims(3); d++) {
    arma::cube cube_tmp(vecArray.begin() + d * arrayDims(0) * arrayDims(1) * arrayDims(2),
                        arrayDims(0), arrayDims(1), arrayDims(2));
    fieldArray(d) = cube_tmp;
  }
  return fieldArray;
}

/////////////////
// END HELPERS //
/////////////////

/////////////////////
// BEGIN FUNCTIONS //
/////////////////////

arma::cube MakeMubar(const arma::cube w_cube,
                     const arma::cube mu_cube) {
  //
  // Args:
  //    w: 3D array (n_animals, n_animals, T)
  //    mu: 3D array (n_animals, 2 (lat, long), T)
  //
  // Return:
  //    3D array same dim as mu
  //
  arma::vec one_vector(w_cube.n_rows, fill::ones);
  arma::mat identity_n_animals = diagmat(one_vector);
  int T = w_cube.n_slices;
  if (T != mu_cube.n_cols)
    throw std::domain_error("time dimensions do not match");
  arma::cube out(mu_cube.n_rows, mu_cube.n_cols, mu_cube.n_slices);
  for (int t = 0; t < T; t++) {
    arma::mat w_t = w_cube.slice(t);
    arma::mat mu_t = mu_cube(0, t, 0, size(mu_cube.n_rows, 1, mu_cube.n_slices));
    arma::mat w_t_rowsums(mu_cube.n_rows, mu_cube.n_slices, fill::zeros);
    for (int d = 0; d < mu_cube.n_slices; d++) {
      w_t_rowsums.col(d) = sum((identity_n_animals + w_cube.slice(t)), 1);
    }
    out(0, t, 0, size(mu_cube.n_rows, 1, mu_cube.n_slices)) =
      (identity_n_animals + w_t) * mu_t / w_t_rowsums;
  }
  return out;
}

arma::cube MakeMuTilde(const arma::cube w_cube,
                       const arma::cube mu_cube,
                       const arma::cube mubar_cube) {
  //
  // Args:
  //    w: cube (n_animals, n_animals, T)
  //    mu: cube (n_animals, 2 (lat, long), T)
  //
  // Return:
  //    cube same dim as mu
  //
  arma::cube difference = mubar_cube - mu_cube;
  arma::mat norm(mu_cube.n_rows, mu_cube.n_cols, fill::zeros);
  for (int d = 0; d < mu_cube.n_slices; d++) {
    norm = norm + pow(difference.slice(d), 2);
  }
  for (int row = 0; row < norm.n_rows; row++) {
    for (int col = 0; col < norm.n_cols; col++) {
      if (norm(row, col) == 0) {
        norm(row, col) = 1;
      }
    }
  }
  norm = sqrt(norm);
  arma::cube norm_cube(difference.n_rows, difference.n_cols, mu_cube.n_slices);
  for (int d = 0; d < mu_cube.n_slices; d++) {
    norm_cube.slice(d) = norm;
  }
  arma::cube out = difference / norm_cube;
  return out;
}

mat ModMuTilde(const arma::mat w_t,
               const arma::mat mu_t,
               const int i,
               const int j,
               const int z) {
  //
  // Args:
  //    w_t:  matrix of dim (n_animals, n_animals)
  //    mu_t: matrix of dim (n_animals, 2)
  //
  // Return: matrix (2, n_animals) which is mu_tilde_t with the modification
  //    based on z.
  //
  int n_animals = mu_t.n_rows;
  mat w_t_mod = w_t;
  w_t_mod(i - 1, j - 1) = z;
  w_t_mod(j - 1, i - 1) = z;
  rowvec rowsums;
  mat w_t_mod_rowsums(n_animals, mu_t.n_cols);
  for (int d = 0; d < mu_t.n_cols; d++) {
    w_t_mod_rowsums.col(d) = sum(w_t_mod, 1) + 1;
  }
  mat identity_1(n_animals, mu_t.n_cols, fill::ones);
  w_t_mod_rowsums = max(w_t_mod_rowsums, identity_1);
  mat identity_diag(n_animals, n_animals, fill::eye);
  mat mubar_mod_t = (w_t_mod + identity_diag) * mu_t / w_t_mod_rowsums;
  mat difference_t = mubar_mod_t - mu_t;
  mat norm(n_animals, mu_t.n_cols);
  for (int d = 0; d < mu_t.n_cols; d++) {
    for (int ii = 0; ii < n_animals; ii++) {
      norm(ii, d) = pow(pow(difference_t(ii, 0), 2) + pow(difference_t(ii, 1), 2), 0.5);
      if (norm(ii, d) == 0) {
        norm(ii, d) = 1;
      }
    }
  }
  mat mutil_mod = difference_t/norm;
  return(mutil_mod);
}

arma::cube MakeHbar(const arma::cube w_cube,
                    const arma::cube mu_cube,
                    const arma::cube mutil_cube,
                    const double beta_iter) {
  //
  // Args:
  //    w_cube: cube of dim (n_animals, n_animals, T)
  //    mu_cube: cube of dim (n_animals, 2 (lat, long), T)
  //    mutil_cube: cube of dim (n_animals, 2 (lat, long), T)
  //    iter: iteration in the MCMC chain min value is 1 NOT 0
  //
  // Return:
  //    cube of dim (n_animals, T - 1, 2)
  //
  int T = w_cube.n_slices;
  arma::cube h = mu_cube(0, 1, 0, size(mu_cube.n_rows, T - 1, mu_cube.n_slices)) -
    mu_cube(0, 0, 0, size(mu_cube.n_rows, T - 1, mu_cube.n_slices)) -
    beta_iter * mutil_cube(0, 0, 0, size(mutil_cube.n_rows, T - 1, mutil_cube.n_slices));
  // time index for h_t is now: 0 indexes t = 2
  arma::cube hbar(mu_cube.n_rows, T - 1, mu_cube.n_slices);
  for (int t = 0; t < (T - 1); t++) {
    arma::mat h_t = h(0, t, 0, size(h.n_rows, 1, h.n_slices));
    arma::mat w_t = w_cube.slice(t + 1);
    arma::mat w_t_rowsums(mu_cube.n_rows, mu_cube.n_slices, fill::zeros);
    for (int d = 0; d < mu_cube.n_slices; d++) {
      w_t_rowsums.col(d) = sum(w_t, 1);
    }
    arma::mat matrix_of_ones(mu_cube.n_rows, mu_cube.n_slices, fill::ones);
    w_t_rowsums = max(w_t_rowsums, matrix_of_ones);
    hbar(0, t, 0, size(h.n_rows, 1, h.n_slices)) = w_t * h_t / w_t_rowsums;
  }
  return hbar;
}

arma::mat ModHbarATt(const arma::mat w_t,
		     const arma::mat h_t,
		     const int i,
		     const int j,
		     const int z) {
  //
  // Args:
  //    w_t: matrix of dim (n_animals, n_animals)
  //    h_t: matrix of dim (n_animals, 2)
  //    i: min value is 1 NOT 0
  //    j: min value is 1 NOT z
  //    z: 0 or 1
  //
  // Return:
  //    matrix of dim (n_animals, 2)
  //
  int n_animals = w_t.n_rows;
  arma::mat w_t_mod = w_t;
  w_t_mod(i - 1, j - 1) = z;
  arma::mat hbar_mod(n_animals, 2);
  arma::mat w_t_rowsums(h_t.n_rows, h_t.n_cols, fill::zeros);
  for (int d = 0; d < h_t.n_cols; d++) {
    w_t_rowsums.col(d) = sum(w_t_mod, 1);
  }
  arma::mat matrix_of_ones(h_t.n_rows, h_t.n_cols, fill::ones);
  w_t_rowsums = max(w_t_rowsums, matrix_of_ones);
  hbar_mod = w_t_mod * h_t / w_t_rowsums;
  return hbar_mod;
}


// double MakeNormDensityAT1(const int i,
// 			  const int j,
// 			  const int z,
// 			  const arma::cube w_cube,
// 			  const arma::cube mu_cube,
// 			  const double alpha_iter,
// 			  const double sigsq_iter) {
//   //
//   // WARNING: THIS IS DEPRICATED AND NOT IN USE ANYMORE
//   //
//   // Args:
//   //    i: min value is 1 NOT 0
//   //    j: min value is 1 NOT 0
//   //    z: 0 or 1
//   //    w_cube: cube of dim (n_animals, n_animals, T)
//   //    mu_cube: cube of dim (n_animals, 2 (lat, long), T)
//   //    mubar_cube: cube of dim (n_animals, 2 (lat, long), T). This is the
//   //                centroid location.
//   //    alpha: vector of length n_iterations
//   //    sigsq: vector of length n_iterations
//   //    iter: iteration in the MCMC chain min value is 1 NOT 0
//   //
//   arma::cube w_cube_function = w_cube;
//   w_cube_function(i - 1, j - 1, 0) = z;
//   arma::rowvec gauss_densities(mu_cube.n_slices);
//   arma::cube mubar_cube = MakeMubar(w_cube_function, mu_cube);
//   for (int d = 0; d < mu_cube.n_slices; d++) {
//     double mean = alpha_iter * mubar_cube(i - 1, 0, d);
//     arma::rowvec w_it = w_cube(i - 1, 0, 0, size(1, w_cube.n_cols, 1));
//     int w_it_sum = sum(w_it) - w_it(j - 1) + z;
//     if (w_it_sum == 0)
//       w_it_sum = c_iter;
//     double sd = sqrt(sigsq_iter / w_it_sum);
//     gauss_densities(d) = R::dnorm(mu_cube(i - 1, 0, d), mean, sd, 0);
//     /* boost::math::normal gauss(mean, sd); */
//     /* gauss_densities(d) = pdf(gauss, mu_cube(i - 1, 0, d)); */
//   }
//   return prod(gauss_densities);
// }

double MakeNormDensityATt(const int i,
			  const int j,
			  const int t,
			  const int z,
			  const arma::cube w_cube,
			  const arma::cube mu_cube,
			  const arma::cube mutil_cube,
			  const double alpha_iter,
			  const double beta_iter,
			  const double sigsq_iter,
			  const double c_iter) {
  //
  // Args:
  //    i: min value is 1 NOT 0
  //    j: min value is 1 NOT 0
  //    t: min value is 1 NOT 0
  //    z: 0 or 1
  //    w_cube: cube of dim (n_animals, n_animals, T)
  //    mu_cube: cube of dim (n_animals, T, 2 (lat, long))
  //    mutil_cube: cube of dim (n_animals, T, 2 (lat, long)). This is a vector
  //                pointing to the centroid location from the previous location
  //                at each time point.
  //    alpha_iter: nth_iteration of alpha
  //    beta_iter:
  //    sigsq_iter:
  //    c_iter:
  //
  arma::mat h_t =
    mu_cube(0, t - 1, 0, size(mu_cube.n_rows, 1, mu_cube.n_slices)) -
    mu_cube(0, t - 2, 0, size(mu_cube.n_rows, 1, mu_cube.n_slices)) -
    beta_iter * mutil_cube(0, t - 2, 0, size(mutil_cube.n_rows, 1, mutil_cube.n_slices));
  arma::mat hbar_mod =
    ModHbarATt(w_cube.slice(t - 1), h_t, i, j, z);
  arma::rowvec gauss_densities(mu_cube.n_slices);
  arma::vec mean(2, fill::zeros);
  arma::vec sigma_sq(2, fill::ones);
  arma::vec sum_w_it(2, fill::zeros);
  for (int d = 0; d < mu_cube.n_slices; d++) {
    mean(d) = mu_cube(i - 1, t - 2, d) +
      beta_iter * mutil_cube(i - 1, t - 2, d) +
      alpha_iter * hbar_mod(i - 1, d);
    arma::rowvec w_it = w_cube(i - 1, 0, t - 1, size(1, w_cube.n_cols, 1));
    sum_w_it(d) = as_scalar(sum(w_it) + (z - w_it(j - 1)));
    if (sum_w_it(d) == 0)
      sum_w_it(d) = c_iter;
    sigma_sq(d) = sigsq_iter / sum_w_it(d);
    gauss_densities(d) = CalcNormDensity(mu_cube(i - 1, t - 1, d), mean(d), sigma_sq(d), 0);
  }
  return prod(gauss_densities);
}

double MakeNormDensityATtp1(const int i,
			    const int j,
			    const int t,
			    const int z,
			    const arma::cube w_cube,
			    const arma::cube mu_cube,
			    const arma::cube mutil_cube,
			    const double alpha_iter,
			    const double beta_iter,
			    const double sigsq_iter,
			    const double c_iter) {
  //
  // Args:
  //    i: min value is 1 NOT 0
  //    j: min value is 1 NOT 0
  //    t: min value is 1 NOT 0
  //    z: 0 or 1
  //    w_cube: cube of dim (n_animals, n_animals, T)
  //    mu_cube: cube of dim (n_animals, T, 2 (lat, long))
  //    mutil_cube: cube of dim (n_animals, T, 2 (lat, long)). This is vector
  //                pointing to the centroid location from the previous
  //                location.
  //
  arma::mat mutil_mod =
    ModMuTilde(w_cube.slice(t - 1),
  	       mu_cube(0, t - 1, 0, size(mu_cube.n_rows, 1, mu_cube.n_slices)),
  	       i, j, z);
  arma::cube mutil_mod_cube = mutil_cube;
  mutil_mod_cube(0, t - 1, 0, size(mutil_cube.n_rows, 1, mutil_cube.n_slices)) =
    mutil_mod;
  arma::cube hbar = MakeHbar(w_cube, mu_cube, mutil_mod_cube, beta_iter);
  arma::rowvec gauss_densities(mu_cube.n_slices);
  arma::vec mean(2, fill::zeros);
  arma::vec sigma_sq(2, fill::ones);
  arma::vec sum_w_itp1(2, fill::zeros);
  for (int d = 0; d < mu_cube.n_slices; d++) {
    mean(d) = mu_cube(i - 1, t - 1, d) + beta_iter * mutil_mod(i - 1, d) +
      alpha_iter * hbar(i - 1, t - 1, d);
    arma::rowvec w_itp1 = w_cube(i - 1, 0, t, size(1, w_cube.n_cols, 1));
    sum_w_itp1(d) = as_scalar(sum(w_itp1));
    if (sum_w_itp1(d) == 0)
      sum_w_itp1(d) = c_iter;
    sigma_sq(d) = sigsq_iter / sum_w_itp1(d);
    gauss_densities(d) = CalcNormDensity(mu_cube(i - 1, t, d),
					 mean(d), sigma_sq(d), 0);
  }
  return prod(gauss_densities);
}

arma::mat MakePrecisionKernelATt(const int t,
				 const double alpha_iter,
				 const arma::cube w_cube,
				 const double c = 0.5) {
  //
  // Args:
  //    t: min value is 1 NOT 0
  //    alpha: vector of length n_iterations
  //    w_cube: cube of dim (n_animals, n_animals, T)
  //
  // Return:
  //    matrix of dimension (2*n_animals, 2*n_animals)
  //
  arma::vec one_vector(w_cube.n_rows, fill::ones);
  arma::vec one_vector_2(2, fill::ones);
  arma::mat identity_n_animals = diagmat(one_vector);
  arma::mat identity_2 = diagmat(one_vector_2);
  arma::mat w_t = w_cube.slice(t - 1);
  arma::vec w_t_rowsums = max(c*one_vector, sum(w_t, 1));
  arma::mat sub_Kernel_t = diagmat(w_t_rowsums) - alpha_iter * w_t;
  arma::mat Kernel_t = kron(sub_Kernel_t, identity_2);
  return Kernel_t;
}

arma::cube MakePrecisionKernel(const double alpha_iter,
			       const arma::cube w_cube,
			       const double c = 0.5) {
  //
  // Args:
  //    alpha: vector of length n_iterations
  //    w_cube: cube of dim (n_animals, n_animals, T)
  //
  // Return:
  //    cube of dimension (n_animals, n_animals, T)
  //
  int T = w_cube.n_slices;
  arma::cube Kernel(2 * w_cube.n_rows, 2 * w_cube.n_cols, T);
  for (int t = 0; t < T; t++) {
    Kernel.slice(t) = MakePrecisionKernelATt(t + 1, alpha_iter, w_cube, c);
  }
  return Kernel;
}

arma::mat MakeConditionalMMatrix(const double p1,
                                 const double phi_iter) {
  // Args:
  //
  // Return: matrix with entries:
  //         (0, 0) = p0g0
  //         (0, 1) = p0g1
  //         (1, 0) = p1g0
  //         (1, 1) = p1g1
  //
  arma::mat M(2, 2);
  M(0, 0) = 1 - (1 - phi_iter) * p1;
  M(1, 0) = (1 - phi_iter) * p1;
  M(0, 1) = (1 - phi_iter) * (1 - p1);
  M(1, 1) = 1 - (1 - phi_iter) * (1 - p1);
  return(M);
}

double CalcLogPhiDensity(const double phi,
			 const double p1,
			 const double alpha_phi,
			 const double beta_phi,
			 const arma::cube w_cube) {
  //
  // Args:
  //
  // Return:
  //
  int T = w_cube.n_slices;
  int n_animals = w_cube.n_rows;
  arma::mat conditional_M_mat = MakeConditionalMMatrix(p1, phi);
  arma::vec one_vector(n_animals, fill::ones);
  arma::mat identity_n_animals = diagmat(one_vector);
  arma::vec g11(T - 1);
  arma::vec g10(T - 1);
  arma::vec g01(T - 1);
  arma::vec g00(T - 1);
  for (int t = 0; t < T - 1; t++) {
    arma::vec w_t = vectorise(trimatu(w_cube.slice(t)));
    arma::vec one_minus_w_t =
      vectorise(trimatu(1 - w_cube.slice(t) - identity_n_animals));
    arma::vec w_t_plus_1 = vectorise(trimatu(w_cube.slice(t + 1)));
    arma::vec one_minus_w_t_plus_1 =
      vectorise(trimatu(1 - w_cube.slice(t + 1) - identity_n_animals));
    g11(t) = as_scalar(w_t.t() * w_t_plus_1);
    g01(t) = as_scalar(w_t.t() * one_minus_w_t_plus_1);
    g10(t) = as_scalar(one_minus_w_t.t() * w_t_plus_1);
    g00(t) = as_scalar(one_minus_w_t.t() * one_minus_w_t_plus_1);
  }
  double log_phi_density =
    sum(g11) * log(conditional_M_mat(1, 1)) +
    sum(g01) * log(conditional_M_mat(0, 1)) +
    sum(g10) * log(conditional_M_mat(1, 0)) +
    sum(g00) * log(conditional_M_mat(0, 0)) +
    (alpha_phi - 1) * log(phi) + (beta_phi - 1) * log(1 - phi);
  return as_scalar(log_phi_density);
}

double CalcLogAlphaDensity(const double alpha_iter,
			   const double alpha_alpha,
			   const double beta_alpha,
			   const arma::cube w_cube,
			   const arma::cube mu_cube,
			   const arma::cube mutil_cube,
			   const arma::cube K_cube,
			   const double beta_iter,
			   const double sigsq_iter) {
  //
  //  Args:
  //
  //
  double log_alpha_density = -1E8;
  if (alpha_iter >= -1 && alpha_iter <= 1) {
    int n_animals = w_cube.n_rows;
    int T = w_cube.n_slices;
    int d = mu_cube.n_slices;
    arma::cube hbar_cube = MakeHbar(w_cube, mu_cube, mutil_cube, beta_iter);
    arma::cube h_cube = mu_cube(0, 1, 0, size(n_animals, T - 1, d)) -
      mu_cube(0, 0, 0, size(n_animals, T - 1, d)) -
      beta_iter * mutil_cube(0, 0, 0, size(n_animals, T - 1, d));
    double alpha_sum = 0;
    for (int i = 0; i < n_animals; i++) {
      for (int t = 0; t < (T - 1); t++) {
	arma::vec h_hbar_it = h_cube(i, t, 0, size(1, 1, d)) -
	  alpha_iter * hbar_cube(i, t, 0, size(1, 1, d));
	double h_hbar_inner_prod = as_scalar(h_hbar_it.t() * h_hbar_it);
	double new_term = K_cube(2*i, 2*i, t + 1) * h_hbar_inner_prod;
	alpha_sum = alpha_sum + new_term;
      }
    }
    log_alpha_density = - alpha_sum / (2 * sigsq_iter) +
      (alpha_alpha - 1) * log(1 + alpha_iter) +
      (beta_alpha - 1) * log(1 - alpha_iter);
  }
  return log_alpha_density;
}

double CalcLogP1Density(const double p1,
			const double phi,
			const arma::cube w_cube,
			const double alpha_p1,
			const double beta_p1) {
  int n_animals = w_cube.n_rows;
  int T = w_cube.n_slices;
  arma::cube omwtmo_wt =
    (1 - w_cube(0, 0, 0, size(n_animals, n_animals, T - 1))) %
    w_cube(0, 0, 1, size(n_animals, n_animals, T - 1));
  arma::cube wtmo_omwt =
    w_cube(0, 0, 0, size(n_animals, n_animals, T - 1)) %
    (1 - w_cube(0, 0, 1, size(n_animals, n_animals, T - 1)));
  arma::cube omwtmo_omwt =
    (1 - w_cube(0, 0, 0, size(n_animals, n_animals, T - 1))) %
    (1 - w_cube(0, 0, 1, size(n_animals, n_animals, T - 1)));
  arma::cube wtmo_wt =
    w_cube(0, 0, 0, size(n_animals, n_animals, T - 1)) %
    w_cube(0, 0, 1, size(n_animals, n_animals, T - 1));
  double sum_p1 = 0.5 * (accu(w_cube.slice(0)) + accu(omwtmo_wt)) +
    alpha_p1 - 1;
  double sum_omp1 = 0.5 * (accu(1 - w_cube.slice(0)) + accu(wtmo_omwt)) +
    beta_p1 - 1;
  double sum_p0g0 = 0.5 * accu(omwtmo_omwt);
  double sum_p1g1 = 0.5 * accu(wtmo_wt);
  double log_density = log(p1) * sum_p1 + log(1 - p1) * sum_omp1 +
    log(1 - (1 - phi)*p1) * sum_p0g0 + log(1 - (1 - phi)*(1 - p1)) * sum_p1g1;
  return log_density;
}

double SamplePhi(const arma::cube w_cube,
		 const double beta_phi_tune,
		 const double p1,
		 const double alpha_p,
		 const double beta_p,
		 const double phi_iter) {
  //
  // Args:
  //
  //
  NumericVector phi_options(2);
  /* previous value */
  phi_options(0) = phi_iter;
  /* proposed value */
  phi_options(1) = R::rbeta(beta_phi_tune * phi_options(0)/(1 - phi_options(0)), beta_phi_tune);
  double a1 =
    exp(CalcLogPhiDensity(phi_options(1), p1, alpha_p, beta_p, w_cube) -
	CalcLogPhiDensity(phi_options(0), p1, alpha_p, beta_p, w_cube));
  double a2 =
    exp(R::dbeta(phi_options(0), beta_phi_tune * phi_options(1)/(1 - phi_options(1)), beta_phi_tune, 1) -
	R::dbeta(phi_options(1), beta_phi_tune * phi_options(0)/(1 - phi_options(0)), beta_phi_tune, 1));
  NumericVector probs(2);
  arma::vec phi_out(1);
  phi_out(0) = phi_options(1);
  if (a1*a2 < 0)
    throw std::domain_error("a1*a2 < 0");
  if (a1*a2 < 1) {
    probs(1) = a1*a2;
    probs(0) = 1 - a1*a2;
    phi_out = Rcpp::RcppArmadillo::sample(phi_options, 1, false, probs);
  }
  return as_scalar(phi_out);
}

double SampleSigsq(const arma::cube w_cube,
		   const arma::cube mu_cube,
		   const arma::cube mutil_cube,
		   const arma::cube K_cube,
		   const double beta_iter,
		   const double a_sigsq,
		   const double b_sigsq) {
  //
  // Args:
  //
  //
  int n_animals = w_cube.n_rows;
  int T = w_cube.n_slices;
  // changed distribution at time t = 1
  // arma::mat mu_1 =
  //   mu_cube(0, 0, 0, size(mu_cube.n_rows, 1, mu_cube.n_slices));
  // arma::rowvec mu_1_row = vectorise(mu_1, 1);
  // double b_addon_1 = as_scalar(mu_1_row * mu_1_row.t());
  double b_addon_t = 0;
  for (int t = 0; t < (T - 1); t++) {
    arma::mat mu_difference =
      mu_cube(0, t + 1, 0, size(mu_cube.n_rows, 1, mu_cube.n_slices)) -
      mu_cube(0, t, 0, size(mu_cube.n_rows, 1, mu_cube.n_slices)) -
      beta_iter *
      mutil_cube(0, t, 0, size(mutil_cube.n_rows, 1, mutil_cube.n_slices));
    arma::rowvec mu_difference_row = vectorise(mu_difference, 1);
    b_addon_t = b_addon_t + as_scalar(mu_difference_row *
				      K_cube.slice(t + 1) *
				      mu_difference_row.t());
  }
  double scale = pow(b_sigsq + 0.5 * b_addon_t, -1);
  double sigsq_inv = R::rgamma(a_sigsq + (n_animals * (T - 1)), scale);
  return pow(sigsq_inv, -1);
}

double SampleAlpha(const arma::cube w_cube,
		   const arma::cube mu_cube,
		   const arma::cube mutil_cube,
		   const arma::cube K_cube,
		   const double alpha_iter,
		   const double beta_iter,
		   const double sigsq_iter,
		   const double alpha_alpha,
		   const double beta_alpha,
		   const double sigsq_alpha_tune) {
  //
  // Args:
  //
  NumericVector alpha_options(2);
  alpha_options(0) = alpha_iter; /* previous value */
  alpha_options(1) = R::rnorm(alpha_iter,
			      sqrt(sigsq_alpha_tune)); /* proposed value */
  double a1 =
    exp(CalcLogAlphaDensity(alpha_options(1), alpha_alpha, beta_alpha,
			    w_cube, mu_cube, mutil_cube,
			    K_cube, beta_iter, sigsq_iter) -
	CalcLogAlphaDensity(alpha_options(0), alpha_alpha, beta_alpha,
			    w_cube, mu_cube, mutil_cube,
			    K_cube, beta_iter, sigsq_iter));
  NumericVector probs(2);
  arma::vec alpha_out(1);
  alpha_out(0) = alpha_options(1);
  if (a1 < 0)
    throw std::domain_error("a1 < 0");
  if (a1 < 1) {
    probs(1) = a1;
    probs(0) = 1 - a1;
    alpha_out = Rcpp::RcppArmadillo::sample(alpha_options, 1, false, probs);
  }
  return as_scalar(alpha_out);
}

// double SampleAlphaOLD(const arma::cube w_cube,
// 		   const arma::cube mu_cube,
// 		   const arma::cube mubar_cube,
// 		   const arma::cube mutil_cube,
// 		   const arma::cube K_cube,
// 		   const double beta_iter,
// 		   const double sigsq_iter,
// 		   const double sigsq_alpha,
// 		   const double mu_alpha) {
//   //
//   // WARNING: THIS IS DEPRICATED AND NOT IN USE ANYMORE
//   //
//   // Args:
//   //
//   //
//   int n_animals = w_cube.n_rows;
//   int T = w_cube.n_slices;
//   int d = mu_cube.n_slices;
//   arma::cube hbar_cube = MakeHbar(w_cube, mu_cube, mutil_cube, beta_iter);
//   arma::cube h_cube = mu_cube(0, 1, 0, size(n_animals, T - 1, d)) -
//     mu_cube(0, 0, 0, size(n_animals, T - 1, d)) -
//     beta_iter * mutil_cube(0, 0, 0, size(n_animals, T - 1, d));
//   double sigsq_alpha_sum = 0;
//   double mu_alpha_sum = 0;
//   for (int i = 0; i < n_animals; i++) {
//     arma::vec mubar_i1 = mubar_cube(i, 0, 0, size(1, 1, d));
//     double mubar_norm = as_scalar(mubar_i1.t() * mubar_i1);
//     double sigsq_new_term_1 = (K_cube(2*i, 2*i, 0) / sigsq_iter) * mubar_norm;
//     sigsq_alpha_sum = sigsq_alpha_sum + sigsq_new_term_1;
//     arma::vec mu_i1 = mu_cube(i, 0, 0, size(1, 1, d));
//     double mubar_mu = as_scalar(mubar_i1.t() * mu_i1);
//     double mu_new_term_1 = (K_cube(2*i, 2*i, 0) / sigsq_iter) * mubar_mu;
//     mu_alpha_sum = mu_alpha_sum + mu_new_term_1;
//     for (int t = 1; t < T; t++) {
//       arma::vec hbar_it = hbar_cube(i, t - 1, 0, size(1, 1, d));
//       double hbar_norm = as_scalar(hbar_it.t() * hbar_it);
//       double sigsq_new_term = (K_cube(2*i, 2*i, t) / sigsq_iter) * hbar_norm;
//       sigsq_alpha_sum = sigsq_alpha_sum + sigsq_new_term;
//       arma::vec h_it = h_cube(i, t - 1, 0, size(1, 1, d));
//       double hbar_h = as_scalar(hbar_it.t() * h_it);
//       double mu_new_term = (K_cube(2*i, 2*i, t) / sigsq_iter) * hbar_h;
//       mu_alpha_sum = mu_alpha_sum + mu_new_term;
//     }
//   }
//   double sigsq_alpha_star = pow((sigsq_alpha_sum + (pow(sigsq_alpha, -1))), -1);
//   double sigma_alpha_star = pow((sigsq_alpha_sum + (pow(sigsq_alpha, -1))), -0.5);
//   double mu_alpha_star = sigsq_alpha_star * (mu_alpha_sum + mu_alpha);
//   double alpha = SampleTruncatedNormal(mu_alpha_star, sigma_alpha_star, -1, 1);
//   return alpha;
//   /* return mu_alpha_star; */
// }

double SampleBeta(const arma::cube mu_cube,
		  const arma::cube mutil_cube,
		  const arma::cube K_cube,
		  const double sigsq_iter,
		  const double sigsq_beta,
		  const double mu_beta) {
  //
  // Args:
  //
  //
  int n_animals = mu_cube.n_rows;
  int T = mu_cube.n_cols;
  int d = mu_cube.n_slices;
  arma::cube Q_cube = (pow(sigsq_iter, -1)) * K_cube;
  double sigsq_beta_sum = 0;
  double mu_beta_sum = 0;
  for (int t = 1; t < T; t++) {
    arma::mat mutil_t_mat = mutil_cube(0, t - 1, 0, size(n_animals, 1, d));
    arma::rowvec mutil_t_rowvec = vectorise(mutil_t_mat, 1);
    arma::mat mu_difference_t_mat =
      mu_cube(0, t, 0, size(n_animals, 1, d)) -
      mu_cube(0, t - 1, 0, size(n_animals, 1, d));
    arma::rowvec mu_difference_t_rowvec = vectorise(mu_difference_t_mat, 1);
    double sigsq_new_term =
      as_scalar(mutil_t_rowvec * Q_cube.slice(t) * mutil_t_rowvec.t());
    sigsq_beta_sum = sigsq_beta_sum + sigsq_new_term;
    double mu_new_term =
      as_scalar(mutil_t_rowvec * Q_cube.slice(t) * mu_difference_t_rowvec.t());
    mu_beta_sum = mu_beta_sum + mu_new_term;
  }
  double sigsq_beta_star = pow((sigsq_beta_sum + (pow(sigsq_beta, -1))), -1);
  double sigma_beta_star = sqrt(sigsq_beta_star);
  double mu_beta_star = sigsq_beta_star * (mu_beta_sum + mu_beta);
  double out = R::rnorm(mu_beta_star, sigma_beta_star);
  return out;
}

int SampleWij1(const int i,
	       const int j,
	       const arma::cube w_cube,
	       const arma::cube mu_cube,
	       const arma::cube mubar_cube,
	       const arma::cube mutil_cube,
	       const double alpha_iter,
	       const double beta_iter,
	       const double sigsq_iter,
	       const double c_iter,
	       const arma::mat conditional_M_mat,
	       const double p1) {
  //
  // Args:
  //
  //
  double Probability_1_1 =
    (MakeNormDensityATtp1(i, j, 1, 1, w_cube, mu_cube, mutil_cube,
    			  alpha_iter, beta_iter, sigsq_iter, c_iter)) *
    (MakeNormDensityATtp1(j, i, 1, 1, w_cube, mu_cube, mutil_cube,
    			  alpha_iter, beta_iter, sigsq_iter, c_iter)) *
    (pow(conditional_M_mat(1, 1), w_cube(i - 1, j - 1, 1))) *
    (pow(conditional_M_mat(0, 1), 1 - w_cube(i - 1, j - 1, 1))) * p1;
  double Probability_1_0 =
    (MakeNormDensityATtp1(i, j, 1, 0, w_cube, mu_cube, mutil_cube,
    			  alpha_iter, beta_iter, sigsq_iter, c_iter)) *
    (MakeNormDensityATtp1(j, i, 1, 0, w_cube, mu_cube, mutil_cube,
    			  alpha_iter, beta_iter, sigsq_iter, c_iter)) *
    (pow(conditional_M_mat(1, 0), w_cube(i - 1, j - 1, 1))) *
    (pow(conditional_M_mat(0, 0), 1 - w_cube(i - 1, j - 1, 1))) * (1 - p1);
  if (Probability_1_1 == Probability_1_0) {
    Probability_1_1 = 0.5;
    Probability_1_0 = 0.5;
  }
  arma::vec w_out(1);
  NumericVector w_options(2);
  w_options(0) = 0;
  w_options(1) = 1;
  NumericVector probs(2);
  probs(0) = Probability_1_0;
  probs(1) = Probability_1_1;
  w_out = Rcpp::RcppArmadillo::sample(w_options, 1, false, probs);
  return as_scalar(w_out);
}

int SampleWijt(const int i,
	       const int j,
	       const int t,
	       const arma::cube w_cube,
	       const arma::cube mu_cube,
	       const arma::cube mutil_cube,
	       const double alpha_iter,
	       const double beta_iter,
	       const double sigsq_iter,
	       const double c_iter,
	       const arma::mat conditional_M_mat) {
  //
  // Args:
  //    iter: iteration in the MCMC chain min value is 1 NOT 0
  //    i: min value is 1 NOT 0
  //    j: min value is 1 NOT 0
  //    t: min value is 1 NOT 0
  //    w_cube: cube of dim (n_animals, n_animals, T)
  //    mu_cube: cube of dim (n_animals, 2 (lat, long), T)
  //    mutil_cube: cube of dim (n_animals, 2 (lat, long), T). This is vector
  //                pointing to the centroid location from the previous
  //                location.
  //    alpha_iter: vector of length n_iterations
  //    beta_iter:
  //    sigsq_iter:
  //    c_iter:
  //    conditional_M_mat: see description in MakeConditionalPMatrix()
  //
  // Return: either 0 (no edge) or 1 (edge)
  //
  double Probability_t_1 =
    (MakeNormDensityATt(i, j, t, 1, w_cube, mu_cube, mutil_cube,
			alpha_iter, beta_iter, sigsq_iter, c_iter)) *
    (MakeNormDensityATt(j, i, t, 1, w_cube, mu_cube, mutil_cube,
			alpha_iter, beta_iter, sigsq_iter, c_iter)) *
    (MakeNormDensityATtp1(i, j, t, 1, w_cube, mu_cube, mutil_cube,
			  alpha_iter, beta_iter, sigsq_iter, c_iter)) *
    (MakeNormDensityATtp1(j, i, t, 1, w_cube, mu_cube, mutil_cube,
			  alpha_iter, beta_iter, sigsq_iter, c_iter)) *
    (pow(conditional_M_mat(1, 0), 1 - w_cube(i - 1, j - 1, t - 2))) *
    (pow(conditional_M_mat(1, 1), w_cube(i - 1, j - 1, t) +
	 w_cube(i - 1, j - 1, t - 2))) *
    (pow(conditional_M_mat(0, 1), 1 - w_cube(i - 1, j - 1, t)));
  double Probability_t_0 =
    (MakeNormDensityATt(i, j, t, 0, w_cube, mu_cube, mutil_cube,
			alpha_iter, beta_iter, sigsq_iter, c_iter)) *
    (MakeNormDensityATt(j, i, t, 0, w_cube, mu_cube, mutil_cube,
			alpha_iter, beta_iter, sigsq_iter, c_iter)) *
    (MakeNormDensityATtp1(i, j, t, 0, w_cube, mu_cube, mutil_cube,
			  alpha_iter, beta_iter, sigsq_iter, c_iter)) *
    (MakeNormDensityATtp1(j, i, t, 0, w_cube, mu_cube, mutil_cube,
			  alpha_iter, beta_iter, sigsq_iter, c_iter)) *
    (pow(conditional_M_mat(1, 0), w_cube(i - 1, j - 1, t))) *
    (pow(conditional_M_mat(0, 0), 1 - w_cube(i - 1, j - 1, t) +
	 1 - w_cube(i - 1, j - 1, t - 2))) *
    (pow(conditional_M_mat(0, 1), w_cube(i - 1, j - 1, t - 2)))
    ;
  if (Probability_t_1 == Probability_t_0) {
    Probability_t_1 = 0.5;
    Probability_t_0 = 0.5;
  }
  arma::vec w_out(1);
  NumericVector w_options(2);
  w_options(0) = 0;
  w_options(1) = 1;
  NumericVector probs(2);
  probs(0) = Probability_t_0;
  probs(1) = Probability_t_1;
  w_out = Rcpp::RcppArmadillo::sample(w_options, 1, false, probs);
  return as_scalar(w_out);
}

int SampleWijT(const int i,
	       const int j,
	       const arma::cube w_cube,
	       const arma::cube mu_cube,
	       const arma::cube mutil_cube,
	       const double alpha_iter,
	       const double beta_iter,
	       const double sigsq_iter,
	       const double c_iter,
	       const arma::mat conditional_M_mat) {
  //
  // Args:
  //    iter: iteration in the MCMC chain min value is 1 NOT 0
  //    i: min value is 1 NOT 0
  //    j: min value is 1 NOT 0
  //    t: min value is 1 NOT 0
  //    w_cube: cube of dim (n_animals, n_animals, T)
  //    mu_cube: cube of dim (n_animals, 2 (lat, long), T)
  //    mutil_cube: cube of dim (n_animals, 2 (lat, long), T). This is vector
  //                pointing to the centroid location from the previous
  //                location.
  //    alpha: vector of length n_iterations
  //    beta: vector of length n_iterations
  //    sigsq: vector of length n_iterations
  //    phi: vector of length n_iterations
  //
  // Return: either 0 (no edge) or 1 (edge)
  //
  int T = mu_cube.n_cols;
  double Probability_T_1 =
    // (MakeNormDensityATt(i, j, T, 1, w_cube, mu_cube, mutil_cube,
    // 			alpha_iter, beta_iter, sigsq_iter, c_iter)) *
    // (MakeNormDensityATt(j, i, T, 1, w_cube, mu_cube, mutil_cube,
    // 			alpha_iter, beta_iter, sigsq_iter, c_iter)) *
    (pow(conditional_M_mat(1, 1), w_cube(i - 1, j - 1, T - 2))) *
    (pow(conditional_M_mat(1, 0), 1 - w_cube(i - 1, j - 1, T - 2)));
  double Probability_T_0 =
    // (MakeNormDensityATt(i, j, T, 0, w_cube, mu_cube, mutil_cube,
    // 			alpha_iter, beta_iter, sigsq_iter, c_iter)) *
    // (MakeNormDensityATt(j, i, T, 0, w_cube, mu_cube, mutil_cube,
    // 			alpha_iter, beta_iter, sigsq_iter, c_iter)) *
    (pow(conditional_M_mat(0, 0), 1 - w_cube(i - 1, j - 1, T - 2))) *
    (pow(conditional_M_mat(0, 1), w_cube(i - 1, j - 1, T - 2)));
  // if (Probability_T_1 == Probability_T_0) {
  //   Probability_T_1 = 0.5;
  //   Probability_T_0 = 0.5;
  // }
  vec w_out(1);
  NumericVector w_options(2);
  w_options(0) = 0;
  w_options(1) = 1;
  NumericVector probs(2);
  probs(0) = Probability_T_0;
  probs(1) = Probability_T_1;
  w_out = Rcpp::RcppArmadillo::sample(w_options, 1, false, probs);
  return as_scalar(w_out);
}

double SampleC(const double a_c,
	       const double b_c,
	       const arma::cube w_cube,
	       const arma::cube mu_cube,
	       const arma::cube mutil_cube,
	       const double beta_iter,
	       const double sigsq_iter) {
  int T = mu_cube.n_cols;
  int n_animals = mu_cube.n_rows;
  int d = mu_cube.n_slices;
  int a_sum = 0;
  double b_sum = 0;
  for (int t = 1; t < T; t++) {
    for (int i = 0; i < n_animals; i++) {
      arma::vec w_i_t = reshape(w_cube(i, 0, t, size(1, n_animals, 1)),
				n_animals, 1, 1);
      if (sum(w_i_t) == 0) {
	a_sum = a_sum + 1;
	arma::cube h_i_t_cube =
	  mu_cube(i, t, 0, size(1, 1, d)) -
	  mu_cube(i, t - 1, 0, size(1, 1, d)) -
	  beta_iter * mutil_cube(i, t - 1, 0, size(1, 1, d));
	arma::rowvec h_i_t = reshape(h_i_t_cube, 1, d, 1);
	b_sum = b_sum + as_scalar(h_i_t * h_i_t.t());
      }
    }
  }
  double scale = pow(b_c + 0.5 * b_sum / sigsq_iter, -1);
  double c = R::rgamma(a_c + a_sum, scale);
  return c;
}

double SampleP1(const double p1_iter,
		const double phi_iter,
		const arma::cube w_cube,
		const double alpha_p1,
		const double beta_p1,
		const double beta_p1_tune) {
  NumericVector p1_options(2);
    /* previous value */
  p1_options(0) = p1_iter;
  /* proposed value */
  p1_options(1) = R::rbeta(beta_p1_tune * p1_options(0)/(1 - p1_options(0)),
			   beta_p1_tune);
  double a1 =
    exp(CalcLogP1Density(p1_options(1), phi_iter, w_cube, alpha_p1, beta_p1) -
	CalcLogP1Density(p1_options(0), phi_iter, w_cube, alpha_p1, beta_p1));
  double a2 =
    exp(R::dbeta(p1_options(0),
		 beta_p1_tune * p1_options(1)/(1 - p1_options(1)),
		 beta_p1_tune, 1) -
	R::dbeta(p1_options(1),
		 beta_p1_tune * p1_options(0)/(1 - p1_options(0)),
		 beta_p1_tune, 1));
  NumericVector probs(2);
  arma::vec p1_out(1);
  p1_out(0) = p1_options(1);
  if (a1*a2 < 0)
    throw std::domain_error("a1*a2 < 0");
  if (a1*a2 < 1) {
    probs(1) = a1*a2;
    probs(0) = 1 - a1*a2;
    p1_out = Rcpp::RcppArmadillo::sample(p1_options, 1, false, probs);
  }
  return as_scalar(p1_out);
}

///////////////////
// END FUNCTIONS //
///////////////////

////////////////
// BEGIN LOOP //
////////////////

//' fit the model using mcmc
//'
//' @param N_iterations if fitting all parameters, this takes approximately 1
//' minutes per 100 iterations on my laptop.
//' @param mu this is an array of dimension \code{c(n.indiv, T, 2, K)} where
//' \code{K} is some number of draws from the distribution of \eqn{s|\mu}.
//' @param phi vector of length \code{N_iterations}
//' @param alpha vector of length \code{N_iterations}
//' @param beta vector of length \code{N_iterations}
//' @param p1 vector of length \code{N_iterations}
//' @param sigsq vector of length \code{N_iterations}
//' @param c vector of length \code{N_iterations}
//' @param w array of dimension \code{(n.indiv, n.indiv, T, N_iterations)}
//' @param alpha_p1 hyperparameters
//' @param beta_p1 hyperparameters
//' @param alpha_phi hyperparameters
//' @param beta_phi hyperparameters
//' @param alpha_alpha hyperparameters
//' @param beta_alpha hyperparameters
//' @param mu_beta hyperparameters
//' @param sigsq_beta hyperparameters
//' @param a_sigsq hyperparameters
//' @param b_sigsq hyperparameters
//' @param a_c hyperparameters
//' @param b_c hyperparameters
//' @param beta_p1_tune tuning parameter for adjusted beta proposal
//' @param beta_phi_tune tuning parameter for beta proposal
//' @param sigsq_alpha_tune tuning parameter for normal proposal
//' @return list of chains for each parameter
//' @useDynLib mingle
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
List fit_mingle(const int N_iterations,
		     const NumericVector mu, //mean process|data
		     arma::vec phi, //storage
		     arma::vec alpha,
		     arma::vec beta,
		     arma::vec p1,
		     arma::vec sigsq,
		     arma::vec c,
		     NumericVector w,
		     const double alpha_p1, //prior parameters
		     const double beta_p1,
		     const double alpha_phi,
		     const double beta_phi,
		     const double alpha_alpha,
		     const double beta_alpha,
		     const double mu_beta,
		     const double sigsq_beta,
		     const double a_sigsq,
		     const double b_sigsq,
		     const double a_c,
		     const double b_c,
		     const double beta_p1_tune, //tuning parameters
		     const double beta_phi_tune,
		     const double sigsq_alpha_tune){
  arma::field<cube> mu_field = ReadField(mu);
  int n_animals = mu_field(1).n_rows;
  int T = mu_field(1).n_cols;
  arma::field<cube> w_field = ReadField(w);
  arma::mat conditional_M_mat = MakeConditionalMMatrix(p1(1), phi(1));
  arma::cube K_cube = MakePrecisionKernel(alpha(1), w_field(1), c(1));
  free(w);
  int mu_index = floor(R::runif(0, mu_field.n_elem));
  arma::cube mu_cube = mu_field(mu_index);
  arma::cube w_cube = w_field(1);
  arma::cube mubar_cube = MakeMubar(w_cube, mu_cube);
  arma::cube mutil_cube = MakeMuTilde(w_cube, mu_cube, mubar_cube);
  for (int iter = 2; iter < N_iterations; iter++){
    mu_index = floor(R::runif(0, mu_field.n_elem));
    mu_cube = mu_field(mu_index);
    w_cube = w_field(iter - 1);
    mubar_cube = MakeMubar(w_cube, mu_cube);
    mutil_cube = MakeMuTilde(w_cube, mu_cube, mubar_cube);
    //
    // c
    //
    if (R_IsNA(c(0))) {
      c(iter) = c(1);
    } else {
      c(iter) = SampleC(a_c, b_c, w_cube, mu_cube, mutil_cube,
        beta(iter - 1), sigsq(iter - 1));
    }
    K_cube = MakePrecisionKernel(alpha(iter - 1), w_cube, c(iter));
    //
    // p1
    //
    if (R_IsNA(p1(0))) {
      p1(iter) = p1(1);
    } else {
      p1(iter) = SampleP1(p1(iter - 1), phi(iter - 1),
         w_cube, alpha_p1, beta_p1, beta_p1_tune);
    }
    //
    // phi
    //
    if (R_IsNA(phi(0))) {
      phi(iter) = phi(1);
    } else {
      phi(iter) = SamplePhi(w_cube, beta_phi_tune, p1(iter),
          alpha_phi, beta_phi, phi(iter - 1));
    }
    conditional_M_mat = MakeConditionalMMatrix(p1(iter), phi(iter));
    //
    // sigsq
    //
    if (R_IsNA(sigsq(0))) {
      sigsq(iter) = sigsq(1);
    } else {
      sigsq(iter) = SampleSigsq(w_cube, mu_cube, mutil_cube, K_cube,
            beta(iter - 1), a_sigsq, b_sigsq);
    }
    //
    // beta
    //
    if (R_IsNA(beta(0))) {
      beta(iter) = beta(1);
    } else {
      beta(iter) = SampleBeta(mu_cube, mutil_cube, K_cube,
           sigsq(iter), sigsq_beta, mu_beta);
    }
    //
    // alpha
    //
    if (R_IsNA(alpha(0))) {
      alpha(iter) = alpha(1);
    } else {
      alpha(iter) = SampleAlpha(w_cube, mu_cube, mutil_cube, K_cube,
            alpha(iter - 1), beta(iter), sigsq(iter),
            alpha_alpha, beta_alpha, sigsq_alpha_tune);
    }
    K_cube = MakePrecisionKernel(alpha(iter), w_cube, c(iter));
    //
    // w
    //
    if (R_IsNA(w_field(0)(0, 0, 0))) {
      w_cube = w_field(1);
    } else {
      for (int i = 0; i < n_animals - 1; i++) {
        for (int j = i + 1; j < n_animals; j++) {
          //
          // wij1
          //
          int wij1 = SampleWij1(i + 1, j + 1, w_cube, mu_cube, mubar_cube, mutil_cube,
                                alpha(iter), beta(iter), sigsq(iter), c(iter),
                                conditional_M_mat, p1(iter));
          w_cube(j, i, 0) = wij1;
          w_cube(i, j, 0) = wij1;
          for (int t = 1; t < (T - 1); t++) {
            //
            // wijt
            //
            int wijt = SampleWijt(i + 1, j + 1, t + 1, w_cube, mu_cube, mutil_cube,
                                  alpha(iter), beta(iter), sigsq(iter), c(iter),
                                  conditional_M_mat);
            w_cube(i, j, t) = wijt;
            w_cube(j, i, t) = wijt;
          }
          //
          // wijT
          //
          int wijT = SampleWijT(i + 1, j + 1, w_cube, mu_cube, mutil_cube,
                                alpha(iter), beta(iter), sigsq(iter), c(iter),
                                conditional_M_mat);
          w_cube(i, j, T - 1) = wijT;
          w_cube(j, i, T - 1) = wijT;
        }
      }
    }
    w_field(iter) = w_cube;
  }
  return List::create(Named("c") = c,
                      Named("p1") = p1,
                      Named("phi") = phi,
                      Named("sigsq") = sigsq,
                      Named("beta") = beta,
                      Named("alpha") = alpha,
                      Named("W") = w_field);
}

//////////////
// END LOOP //
//////////////
