#ifndef LAMBDA_LANCZOS_H_
#define LAMBDA_LANCZOS_H_

#include <vector>
#include <functional>
#include <cassert>
#include <limits>
#include <cmath>
#include <numeric>
#include <random>
#include <lambda_lanczos/lambda_lanczos_util.hpp>


namespace lambda_lanczos {

namespace {
template<typename T>
using vector = std::vector<T>;

template<typename T>
using complex = std::complex<T>;

using std::cout;
using std::endl;
using std::fill;
using std::abs;
using std::begin;
using std::end;
using namespace lambda_lanczos_util;
}



/*
 * "Partially specialization of function" is not allowed, 
 * so here it is mimicked by wrapping the "init" function with a class template.
 */
template <typename T>
struct VectorRandomInitializer {
public:
  static void init(vector<T>&);
};


template <typename T>
struct VectorRandomInitializer<complex<T>> {
public:
  static void init(vector<complex<T>>&);
};


template <typename T>
class LambdaLanczos {
public:
  LambdaLanczos(std::function<void(const vector<T>&, vector<T>&)> mv_mul, int matrix_size, bool find_maximum);
  LambdaLanczos(std::function<void(const vector<T>&, vector<T>&)> mv_mul, int matrix_size) : LambdaLanczos(mv_mul, matrix_size, false) {}
  
  int matrix_size;
  int max_iteration;
  real_t<T> eps = minimum_effective_decimal<real_t<T>>() * 1e3;
  real_t<T> tridiag_eps_ratio = 1e-1;
  int initial_vector_size = 200;
  bool find_maximum = false;
  real_t<T> eigenvalue_offset = 0.0;

  std::function<void(const vector<T>&, vector<T>&)> mv_mul;
  std::function<void(vector<T>&)> init_vector = VectorRandomInitializer<T>::init;

  int run(real_t<T>&, vector<T>&) const;

private:
  static void schmidt_orth(vector<T>&, const vector<vector<T>>&);
  real_t<T> find_minimum_eigenvalue(const vector<real_t<T>>&,
				    const vector<real_t<T>>&) const;
  real_t<T> find_maximum_eigenvalue(const vector<real_t<T>>&,
				    const vector<real_t<T>>&) const;
  static real_t<T> tridiagonal_eigen_limit(const vector<real_t<T>>&,
					   const vector<real_t<T>>&);
  static int num_of_eigs_smaller_than(real_t<T>,
				      const vector<real_t<T>>&,
				      const vector<real_t<T>>&);
};




/* Implementation */

template <typename T>
inline LambdaLanczos<T>::LambdaLanczos(std::function<void(const vector<T>&, vector<T>&)> mv_mul,
				int matrix_size, bool find_maximum) {
  this->mv_mul = mv_mul;
  this->matrix_size = matrix_size;
  this->max_iteration = matrix_size;
  this->find_maximum = find_maximum;
}


template <typename T>
inline int LambdaLanczos<T>::run(real_t<T>& eigvalue, vector<T>& eigvec) const {
  assert(0 < this->tridiag_eps_ratio && this->tridiag_eps_ratio < 1);
  
  vector<vector<T>> u;     // Lanczos vectors
  vector<real_t<T>> alpha; // Diagonal elements of an approximated tridiagonal matrix
  vector<real_t<T>> beta;  // Subdiagonal elements of an approximated tridiagonal matrix

  const int n = this->matrix_size;

  u.reserve(this->initial_vector_size);
  alpha.reserve(this->initial_vector_size);
  beta.reserve(this->initial_vector_size);

  u.emplace_back(n, 0.0); // Same as u.push_back(vector<T>(n, 0.0))
  
  vector<T> vk(n, 0.0);
  
  real_t<T> alphak = 0.0;
  alpha.push_back(alphak);
  real_t<T> betak = 0.0;
  beta.push_back(betak);
  
  vector<T> uk(n);
  this->init_vector(uk);
  normalize(uk);
  u.push_back(uk);

  real_t<T> ev, pev; // Calculated eigenvalue and previous one
  pev = std::numeric_limits<real_t<T>>::max();

  int itern = this->max_iteration;
  for(int k = 1;k <= this->max_iteration;k++) {
    /* vk = (A + offset*E)uk, here E is the identity matrix */
    for(int i = 0;i < n;i++) {
      vk[i] = uk[i]*this->eigenvalue_offset;
    }
    this->mv_mul(uk, vk);
    
    alphak = std::real(inner_prod(u.back(), vk));
    
    /* The inner product <uk|vk> is real.
     * Proof:
     *     <uk|vk> = <uk|A|uk>
     *   On the other hand its complex conjugate is
     *     <uk|vk>^* = <vk|uk> = <uk|A^*|uk> = <uk|A|uk>
     *   here the condition that matrix A is a symmetric (Hermitian) is used.
     *   Therefore
     *     <uk|vk> = <vk|uk>^*
     *   <uk|vk> is real.
     */
    
    alpha.push_back(alphak);
    
    for(int i = 0;i < n; i++) {
      uk[i] = vk[i] - betak*u[k-1][i] - alphak*u[k][i];
    }
    
    schmidt_orth(uk, u);
    
    betak = norm(uk);
    beta.push_back(betak);

    if(this->find_maximum) {
      ev = find_maximum_eigenvalue(alpha, beta);
    } else {
      ev = find_minimum_eigenvalue(alpha, beta);
    }

    const real_t<T> zero_threshold = minimum_effective_decimal<real_t<T>>()*1e-1;
    if(betak < zero_threshold) {
      u.push_back(uk);
      /* This element will never be accessed,
       * but this "push" guarantees u to always have one more element than 
       * alpha and beta do.
       */
      itern = k;
      break;
    }

    normalize(uk);
    u.push_back(uk);

    if(abs(ev-pev) < std::min(abs(ev), abs(pev))*this->eps) {
      itern = k;
      break;
    } else {
      pev = ev;
    }
  }

  eigvalue = ev - this->eigenvalue_offset;

  int m = alpha.size();
  vector<T> cv(m+1);
  cv[0] = 0.0;
  cv[m] = 0.0;
  cv[m-1] = 1.0;
  
  beta[m-1] = 0.0;

  if(eigvec.size() < n) {
    eigvec.resize(n);
  }
  
  for(int i = 0;i < n;i++) {
    eigvec[i] = cv[m-1]*u[m-1][i];
  }

  for(int k = m-2;k >= 1;k--) {
    cv[k] = ((ev - alpha[k+1])*cv[k+1] - beta[k+1]*cv[k+2])/beta[k];

    for(int i = 0;i < n;i++) {
      eigvec[i] += cv[k]*u[k][i];
    }
  }

  normalize(eigvec);

  return itern;
}


template <typename T>
inline void LambdaLanczos<T>::schmidt_orth(vector<T>& uorth, const vector<vector<T>>& u) {
  /* Vectors in u must be normalized, but uorth doesn't have to be. */
  
  int n = uorth.size();
  
  for(int k = 0;k < u.size();k++) {
    T innprod = inner_prod(uorth, u[k]);
    
    for(int i = 0;i < n;i++) {
      uorth[i] -= innprod * u[k][i];
    }
  }
}


template <typename T>
inline real_t<T> LambdaLanczos<T>::find_minimum_eigenvalue(const vector<real_t<T>>& alpha,
							   const vector<real_t<T>>& beta) const {
  real_t<T> eps = this->eps * this->tridiag_eps_ratio;
  real_t<T> pmid = std::numeric_limits<real_t<T>>::max();
  real_t<T> r = tridiagonal_eigen_limit(alpha, beta);
  real_t<T> lower = -r;
  real_t<T> upper = r;
  real_t<T> mid;
  int nmid; // Number of eigenvalues smaller than the "mid"
  
  while(upper-lower > std::min(abs(lower), abs(upper))*eps) {
    mid = (lower+upper)/2.0;
    nmid = num_of_eigs_smaller_than(mid, alpha, beta);
    if(nmid >= 1) {
      upper = mid;
    } else {
      lower = mid;
    }
    
    if(mid == pmid) {
      /* This avoids an infinite loop due to zero matrix */
      break;
    }
    pmid = mid;
  }

  return lower; // The "lower" almost equals the "upper" here.
}


template <typename T>
inline real_t<T> LambdaLanczos<T>::find_maximum_eigenvalue(const vector<real_t<T>>& alpha,
							   const vector<real_t<T>>& beta) const {
  real_t<T> eps = this->eps * this->tridiag_eps_ratio;
  real_t<T> pmid = std::numeric_limits<real_t<T>>::max();
  real_t<T> r = tridiagonal_eigen_limit(alpha, beta);
  real_t<T> lower = -r;
  real_t<T> upper = r;
  real_t<T> mid;
  int nmid; // Number of eigenvalues smaller than the "mid"

  int m = alpha.size() - 1;  /* Number of eigenvalues of the approximated triangular matrix,
				which equals the rank of it */
  
  
  while(upper-lower > std::min(abs(lower), abs(upper))*eps) {
    mid = (lower+upper)/2.0;
    nmid = num_of_eigs_smaller_than(mid, alpha, beta);
    
    if(nmid < m) {
      lower = mid;
    } else {
      upper = mid;
    }
    
    if(mid == pmid) {
      /* This avoids an infinite loop due to zero matrix */
      break;
    }
    pmid = mid;
  }

  return lower; // The "lower" almost equals the "upper" here.
}  


/*
 * Compute the upper bound of the absolute value of eigenvalues
 * by Gerschgorin theorem. This routine gives a rough upper bound,
 * but it is sufficient because the bisection routine using
 * the upper bound converges exponentially.
 */
template <typename T>
inline real_t<T> LambdaLanczos<T>::tridiagonal_eigen_limit(const vector<real_t<T>>& alpha,
							   const vector<real_t<T>>& beta) {
  real_t<T> r = l1_norm(alpha);
  r += 2*l1_norm(beta);
  
  return r;
}


/*
 * Algorithm from
 * Peter Arbenz et al. / "High Performance Algorithms for Structured Matrix Problems" /
 * Nova Science Publishers, Inc.
 */
template <typename T>
inline int LambdaLanczos<T>::num_of_eigs_smaller_than(real_t<T> c,
						      const vector<real_t<T>>& alpha,
						      const vector<real_t<T>>& beta) {
  real_t<T> q_i = 1.0;
  int count = 0;
  int m = alpha.size();
  
  for(int i = 1;i < m;i++){
    q_i = alpha[i] - c - beta[i-1]*beta[i-1]/q_i;
    if(q_i < 0){
      count++;
    }
    if(q_i == 0){
      q_i = minimum_effective_decimal<real_t<T>>();
    }
  }

  return count;
}


template <typename T>
inline void VectorRandomInitializer<T>::init(vector<T>& v) {
  std::random_device dev;
  std::mt19937 mt(dev());
  std::uniform_real_distribution<T> rand((T)(-1.0), (T)(1.0));

  int n = v.size();
  for(int i = 0;i < n;i++) {
    v[i] = rand(mt);
  }

  normalize(v);
}


template <typename T>
inline void VectorRandomInitializer<complex<T>>::init(vector<complex<T>>& v) {
  std::random_device dev;
  std::mt19937 mt(dev());
  std::uniform_real_distribution<T> rand((T)(-1.0), (T)(1.0));

  int n = v.size();
  for(int i = 0;i < n;i++) {
    v[i] = complex<T>(rand(mt), rand(mt));
  }

  normalize(v);
}

} /* namespace lambda_lanczos */

#endif /* LAMBDA_LANCZOS_H_ */
