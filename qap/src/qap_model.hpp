#pragma once

#include <cassert>
#include <valarray>
#include <algorithm>
#if defined (WIN32)
#  include <random>
#else
#  include <tr1/random>
#endif
#include <metslib/mets.h>

class qap_model : public mets::permutation_problem
{
protected:
  std::valarray< std::valarray<int> > a;
  std::valarray< std::valarray<int> > b;
  int64_t c_m;
  
public:
  qap_model() : permutation_problem(0), a(), b(), c_m(0) {};
  
  /// @brief Returns the objective function value. This value is
  /// updated every time the variable is modified.
  mets::gol_type cost_function() const 
  {
    return (mets::gol_type)c_m;
  }
  
  void copy_from(const mets::feasible_solution& sol)
  {
    const qap_model* o = dynamic_cast<const qap_model*>(&sol);
    if(o)
      {
	mets::permutation_problem::copy_from(sol);
	unsigned int n = o->pi_m.size();
	a.resize(n);
	b.resize(n);
	for(unsigned int ii = 0; ii != n; ++ii)
	  { a[ii].resize(n); b[ii].resize(n); }
	a = o->a;
	b = o->b;
	c_m = o->c_m;
      }
    else
      {
	std::cerr << "Should not happen." << std::endl;
      }
  }

  void random_shuffle(std::tr1::mt19937& rng)
  {
    mets::random_shuffle(*this, rng);
    c_m = cost_calculator();
  }

  void perturbate(int n, std::tr1::mt19937& rng)
  {
    mets::perturbate(*this, n, rng);
  }

  /// @brief: swap move that does delta updates of the objective.
  ///
  /// This is much faster as it runs in O(n) instead of O(n^2) (it
  /// takes 8n multiplications instead of n^2 to compute the objective
  /// function after a swap).
  void
  swap(int i, int j)
  {
    for(unsigned int ii=0; ii != a.size(); ++ii)
      {
	c_m -= a[i][ii] * b[pi_m[i]][pi_m[ii]];
	c_m -= a[ii][i] * b[pi_m[ii]][pi_m[i]];
	c_m -= a[j][ii] * b[pi_m[j]][pi_m[ii]];
	c_m -= a[ii][j] * b[pi_m[ii]][pi_m[j]];
      }
    mets::permutation_problem::swap(i, j);
    for(unsigned int ii=0; ii != a.size(); ++ii)
      {
	c_m += a[i][ii] * b[pi_m[i]][pi_m[ii]];
	c_m += a[ii][i] * b[pi_m[ii]][pi_m[i]];
	c_m += a[j][ii] * b[pi_m[j]][pi_m[ii]];
	c_m += a[ii][j] * b[pi_m[ii]][pi_m[j]];
     }
    
    // c_m = cost_calculator();    
    assert(c_m == cost_calculator());
  }

  friend std::ostream& operator<<(std::ostream& os, const qap_model& p);
  friend std::istream& operator>>(std::istream& is, qap_model& p);
  
protected:
  
  // Straight cost calculator
  int64_t cost_calculator() const
  {
    int64_t sum = 0;
    for(unsigned int ii = 0; ii != pi_m.size(); ++ii)
      for(unsigned int jj = 0; jj != pi_m.size(); ++jj)
	sum += (a[ii][jj]) * b[pi_m[ii]][pi_m[jj]];
    return sum;
  }

};

// Input/Output functions
std::ostream&
operator<<(std::ostream& os, const qap_model& qap)
{
  for(unsigned int ii = 0; ii != qap.pi_m.size(); ++ii) 
    os << (qap.pi_m[ii]+1) << " ";
  return os;
}

std::istream&
operator>>(std::istream& is, qap_model& qap)
{
  unsigned int n;
  is >> n;
  qap.pi_m.resize(n);
  qap.a.resize(n);
  qap.b.resize(n);
  for(unsigned int ii = 0; ii != n; ++ii)
    { qap.pi_m[ii] = ii; qap.a[ii].resize(n); qap.b[ii].resize(n); }

  for(unsigned int ii = 0; ii != n; ++ii)
    for(unsigned int jj = 0; jj != n; ++jj)
      {
	is >> qap.a[ii][jj];
      }

  for(unsigned int ii = 0; ii != n; ++ii)
    for(unsigned int jj = 0; jj != n; ++jj)
      {
	is >> qap.b[ii][jj];
      }
  qap.c_m = qap.cost_calculator();
  return is;
}
