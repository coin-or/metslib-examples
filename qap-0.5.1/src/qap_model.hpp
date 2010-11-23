#pragma once

#include <cassert>
#include <vector>
#include <algorithm>
#include <metslib/mets.hh>

class qap_model : public mets::permutation_problem
{
protected:
  std::vector< std::vector<int> > f_m;
  std::vector< std::vector<int> > d_m;
  
public:
  qap_model() : permutation_problem(0), f_m(), d_m() {};
  
  void copy_from(const mets::copyable& sol)
  {
    const qap_model& o = dynamic_cast<const qap_model&>(sol);
    permutation_problem::copy_from(sol);
    f_m = o.f_m;
    d_m = o.d_m;
  }

  /// @brief: swap move that does delta updates of the objective.
  ///
  /// This is much faster as it runs in O(n) instead of O(n^2) (it
  /// takes 8n multiplications instead of n^2 to compute the objective
  /// function after a swap).
  mets::gol_type
  evaluate_swap(int i, int j) const
  {
    assert(i!=j);
    double delta = 0.0;
    for(unsigned int ii=0; ii != f_m.size(); ++ii)
      {
	delta -= f_m[i][ii] * d_m[pi_m[i]][pi_m[ii]];
	delta -= f_m[ii][i] * d_m[pi_m[ii]][pi_m[i]];
	delta -= f_m[j][ii] * d_m[pi_m[j]][pi_m[ii]];
	delta -= f_m[ii][j] * d_m[pi_m[ii]][pi_m[j]];
	int ni = ii;
	if(ii==i) ni = j; else if(ii==j) ni = i;
	delta += f_m[i][ii] * d_m[pi_m[j]][pi_m[ni]];
	delta += f_m[ii][i] * d_m[pi_m[ni]][pi_m[j]];
	delta += f_m[j][ii] * d_m[pi_m[i]][pi_m[ni]];
	delta += f_m[ii][j] * d_m[pi_m[ni]][pi_m[i]];
      }
    return delta;
  }

  friend std::ostream& operator<<(std::ostream& os, const qap_model& p);
  friend std::istream& operator>>(std::istream& is, qap_model& p);
  
protected:
  
  // Full cost calculation
  mets::gol_type compute_cost() const
  {
    double sum = 0.0;
    for(unsigned int ii = 0; ii != pi_m.size(); ++ii)
      for(unsigned int jj = 0; jj != pi_m.size(); ++jj)
	sum += f_m[ii][jj] * d_m[pi_m[ii]][pi_m[jj]];
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
  qap.f_m.resize(n);
  qap.d_m.resize(n);
  for(unsigned int ii = 0; ii != n; ++ii)
    { qap.pi_m[ii] = ii; qap.f_m[ii].resize(n); qap.d_m[ii].resize(n); }

  for(unsigned int ii = 0; ii != n; ++ii)
    for(unsigned int jj = 0; jj != n; ++jj)
      {
	is >> qap.f_m[ii][jj];
	if(ii==jj) assert(qap.f_m[ii][jj] == 0);
      }
  
  for(unsigned int ii = 0; ii != n; ++ii)
    for(unsigned int jj = 0; jj != n; ++jj)
      {
	is >> qap.d_m[ii][jj];
	if(ii==jj) assert(qap.d_m[ii][jj] == 0);
      }
  qap.update_cost();
  return is;
}
