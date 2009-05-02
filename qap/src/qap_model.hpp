#pragma once

#include <cassert>
#include <valarray>
#include <algorithm>
#include <tr1/random>
#include <metslib/mets.h>

class qap_move;

class qap_model : public mets::feasible_solution
{
protected:
  std::valarray<int> index;
  std::valarray< std::valarray<int> > a;
  std::valarray< std::valarray<int> > b;
  uint64_t c_m;
  
public:
  
  qap_model() : index(), a(), b(), c_m(0) {};
  
  friend class qap_move;
  
  mets::gol_type cost_function() const 
  {
    return 2.*c_m;
  }
  
  virtual qap_model&
  operator=(const qap_model& o)
  {
    unsigned int n = o.index.size();
    index.resize(n);
    a.resize(n);
    b.resize(n);
    for(unsigned int ii = 0; ii != n; ++ii)
      { a[ii].resize(n); b[ii].resize(n); }
    index = o.index;
    a = o.a;
    b = o.b;
    c_m = o.c_m;
    return *this;
  }
  
  // Copy operator
  feasible_solution&
  operator=(const mets::feasible_solution& sol)
  {
    const qap_model* o = dynamic_cast<const qap_model*>(&sol);
    unsigned int n = o->index.size();
    index.resize(n);
    a.resize(n);
    b.resize(n);
    for(unsigned int ii = 0; ii != n; ++ii)
      { a[ii].resize(n); b[ii].resize(n); }
    index = o->index;
    a = o->a;
    b = o->b;
    c_m = o->c_m;
    return *this;
  }
  
  unsigned int size() const { return index.size(); }
  
  /// @brief Generate a random starting point.
  template<typename rng_t>
  void
  random_shuffle(rng_t rng)
  {
    std::tr1::uniform_int<> unigen(0, index.size());
    std::tr1::variate_generator<rng_t, 
      std::tr1::uniform_int<> >gen(rng, unigen);
    std::random_shuffle(&index[0], &index[index.size()], gen);
    c_m = half_cost_calculator();
  }
  
  /// @brief Generate a random starting point.
  template<typename rng_t>
  void
  perturbate(unsigned int n, rng_t rng)
  {
    std::tr1::uniform_int<> int_range;
    for(int ii=0; ii!=n;++ii) {
      int p1 = int_range(rng, a.size());
      int p2 = int_range(rng, a.size());
      while(p1 == p2) 
	p2 = int_range(rng, a.size());
      swap(p1, p2);
    }
  }
  
  void
  swap(int i, int j)
  {
    for(unsigned int ii=0; ii != a.size(); ++ii)
      {
	c_m -= a[i][ii] * b[index[i]][index[ii]];
	c_m -= a[j][ii] * b[index[j]][index[ii]];
      }
    std::swap(index[i],index[j]);
    for(unsigned int ii=0; ii != a.size(); ++ii)
      {
	c_m += a[i][ii] * b[index[i]][index[ii]];
	c_m += a[j][ii] * b[index[j]][index[ii]];
     }
    // assert(c_m == half_cost_calculator());
  }

  friend std::ostream& operator<<(std::ostream& os, const qap_model& p);
  friend std::istream& operator>>(std::istream& is, qap_model& p);
  
protected:
  
  // Straight cost calculator
  uint64_t half_cost_calculator() const
  {
    uint64_t sum = 0;
    for(unsigned int jj = 0; jj != a.size(); ++jj)
      for(unsigned int ii = 0; ii != jj; ++ii)
	sum += a[ii][jj] * b[index[ii]][index[jj]];
    return sum;
  }

};

// Input/Output functions
std::ostream&
operator<<(std::ostream& os, const qap_model& qap)
{
  for(unsigned int ii = 0; ii != qap.index.size(); ++ii) 
    os << (qap.index[ii]+1) << " ";
  return os;
}

std::istream&
operator>>(std::istream& is, qap_model& qap)
{
  unsigned int n;
  is >> n;
  qap.index.resize(n);
  qap.a.resize(n);
  qap.b.resize(n);
  for(unsigned int ii = 0; ii != n; ++ii)
    { qap.index[ii] = ii; qap.a[ii].resize(n); qap.b[ii].resize(n); }

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
  qap.c_m = qap.half_cost_calculator();
  return is;
}
