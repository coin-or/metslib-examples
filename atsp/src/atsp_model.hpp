#pragma once

#include <string>
#include <cassert>
#include <valarray>
#include <algorithm>
#if defined (WIN32)
#  include <random>
#else
#  include <tr1/random>
#endif
#include <metslib/mets.h>

class atsp_model : public mets::permutation_problem
{
protected:
  std::valarray< std::valarray<int> > matrix;
  // int64_t c_m;
  
public:
  atsp_model() : permutation_problem(0), matrix() /* c_m(0) */ {};
  
  /// @brief Returns the objective function value. This value is
  /// updated every time the variable is modified.
  mets::gol_type cost_function() const 
  {
    return (mets::gol_type)cost_calculator();
  }
  
  void copy_from(const mets::feasible_solution& sol)
  {
    const atsp_model* o = dynamic_cast<const atsp_model*>(&sol);
    if(o)
      {
	mets::permutation_problem::copy_from(sol);
	unsigned int n = o->pi_m.size()+1;
	matrix.resize(n);
	for(unsigned int ii = 0; ii != n; ++ii)
	  { matrix[ii].resize(n); }
	matrix = o->matrix;
	// c_m = o->c_m;
      }
    else
      {
	std::cerr << "Should not happen." << std::endl;
	::exit(-1);
      }
  }

  void random_shuffle(std::tr1::mt19937& rng)
  {
    mets::random_shuffle(*this, rng);
    // c_m = cost_calculator();
  }

  void perturbate(int n, std::tr1::mt19937& rng)
  {
    mets::perturbate(*this, n, rng);
  }

  friend std::ostream& operator<<(std::ostream& os, const atsp_model& p);
  friend std::istream& operator>>(std::istream& is, atsp_model& p);
  
protected:
  
  // Straight cost calculator
  int64_t cost_calculator() const
  {
    int64_t sum = 0;
    // opens and closes on pi_m.size()
    sum += matrix[ pi_m.size() ] [ pi_m[0] ];
    for(unsigned int ii(0); ii != pi_m.size()-1; ++ii)
      sum += matrix[ pi_m [ ii ] ] [ pi_m [ ii + 1 ] ];
    sum += matrix[ pi_m[pi_m.size()-1] ] [ pi_m.size() ];
    return sum;
  }

};


/// @brief Generates a the full subsequence inversion neighborhood.
class three_opt_full_neighborhood : public mets::move_manager
{
public:
  three_opt_full_neighborhood(int size) : move_manager()
  {
    for(int ii(0); ii!=size; ++ii)
      for(int jj(0); jj!=size; ++jj)
	for(int kk(0); kk!=size; ++kk)
	  if(ii != jj && ii != kk && jj != kk)
	    {
	      mets::complex_mana_move* m = new mets::complex_mana_move();
	      m->push_back(new mets::invert_subsequence(ii,jj));
	      m->push_back(new mets::invert_subsequence(kk, (kk+1)%size));
	      moves_m.push_back(m);
	    }
  } 
  
  /// @brief Dtor.
  ~three_opt_full_neighborhood() { }
  
  /// @brief Selects a different set of moves at each iteration.
  void refresh(mets::feasible_solution& s) { }
  
};
//________________________________________________________________________

// Input/Output functions
std::ostream&
operator<<(std::ostream& os, const atsp_model& atsp)
{
  for(unsigned int ii = 0; ii != atsp.pi_m.size(); ++ii) 
    os << (atsp.pi_m[ii]+1) << " ";
  return os;
}

std::istream&
operator>>(std::istream& is, atsp_model& atsp)
{
  std::string cmd, name, type, comment, edge_weight_type,
    edge_weight_format;
  int dimension;

  is >> cmd;
  while(cmd != "EDGE_WEIGHT_SECTION" && cmd != "EOF")
    {
      if(cmd == "NAME:")
	std::getline(is, name);
      else if(cmd == "TYPE:")
	is >> type;
      else if(cmd == "COMMENT:")
	std::getline(is, comment);
      else if(cmd == "DIMENSION:")
	is >> dimension;
      else if(cmd == "EDGE_WEIGHT_TYPE:")
	is >> edge_weight_type;
      else if(cmd == "EDGE_WEIGHT_FORMAT:")
	is >> edge_weight_format;
      else
	{
	  std::getline(is, edge_weight_format);
	  std::cerr << "Unrecognised command: " << cmd << std::endl;
	  ::exit(-1);
	}
      is >> cmd;
  }
  
  if(edge_weight_format != "FULL_MATRIX" || 
     edge_weight_type != "EXPLICIT" || type != "ATSP")
    {
      std::cerr << "Unsupported matrix/problem type." 
		<< edge_weight_format << " " 
		<< edge_weight_type << " " 
		<< type 
		<< std::endl;
      ::exit(-1);
    }
  if(cmd == "EDGE_WEIGHT_SECTION")
    {
      atsp.pi_m.resize(dimension-1);
      std::generate(atsp.pi_m.begin(), atsp.pi_m.end(), mets::sequence(0));
      atsp.matrix.resize(dimension);
      for(int ii(0); ii!=dimension; ++ii)
	{
	  atsp.matrix[ii].resize(dimension);
	  for(int jj(0); jj!=dimension; ++jj)
	    {
	      is >> atsp.matrix[ii][jj];
	    }
	}
      is >> cmd;
      if(cmd != "EOF")
	{
	  std::cerr << "Error reading full matrix" << std::endl;
	  ::exit(-1);
	}
    }
  return is;
}
