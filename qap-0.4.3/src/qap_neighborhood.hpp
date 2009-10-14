#pragma once
#include <metslib/mets.h>

#include "qap_move.hpp"

/// @brief Generates a stochastic subset of the neighborhood.
template<typename random_generator = std::tr1::minstd_rand0>
class qap_neighborhood : public mets::move_manager
{
public:  
  qap_neighborhood(random_generator& r, 
		   unsigned int moves, 
		   unsigned int complex_moves)
    : mets::move_manager(), rng(r), int_range(0), n(moves), nc(complex_moves)
  { 
    // n simple moves
    for(unsigned int ii = 0; ii != n; ++ii) 
      moves_m.push_back(new qap_move(0,0));

    // nc double moves
    for(unsigned int ii = 0; ii != nc; ++ii) 
      {
	mets::complex_mana_move& cm = *new mets::complex_mana_move(2);
	cm[0] = new qap_move(0,0);
	cm[1] = new qap_move(0,0);
	moves_m.push_back(&cm);
      }

  }  

  ~qap_neighborhood()
  {
    // delete all moves
    for(iterator ii = begin(); ii != end(); ++ii)
      delete (*ii);
  }

  void refresh(mets::feasible_solution& s)
  {
    qap_model& sol = dynamic_cast<qap_model&>(s);
    iterator ii = begin();

    // the first n are simple qap_moveS
    for(unsigned int cnt = 0; cnt != n; ++cnt)
      {
	qap_move* m = static_cast<qap_move*>(*ii);
	randomize_move(*m, sol.size());
	++ii;
      }

    // the following nc are complex_mana_moves made of 2 qap_moveS
    for(unsigned int cnt = 0; cnt != nc; ++cnt)
      {
	mets::complex_mana_move& cm = 
	  static_cast<mets::complex_mana_move&>(**ii);
	for(int jj = 0; jj != 2; ++jj)
	  {
	    qap_move* m = static_cast<qap_move*>(cm[jj]);
	    randomize_move(*m, sol.size());
	  }
	++ii;
      }

  }

protected:
  random_generator& rng;
  std::tr1::uniform_int<> int_range;
  unsigned int n;
  unsigned int nc;

  void randomize_move(qap_move& m, unsigned int size)
  {
    int p1 = int_range(rng, size);
    int p2 = int_range(rng, size);
    while(p1 == p2) 
      p2 = int_range(rng, size);
    // we are friend, so we know how to 
    // handle nuts&bolts of qap_moves
    m.p1 = std::min(p1,p2); 
    m.p2 = std::max(p1,p2); 
  }

};
