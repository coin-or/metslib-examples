// METSlib tutorial source file - main-tut1.cc                 -*- C++ -*-
//
// Copyright (C) 2009 Mirko Maischberger <mirko.maischberger@gmail.com>
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <iostream>
#include <metslib/mets.h>
#include "tut_model.hh"
#include "tut_moves.hh"
#include "tut_neighborhoods.hh"

using namespace std;

/// @brief Generic progress printer.
struct logger : public mets::search_listener
{
  explicit
  logger(std::ostream& o) 
    : mets::search_listener(), iteration(0), os(o) 
  { }
  
  void 
  update(mets::abstract_search* as) 
  {
    const mets::feasible_solution& p = as->working();
    if(as->step() == mets::abstract_search::MOVE_MADE)
      {
        os << iteration++ << " " << p.cost_function() << "\n";
      }
  }
  
protected:
  int iteration;
  ostream& os;
};

int main() 
{
  
  // simple tabu search for the subset sum problem

  // the set
  int numbers[] = { 475, 382, -202, 351, 296, -362, 336, 117, -319, 
		    416, -304, 364, -386, -9, 391, 389, -457, 261, 
		    -323, -498, 407, -81, 445, -308, 258, -274, 156 };
  vector<int> v(&numbers[0], &numbers[sizeof(numbers)/sizeof(int)]);

  // the starting point 
  subsetsum model(v, 109);
  // the best known solution
  subsetsum best(model);
  // our neighborhood generator
  full_neighborhood neigh(model.size());
  // progress logger to be attached to the algorithm
  logger g(clog);
  // simple tabu list (recency on moves)
  mets::simple_tabu_list tabu_list(5);
  // simple aspiration criteria
  mets::best_ever_criteria aspiration_criteria;
  // stop searching when not improving for 100 times
  mets::noimprove_termination_criteria 
    noimprove(100);
  // chain the previous termination criteria with a threshold (when we
  // reach 0 we have solved our problem)
  mets::threshold_termination_criteria 
    threshold_noimprove(&noimprove, 0);
  mets::tabu_search algorithm(model, 
			      best, 
			      neigh, 
			      tabu_list, 
			      aspiration_criteria, 
			      threshold_noimprove);
  algorithm.attach(g);
  algorithm.search();
  cout << "Best solution: " << best.cost_function()  << endl;
  cout << best << endl;

  return 0;
}
