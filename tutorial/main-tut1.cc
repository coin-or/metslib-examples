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

// Include metslib from the metslib subfolder
#include <metslib/mets.hh>

// Include the tutorial model (and neighborhood) definitions
#include "tut_model.h"
#include "tut_moves.h"
#include "tut_neighborhoods.h"

using namespace std;

/// @brief Generic progress printer.
///
/// This is actually an observer of the algorithm that receives
/// notification on the algorithm advancement. 
///
/// Update is called each time the algorithm has something to say us:
/// a move was made, the aspiration criteria was used, an improvement
/// was achieved and so on.
///
struct logger : public mets::search_listener<full_neighborhood> 
{
  explicit
  logger(std::ostream& o) 
    : mets::search_listener<full_neighborhood> (), iteration(0), os(o) 
  { }
  
  void 
  update(mets::abstract_search<full_neighborhood>  * as) 
  {
    const subsetsum& ss = static_cast<const subsetsum&>(as->working());
    if(as->step() == mets::abstract_search<full_neighborhood>::MOVE_MADE)
      {
        os  << iteration++ << ": " << ss.cost_function() << "/" << ss << "\n";
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

  // Search the subset of v such that it sums up to a number as near
  // to 109 as possible.
  // subsetsum is derived from mets::feasible_solution in tut_model.h 
  subsetsum model(v, 109);

  // storage for the best known solution.
  subsetsum best(model);

  // our neighborhood generator the neighborhood is explored using the
  // following instance derived from mets::move_manager in
  // tut_neighborhoods.h. Each element of the neighborhood is an
  // instance of a subclass mets::mana_move defined in tut_moves.h.
  full_neighborhood neigh(model.size());

  // progress logger to be attached to the algorithm (defined previously)
  logger g(clog);

  // We are done defining the model in the metslib framework, now we
  // can use the toolkit provided classes to try solve our problem.

  // simple tabu list (recency on moves)
  mets::simple_tabu_list tabu_list(7);
  // simple aspiration criteria
  mets::best_ever_criteria aspiration_criteria;

  // stop searching when not improving for 100 times
  mets::noimprove_termination_criteria noimprove(200);
  // chain the previous termination criteria with a threshold (when we
  // reach 0 we have solved our problem). The resulting
  // threshold_noimprove criterion terminate either when the objective
  // function reaches 0 or after 100 non improving moves.
  mets::threshold_termination_criteria 
    threshold_noimprove(&noimprove, 0);

  mets::best_ever_solution best_recorder(best);
  // Create a tabu_search algorithm instance starting from "model",
  // recording the best solution in "best", exploring the neighborhood
  // using "neigh", using the tabu list "tabu_list", the best ever
  // aspiration criteria "aspiration_criteria" and the combined
  // termination criteria "threshold_noimprove".
  mets::tabu_search<full_neighborhood> algorithm(model, 
						 best_recorder, 
						 neigh, 
						 tabu_list, 
						 aspiration_criteria, 
						 threshold_noimprove);
  algorithm.attach(g);
  algorithm.search();
  cout << "Best solution: " << best_recorder.best_ever().cost_function()  << endl;
  cout << (const subsetsum&)best_recorder.best_ever() << endl;

  return 0;
}
