#include <iostream>
#include <cstdlib>
#include <fstream>
#include <vector>

#include <metslib/mets.h>

#include "atsp_model.hpp"

using namespace std;

void usage()
{
  cerr << "atsp tsplib.dat" << endl;
  ::exit(1);
}

struct logger : public mets::search_listener
{
  explicit
  logger(std::ostream& o) 
    : mets::search_listener(), 
      iteration(0), 
      os(o) 
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


int main(int argc, char* argv[]) 
{

  if(argc != 2) usage();
  ifstream in(argv[1]);
  if(!in.is_open()) usage();

  // random number generator from C++ TR1 extension
  std::tr1::mt19937 rng(time(NULL));

  // user define problem
  atsp_model problem_instance;

  // read problem instance from standard input (no check is made)
  in >> problem_instance;

  unsigned int N = problem_instance.size();

  std::tr1::uniform_int<int> tlg(N/2, 3*N);

  // best ever solution 
  atsp_model optimum(problem_instance);

  // Neighborhood made of all possibile subsequence inversions.
  // It's the 2-opt neighborhood
  std::vector<mets::move_manager*> neighborhoods;
  neighborhoods.push_back(new mets::invert_full_neighborhood(N));
  neighborhoods.push_back(new three_opt_full_neighborhood(N));

  // log to standard error
  logger g(clog);

  for(unsigned int starts = 0; starts != 3; ++starts) {
    // generate a random starting point
    problem_instance.random_shuffle(rng);

    // best solution instance (records the best solution of each iteration)
    atsp_model major_best_solution(problem_instance);

    // Do minor iterations with a max no-improve criterion
    mets::noimprove_termination_criteria minor_it_criteria(8);
    while(!minor_it_criteria(major_best_solution))
      {

      // best solution instance (records the best solution of each iteration)
      atsp_model minor_best_solution(problem_instance);

      for(int ii=0; ii!=neighborhoods.size(); ++ii)
	{
	  std::cout << " * Run: " << starts+1 <<
	    "/" << ii << std::flush;
	  // the search algorithm
	  mets::local_search algorithm(problem_instance, 
				       minor_best_solution, 
				       *(neighborhoods[ii]),
				       true);
	  algorithm.attach(g);
	  algorithm.search();
	  if(minor_best_solution.cost_function() 
	     < major_best_solution.cost_function())
	    major_best_solution = minor_best_solution;
	  
	  if(major_best_solution.cost_function() 
	     < optimum.cost_function())
	    optimum = major_best_solution;
	  
	  std::cout << " Cost: " << minor_best_solution.cost_function() 
		    << "/" << major_best_solution.cost_function() 
		    << std::endl;
      }

      
      problem_instance = major_best_solution;
      // perturbate point with random swaps
      problem_instance.perturbate(N/3, rng);
    }
    
    cout << "Best of this run/so far: " 
	 << major_best_solution.cost_function() 
	 << "/" 
	 << optimum.cost_function()  << endl;
  }
  cout << "Best ever: " << optimum.cost_function()  << endl;
  // write solution to standard output
  cout << optimum.size() << " " 
       << optimum.cost_function()  << endl;
  cout << optimum << endl;
}
