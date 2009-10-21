#include <iostream>
#include <cstdlib>
#include <fstream>
#include <vector>

#include <metslib/mets.h>

#include "qap_model.hpp"

using namespace std;

void usage()
{
  cerr << "qap qaplib.dat" << endl;
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
  qap_model problem_instance;

  // read problem instance from standard input (no check is made)
  in >> problem_instance;

  unsigned int N = problem_instance.size();

  std::tr1::uniform_int<int> tlg(7, N*7);

  // best solution instance for recording
  qap_model optimum(problem_instance);

  // A library provided neighborhood make of random swaps and double
  // swaps.
  mets::swap_neighborhood<std::tr1::mt19937> 
    neighborhood(rng, N*12, N*6);

  // log to standard error
  logger g(clog);

  for(unsigned int starts = 0; starts != N/5; ++starts) {
    // generate a random starting point
    problem_instance.random_shuffle(rng);

    // record the best solution of each major iteration
    qap_model major_best_solution(problem_instance);

    // use framework provided strategies
    mets::simple_tabu_list tabu_list(tlg(rng));
    mets::best_ever_criteria aspiration_criteria;
    
    // Do minor iterations with a max no-improve criterion
    mets::noimprove_termination_criteria 
      minor_it_criteria(20);

    while(!minor_it_criteria(major_best_solution)) {
      
      // record best solution of minor iteration
      qap_model minor_best_solution(problem_instance);

      // random tabu list tenure
      tabu_list.tenure(tlg(rng));

      // fixed number of non improving moves before termination
      mets::noimprove_termination_criteria 
	termination_criteria(750);
      
      // the search algorithm
      mets::tabu_search algorithm(problem_instance, 
				  minor_best_solution, 
				  neighborhood, 
				  tabu_list, 
				  aspiration_criteria, 
				  termination_criteria);
      
      algorithm.attach(g);
      std::cout << "New iteration with tenure: " 
		<< tabu_list.tenure() << std::endl;
      algorithm.search();
      
      if(minor_best_solution.cost_function() 
	 < major_best_solution.cost_function())
	major_best_solution = minor_best_solution;
      
      if(major_best_solution.cost_function() 
	 < optimum.cost_function())
	optimum = major_best_solution;

      problem_instance = major_best_solution;
      // perturbate point with N/5 random swaps
      problem_instance.perturbate(N/4, rng);
    }
    
    cout << "Best of this run/so far: " 
	 << major_best_solution.cost_function()  
	 << "/"
	 << optimum.cost_function() << endl;

  }
  // write solution to standard output
  cout << optimum.size() << " " 
       << optimum.cost_function()  << endl;
  cout << optimum << endl;
}
