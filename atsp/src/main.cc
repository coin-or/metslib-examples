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

  std::tr1::uniform_int<int> tlg(3, N);

  // best solution instance for recording
  atsp_model best_solution(problem_instance);

  cout << best_solution.cost_function() << " " << best_solution << endl;
  cout << problem_instance.cost_function() << " " << problem_instance << endl;
  
  // user defined neighborhood generator 
  mets::swap_neighborhood<std::tr1::mt19937> neighborhood(rng, N*(N-1)/4, N*2);

  // log to standard error
  logger g(clog);

  for(unsigned int starts = 0; starts != 5; ++starts) {
    // generate a random starting point
    problem_instance.random_shuffle(rng);
    
    // use framework provided strategies
    mets::simple_tabu_list tabu_list(tlg(rng));
    mets::best_ever_criteria aspiration_criteria;
    
    for(int iteration=0; iteration!=10; ++iteration) {
      // random tabu list tenure
      tabu_list.tenure(tlg(rng));

      // fixed number of non improving moves before termination
      mets::noimprove_termination_criteria 
	termination_criteria(100);
      
      // the search algorithm
      mets::tabu_search algorithm(problem_instance, 
				  best_solution, 
				  neighborhood, 
				  tabu_list, 
				  aspiration_criteria, 
				  termination_criteria);
      
      algorithm.attach(g);
      std::cout << "New iteration with tenure: " 
		<< tabu_list.tenure() << std::endl;
      algorithm.search();
      
      // restore best solution
      problem_instance = best_solution;
      // perturbate point with random swaps
      problem_instance.perturbate(N/10, rng);
    }
    
    cout << "Best so far: " << best_solution.cost_function()  << endl;
  }
  // write solution to standard output
  cout << best_solution.size() << " " 
       << best_solution.cost_function()  << endl;
  cout << best_solution << endl;
}
