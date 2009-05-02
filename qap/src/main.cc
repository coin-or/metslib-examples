#include <iostream>
#include <cstdlib>
#include <fstream>
#include <vector>

#include <metslib/mets.h>

#include "qap_model.hpp"
#include "qap_move.hpp"
#include "qap_neighborhood.hpp"

using namespace std;

void usage()
{
  cerr << "qap qaplib.dat" << endl;
  ::exit(1);
}

struct grapher : public mets::search_listener
{
  explicit
  grapher(std::ostream& o) 
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
	iteration++;
      }
    if(as->step() == mets::abstract_search::IMPROVEMENT_MADE)
      {
      	os << iteration << ": " << p.cost_function() << std::endl;
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

  // user define problem
  qap_model problem_instance;

  // read problem instance from standard input (no check is made)
  in >> problem_instance;

  unsigned int N = problem_instance.size();

  std::tr1::mt19937 rng(time(NULL));
  problem_instance.random_shuffle<std::tr1::mt19937>(rng);

  qap_model best_solution(problem_instance);

  // user defined neighborhood generator 
  // (give specific random generator)
  qap_neighborhood<std::tr1::mt19937> neighborhood(rng, N*12,N*6);

  std::tr1::uniform_int<int> tlg(1, 150);
  
  // framework provided strategies
  mets::simple_tabu_list tabu_list(tlg(rng));
  mets::best_ever_criteria aspiration_criteria;
  
  grapher g(clog);
  for(int run=0; run!=12; ++run)
    {
      tabu_list.tenure(tlg(rng));
      mets::noimprove_termination_criteria 
	termination_criteria(2000);
      
      // the search algorithm
      mets::tabu_search algorithm(problem_instance, 
				  best_solution, 
				  neighborhood, 
				  tabu_list, 
				  aspiration_criteria, 
				  termination_criteria);
      
      algorithm.attach(g);
      std::clog << "New iteration with tenure: " 
		<< tabu_list.tenure() << std::endl;
      algorithm.search();
      problem_instance = best_solution;
      problem_instance.perturbate(N/5, rng);
    }
  cout << best_solution.size() << " " 
       << best_solution.cost_function()  << endl;
  cout << best_solution << endl;
  
}
