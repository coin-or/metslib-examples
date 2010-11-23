#include <iostream>
#include <cstdlib>
#include <fstream>
#include <vector>

#include <metslib/mets.hh>

#include "qap_model.hpp"

using namespace std;

void usage()
{
  cerr << "qap qaplib.dat" << endl;
  ::exit(1);
}

typedef mets::swap_neighborhood<std::tr1::mt19937> swap_neighborhood_t;

struct logger : public mets::search_listener<swap_neighborhood_t>
{
  explicit
  logger(std::ostream& o) 
    : mets::search_listener<swap_neighborhood_t>(), 
      iteration(0), 
      os(o) 
  { }
  
  void 
  update(mets::abstract_search<swap_neighborhood_t>* as) 
  {
    const mets::feasible_solution& p = as->working();
    if(as->step() == mets::abstract_search<swap_neighborhood_t>::MOVE_MADE)
      {
	os << iteration++ << " " 
	   << static_cast<const mets::evaluable_solution&>(p).cost_function() << "\n";
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

  // user defined problem
  qap_model problem_instance;
  in >> problem_instance;
  unsigned int N = problem_instance.size();

  // best solution instance used to record the best known solution.
  qap_model incumbent_solution(problem_instance);
  mets::best_ever_solution incumbent_recorder(incumbent_solution);

  // A neighborhood made of 2N random swaps
  swap_neighborhood_t neighborhood(rng, sqrt(N)*N);

  // generate a random starting point
  mets::random_shuffle(problem_instance, rng);

  // use framework provided strategies
  mets::simple_tabu_list tabu_list(N*sqrt(N));
  mets::best_ever_criteria aspiration_criteria;
      
  // fixed number of non improving moves before termination
  mets::noimprove_termination_criteria termination_criteria(1000);
	  
  // the search algorithm
  mets::tabu_search<swap_neighborhood_t> algorithm(problem_instance, 
						   incumbent_recorder, 
						   neighborhood, 
						   tabu_list, 
						   aspiration_criteria, 
						   termination_criteria);
  
  // log to standard error
  logger g(clog);
  algorithm.attach(g);
  algorithm.search();
	  
  // write solution to standard output
  cout << fixed << N << " " <<  incumbent_solution.cost_function() << endl
       << incumbent_solution << endl;
}
