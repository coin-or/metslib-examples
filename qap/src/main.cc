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

typedef mets::swap_neighborhood<std::tr1::mt19937> neighborhood_t;

struct logger : public mets::search_listener<neighborhood_t>
{
  explicit
  logger(std::ostream& o) 
    : mets::search_listener<neighborhood_t>(), 
      iteration(0), 
      os(o) 
  { }
  
  void 
  update(mets::abstract_search<neighborhood_t>* as) 
  {
    const mets::feasible_solution& p = as->working();
    if(as->step() == mets::abstract_search<neighborhood_t>::MOVE_MADE)
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

  // user define problem
  qap_model problem_instance;

  // read problem instance from standard input (no check is made)
  in >> problem_instance;

  unsigned int N = problem_instance.size();

  std::tr1::uniform_int<int> tlg(7, N*7);
  std::tr1::uniform_int<int> psg(7, N/2);

  // best solution instance for recording
  // storage for the best known solution.
  qap_model incumbent_solution(problem_instance);
  mets::best_ever_solution incumbent_recorder(incumbent_solution);

  // A neighborhood made of random swaps
  neighborhood_t
    neighborhood(rng, N*12);

  // log to standard error
  logger g(clog);

  for(unsigned int starts = 0; starts != int(sqrt(N)); ++starts) 
    {
      // generate a random starting point
      mets::random_shuffle(problem_instance, rng);

      // best solution instance for recording storage for the best
      // known solution of the major iteration.
      qap_model majorit_solution(problem_instance);
      mets::best_ever_solution majorit_recorder(majorit_solution);

      // use framework provided strategies
      mets::simple_tabu_list tabu_list(tlg(rng));
      mets::best_ever_criteria aspiration_criteria;
      
      // Do minor iterations with a max no-improve criterion
      mets::noimprove_termination_criteria 
	minor_it_criteria(100);
      
      int perturbation_size = N;

      while(!minor_it_criteria(majorit_recorder.best_seen())) 
	{
	  
	  // best solution instance for recording storage for the best
	  // known solution of the major iteration.
	  qap_model minorit_solution(problem_instance);
	  mets::best_ever_solution minorit_recorder(minorit_solution);
	  
	  // random tabu list tenure
	  tabu_list.tenure(tlg(rng));
	  
	  // fixed number of non improving moves before termination
	  mets::noimprove_termination_criteria 
	    termination_criteria(200);
	  
	  // the search algorithm
	  mets::tabu_search<neighborhood_t> algorithm(problem_instance, 
						      minorit_recorder, 
						      neighborhood, 
						      tabu_list, 
						      aspiration_criteria, 
						      termination_criteria);
	  
	  algorithm.attach(g);
	  std::cout << "New iteration with tenure: " 
		    << tabu_list.tenure() << std::endl;

	  algorithm.search();
	  
	  majorit_recorder.accept(minorit_recorder.best_seen());
	  problem_instance.copy_from(majorit_recorder.best_seen());


	  perturbation_size = psg(rng);
	  // perturbate point with random swaps
	  mets::perturbate(problem_instance, perturbation_size, rng);
	}
      
      incumbent_recorder.accept(majorit_recorder.best_seen());
      
      cout << "Best of this run/so far: " 
	   << majorit_solution.cost_function()  
	   << "/"
	   << incumbent_solution.cost_function() << endl;
      
    }
  // write solution to standard output
  cout << N << " " <<  incumbent_solution.cost_function() << endl
       << incumbent_solution << endl;
}
