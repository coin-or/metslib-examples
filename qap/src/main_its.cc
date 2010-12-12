#include <iostream>
#include <cstdlib>
#include <fstream>
#include <vector>

#include <metslib/mets.hh>

#include "qap_model.hpp"

using namespace std;

/// a typedef for the neighbourhood we will use
typedef mets::swap_neighborhood<std::tr1::mt19937> neighborhood_t;

// print usage message and exit program
void usage()
{
  cerr << "qap qaplib.dat" << endl;
  ::exit(1);
}

// the main function implementing the ITS algorithm 
//
// the algorithm does multistart, at each start the solution is
// iteratively improved using a TS strategy and perturbated by random
// swaps (using an Iterated Local Search strategy).
//  
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

  assert(N>10);

  std::tr1::uniform_int<int> tlg(5, N*sqrt(N)/2);
  std::tr1::uniform_int<int> psg(N/5, N/2);

  // best solution instance for recording
  // storage for the best known solution.
  qap_model ils_solution(problem_instance);
  mets::best_ever_solution ils_recorder(ils_solution);

  // A neighborhood made of random swaps
  neighborhood_t neighborhood(rng, N*sqrt(N));

  // log to standard error
  clog << fixed;
  mets::improvement_logger<neighborhood_t> g(clog);

  // generate a random starting point
  mets::random_shuffle(problem_instance, rng);

  // Do minor iterations with a max no-improve criterion
  mets::noimprove_termination_criteria 
    ils_stop(200);
  
  int perturbation_size = N/2;

  while(!ils_stop(ils_recorder.best_seen())) 
    {

      // best solution instance for recording storage for the best
      // known solution of the ils iteration.
      qap_model ts_solution(problem_instance);
      mets::best_ever_solution ts_recorder(ts_solution);

      // use framework provided strategies
      mets::simple_tabu_list tabu_list(tlg(rng));
      mets::best_ever_criteria aspiration_criteria;
      
      // random tabu list tenure
      tabu_list.tenure(tlg(rng));
	  
      // fixed number of non improving moves before termination
      mets::noimprove_termination_criteria 
	ts_stop(750);
      
      // the search algorithm
      mets::tabu_search<neighborhood_t> algorithm(problem_instance, 
						  ts_recorder, 
						  neighborhood, 
						  tabu_list, 
						  aspiration_criteria, 
						  ts_stop);
	  
      algorithm.attach(g);
      algorithm.search();	
      ils_recorder.accept(ts_recorder.best_seen());
      problem_instance.copy_from(ils_recorder.best_seen());
      perturbation_size = psg(rng);
      // perturbate point with random swaps
      // mets::random_shuffle(problem_instance, rng);
      mets::perturbate(problem_instance, perturbation_size, rng);
    }
  // write solution to standard output
  cout << fixed << N << " " <<  ils_solution.cost_function() << endl
       << ils_solution << endl;
}
