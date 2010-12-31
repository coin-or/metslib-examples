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

typedef mets::swap_full_neighborhood swap_neighborhood_t;

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
  mets::best_ever_solution best_recorder(incumbent_solution);

  // A neighborhood made of 2N random swaps
  swap_neighborhood_t neighborhood(N);

  mets::iteration_termination_criteria rrls_stop(200);
  mets::iteration_logger<swap_neighborhood_t> g1(clog);
  mets::improvement_logger<swap_neighborhood_t> g2(clog);
  while(!rrls_stop(best_recorder.best_seen())) {
    // generate a random starting point
    mets::random_shuffle(problem_instance, rng);
    // the search algorithm
    mets::local_search<swap_neighborhood_t> algorithm(problem_instance, 
						      best_recorder, 
						      neighborhood);
    // log to standard error
    algorithm.attach(g1);
    algorithm.attach(g2);
    algorithm.search();
  }
	  
  // write solution to standard output
  cout << fixed << N << " " <<  incumbent_solution.cost_function() << endl
       << incumbent_solution << endl;
}
