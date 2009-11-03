#include <vector>
#include <fstream>
#include <iostream>

#include <metslib/mets.h>
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/random.hpp>

int g_colors;

class vcp_neighborhood2;

class vcp: public mets::feasible_solution {

public:
  typedef boost::adjacency_list<boost::vecS, 
				boost::vecS, 
				boost::bidirectionalS> graph_type;
  typedef graph_type::vertex_iterator vertex_iterator;
  typedef graph_type::out_edge_iterator out_edge_iterator;
  typedef graph_type::in_edge_iterator in_edge_iterator;
  typedef graph_type::edge_iterator edge_iterator;
  
  vcp(int nodes, const std::vector< std::pair<int, int> >& edges) 
    : cost_m(0), g_m(edges.begin(), edges.end(), nodes), 
      color_m(nodes) 
  { update_cost(); }

  mets::gol_type cost_function() const
  {
    // update_cost();
    return cost_m;
  }

  size_t size() { return boost::num_vertices(g_m); }

  void copy_from(const mets::feasible_solution& other)
  {
    const vcp& o = static_cast<const vcp&>(other);
    cost_m = o.cost_m;
    g_m = o.g_m;
    color_m = o.color_m;
  }

  int color(int i) const 
  { return color_m[i]; }

  int color(int i, int c) 
  { 
    int oldc = color_m[i]; 

    if(oldc == c) return c;
    {
      out_edge_iterator ei, ee;
      for(boost::tie(ei, ee) = boost::out_edges(i, g_m); ei!=ee; ++ei)
	{
	  if(oldc == color_m[boost::target(*ei,g_m)])
	    {
	      --cost_m;
	    }
	  else if(c == color_m[boost::target(*ei,g_m)])
	    {
	      ++cost_m;
	    }
	}
    }
    {
      in_edge_iterator ei, ee;
      for(boost::tie(ei, ee) = boost::in_edges(i, g_m); ei!=ee; ++ei)
	{
	  if(oldc == color_m[boost::source(*ei,g_m)])
	    {
	      --cost_m;
	    }
	  else if(c == color_m[boost::source(*ei,g_m)])
	    {
	      ++cost_m;
	    }
	}
    }
    color_m[i] = c; 
    return oldc; 
  }

  template<typename generator>
  void randomize(int colors, generator& gen)
  {
    boost::uniform_int<> color_dist(0, colors-1);
    boost::variate_generator<generator&, boost::uniform_int<> >
      colorgen(gen, color_dist);
    for(std::vector<int>::iterator ii = color_m.begin(); 
	ii != color_m.end(); 
	++ii)
      {
	*ii = colorgen();
      }
    update_cost();
  }

  template<typename generator>
  void perturbate(int colors, int qty, generator& gen)
  {
    boost::uniform_int<> color_dist(0, colors-1);
    boost::uniform_int<> node_dist(0, boost::num_vertices(g_m)-1);
    boost::variate_generator<generator&, boost::uniform_int<> >
      colorgen(gen, color_dist);
    boost::variate_generator<generator&, boost::uniform_int<> >
      nodegen(gen, node_dist);
    
    for(int ii(0); ii!=qty; ++ii)
      {
	color_m[nodegen()] = colorgen();
      }
    update_cost(); 
  }

  void print(std::ostream& os)
  {
    for(int ii(0); ii!=boost::num_vertices(g_m); ++ii)
      {
	os << color_m[ii] << "\n";
      }
  }

  void print_edges(std::ostream& os)
  {
    edge_iterator ei, ee;
    for(boost::tie(ei, ee) = boost::edges(g_m); ei!=ee; ++ei)
      {
	os << "e " 
	   << boost::source(*ei, g_m) << " " 
	   << boost::target(*ei, g_m) << " ("
	   << color_m[boost::source(*ei,g_m)] << " " 
	   << color_m[boost::target(*ei,g_m)] << ")\n";
      }
  }

  void display_conflicts(std::ostream& os)
  {
    vertex_iterator vi, ve;
    cost_m = 0;
    edge_iterator ei, ee;
    for(boost::tie(ei, ee) = boost::edges(g_m); ei!=ee; ++ei)
      {
	if(color_m[boost::source(*ei,g_m)] == 
	   color_m[boost::target(*ei,g_m)])
	  {
	    os << boost::source(*ei,g_m) << " " 
	       << boost::target(*ei,g_m) << " ("
	       << color_m[boost::source(*ei,g_m)] << ")"
	       << std::endl;
	  }
      } 
  }

protected:
  mutable int cost_m;
  graph_type g_m;
  std::vector<int> color_m;

  friend class vcp_neighborhood;

  void update_cost() const
  {
    cost_m = 0;
    edge_iterator ei, ee;
    for(boost::tie(ei, ee) = boost::edges(g_m); ei!=ee; ++ei)
      {
	if(color_m[boost::source(*ei,g_m)] == 
	   color_m[boost::target(*ei,g_m)])
	  {
	    cost_m++;
	  }
      }
  }

};

class vcp_set : public mets::mana_move
{
public:

  vcp_set(int i, int c) : i_m(i), c_m(c) {}

  void set(int i, int c) 
  {i_m = i; c_m = c;}

  int node() const
  { return i_m; }

  int color() const
  { return c_m; }

  void apply(mets::feasible_solution& sol)
  {
    vcp& v = static_cast<vcp&>(sol);
    oldc_m = v.color(i_m, c_m);
  }

  void unapply(mets::feasible_solution& sol)
  {
    vcp& v = static_cast<vcp&>(sol);
    v.color(i_m, oldc_m);
  }

  mets::mana_move* clone() const { return new vcp_set(i_m, c_m); }

  bool operator==(const mets::mana_move& other) const
  {
    const vcp_set& v = static_cast<const vcp_set&>(other);
    return (v.i_m == i_m && v.c_m == c_m);
  }

  size_t hash() const 
  { return i_m << 7 | c_m; }

  int i_m;
  int c_m;
  int oldc_m;
};

class vcp_neighborhood : public mets::move_manager
{
public:
  vcp_neighborhood(int n, int ki) : nodes_m(n), ki_m(ki) 
  {
    for(int i(0); i!=nodes_m; ++i)
      for(int c(0); c!=ki_m-1; ++c)
	moves_m.push_back(new vcp_set(i, 0));
  }
  
  void refresh(mets::feasible_solution& other)
  { 
    vcp& v = static_cast<vcp&>(other);
    iterator mi = moves_m.begin();
    for(int i(0); i!=nodes_m; ++i)
      for(int c(0); c!=ki_m; ++c)
	{
	  vcp_set& m = static_cast<vcp_set&>(**mi);
	  if(c == v.color(i)) ++c;
	  if(c == ki_m) break;
	  m.set(m.node(), c);
	  ++mi;
	}
  }
  int nodes_m;
  int ki_m;
};

class vcp_neighborhood2 : public mets::move_manager
{
public:
  vcp_neighborhood2(boost::mt19937& rng, int ki, int sm, int dm) : 
    rng_m(rng), ki_m(ki), sm_m(sm), dm_m(dm)
  {
    for(int ii(0); ii!=sm_m; ++ii)
      {
	moves_m.push_back(new vcp_set(0, 0));
      }
    for(int ii(0); ii!=dm_m; ++ii)
      {
	mets::complex_mana_move& m = *new mets::complex_mana_move(0);
	m.push_back(new vcp_set(0, 0));
	m.push_back(new vcp_set(0, 0));
	moves_m.push_back(&m);
      }
  }
  
  void refresh(mets::feasible_solution& other)
  { 
    vcp& v = static_cast<vcp&>(other);
    boost::uniform_int<> color_dist(0, ki_m);
    boost::uniform_int<> node_dist(0, v.size()-1);
    boost::variate_generator<boost::mt19937&, boost::uniform_int<> >
      colorgen(rng_m, color_dist);
    boost::variate_generator<boost::mt19937&, boost::uniform_int<> >
      nodegen(rng_m, node_dist);
    iterator m = moves_m.begin();
    for(int ii(0); ii!=sm_m; ++ii)
      {
	vcp_set* vm = dynamic_cast<vcp_set*>(*m);
	vm->set(nodegen(), colorgen());
	++m;
      }
    for(int ii(0); ii!=dm_m; ++ii)
      {
	mets::complex_mana_move& cm = 
	  dynamic_cast<mets::complex_mana_move&>(**m);
	dynamic_cast<vcp_set*>(cm[0])->set(nodegen(), colorgen());
	dynamic_cast<vcp_set*>(cm[1])->set(nodegen(), colorgen());
	++m;
      }
  }
  boost::mt19937 rng_m;
  int sm_m;
  int dm_m;
  int ki_m;
};

struct logger : public mets::search_listener
{
  explicit
  logger(std::ostream& o) 
    : mets::search_listener(), 
      iteration(0), 
      os(o), 
      gen(time(NULL))
  { }
  
  void 
  update(mets::abstract_search* as) 
  {
    vcp& v = static_cast<vcp&>(as->working());
    if(as->step() == mets::abstract_search::MOVE_MADE)
      {
	iteration++;
	os << v.cost_function() << "\t";
	if(iteration % 20 == 0) os << std::flush;
	if(iteration % 2500 == 0) v.perturbate(g_colors, v.size()/2, gen);
      }
  }
  
protected:
  int iteration;
  std::ostream& os;
  boost::mt19937 gen;
};


using namespace std;

int main(int argc, char* argv[])
{
  if(argc != 4) {
    clog << "vcp file.col colors tenure" << endl;
    exit(1);
  }
  std::ifstream in(argv[1]);
  g_colors = ::atoi(argv[2]);
  int tenure = ::atoi(argv[3]);
  vector< pair<int, int> > edges;
  int n, m;
  while(in.good())
    {
      char h;
      in >> h;
      if(!in.good()) break;
      switch(h)
	{
	case 'c':
	  {
	    string comment;
	    getline(in, comment);
	    clog << comment << endl;
	    break;
	  }
	case 'p':
	  {
	    string problem;
	    in >> problem;
	    if(problem != "edge")
	      {
		clog << "Unknown file format." << endl;
		exit(1);
	      }
	    in >> n >> m;
	    break;
	  }
	case 'e':
	  {
	    int v1,v2;
	    in >> v1 >> v2;
	    edges.push_back(std::make_pair(v1-1,v2-1));
	  }
	}
    }

  vcp point(n, edges);
  // storage for the best known solution.
  vcp incumbent(point);

  boost::mt19937 gen(time(NULL));
  point.randomize(g_colors, gen);

  // our neighborhood generator the neighborhood is explored using the
  // following instance derived from mets::move_manager in
  // tut_neighborhoods.h. Each element of the neighborhood is an
  // instance of a subclass mets::mana_move defined in tut_moves.h.
  g_colors -= 5;
  vcp_neighborhood2 neigh(gen, g_colors, n*g_colors/4, n*g_colors/8);
  // vcp_neighborhood neigh(n, g_colors);
  // We are done defining the model in the metslib framework, now we
  // can use the toolkit provided classes to try solve our problem.

  // simple tabu list (recency on moves)
  mets::simple_tabu_list tabu_list(tenure);
  // simple aspiration criteria
  mets::best_ever_criteria aspiration_criteria;

  // stop searching when not improving for N times
  mets::noimprove_termination_criteria noimprove(2000);
  // chain the previous termination criteria with a threshold (when we
  // reach 0 we have solved our problem). The resulting
  // threshold_noimprove criterion terminate either when the objective
  // function reaches 0 or after 100 non improving moves.
  mets::threshold_termination_criteria 
    threshold_noimprove(&noimprove, 0);

  // Create a tabu_search algorithm instance starting from "model",
  // recording the best solution in "best", exploring the neighborhood
  // using "neigh", using the tabu list "tabu_list", the best ever
  // aspiration criteria "aspiration_criteria" and the combined
  // termination criteria "threshold_noimprove".
  mets::tabu_search algorithm(point, 
			      incumbent, 
			      neigh, 
			      tabu_list, 
			      aspiration_criteria, 
			      threshold_noimprove);
  //true);
  // algorithm.attach(g);
  logger log(cout);
  algorithm.attach(log);
  algorithm.search();

  point = incumbent;

  clog << "Best solution: " << incumbent.cost_function()  << endl;

  g_colors += 5;
  vcp_neighborhood neigh2(n, g_colors);
  mets::local_search algo2(point, 
			  incumbent, 
			  neigh2, 
			  false);

  algo2.search();

  clog << "Best solution: " << incumbent.cost_function()  << endl;

  ofstream sol("solution.txt");
  incumbent.print(sol);
  sol.close();

  ofstream dbg("debug");
  incumbent.print_edges(dbg);
  dbg.close();

  incumbent.display_conflicts(clog);

  return 0;
}
