#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>

#include <metslib/mets.hh>
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/random.hpp>

int g_colors;

class vcp_neighborhood;

class vcp: public mets::evaluable_solution {

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
      color_m(nodes), conflicts_m(nodes)
  { update_cost(); }

  mets::gol_type cost_function() const { return cost_m; }

  size_t size() { return boost::num_vertices(g_m); }

  void copy_from(const mets::copyable& other)
  {
    const vcp& o = static_cast<const vcp&>(other);
    cost_m = o.cost_m;
    g_m = o.g_m;
    color_m = o.color_m;
    conflicts_m = o.conflicts_m;
  }

  int conflicts(int i) const
  { return conflicts_m[i]; }

  int color(int i) const 
  { return color_m[i]; }

  double evaluate(int i, int c) const
  {
    double eval = cost_m;
    int oldc = color_m[i]; 
    {
      out_edge_iterator ei, ee;
      for(boost::tie(ei, ee) = boost::out_edges(i, g_m); ei!=ee; ++ei)
	{
	  if(oldc == color_m[boost::target(*ei,g_m)])
	    {
	      --eval;
	    }
	  else if(c == color_m[boost::target(*ei,g_m)])
	    {
	      ++eval;
	    }
	}
    }
    {
      in_edge_iterator ei, ee;
      for(boost::tie(ei, ee) = boost::in_edges(i, g_m); ei!=ee; ++ei)
	{
	  if(oldc == color_m[boost::source(*ei,g_m)])
	    {
	      --eval;
	    }
	  else if(c == color_m[boost::source(*ei,g_m)])
	    {
	      ++eval;
	    }
	}
    }
    return eval;
  }

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
	      --conflicts_m[boost::source(*ei,g_m)];
	      --conflicts_m[boost::target(*ei,g_m)];
	      --cost_m;
	    }
	  else if(c == color_m[boost::target(*ei,g_m)])
	    {
	      ++conflicts_m[boost::source(*ei,g_m)];
	      ++conflicts_m[boost::target(*ei,g_m)];
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
	      --conflicts_m[boost::source(*ei,g_m)];
	      --conflicts_m[boost::target(*ei,g_m)];
	      --cost_m;
	    }
	  else if(c == color_m[boost::source(*ei,g_m)])
	    {
	      ++conflicts_m[boost::source(*ei,g_m)];
	      ++conflicts_m[boost::target(*ei,g_m)];
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
    for(unsigned int ii(0); ii!=boost::num_vertices(g_m); ++ii)
      {
	os << color_m[ii] << " ";
      }
  }

  void print_dot(int ki, std::ostream& os)
  {
    os << "graph VCP { " << std::endl;
    for(unsigned int ii(0); ii!=boost::num_vertices(g_m); ++ii)
      {
	os << "  n" << ii << " [label=\"" 
	   << ii << "\",style=\"filled\",fillcolor=\"";
	os.fill('0');
	os << std::hex;
	os << "#" << std::setw(6) 
	   << (uint32_t)(((1+color_m[ii])*65087*131) % 0xffffff) << " ";
	os << std::dec;
	os.fill(' ');
	os << "\"];" << std::endl;
      }
    edge_iterator ei, ee;
    for(boost::tie(ei, ee) = boost::edges(g_m); ei!=ee; ++ei)
      {
	os << "  n" 
	   << boost::source(*ei, g_m) << " -- n" 
	   << boost::target(*ei, g_m) << " ";

	if(color_m[boost::source(*ei,g_m)] == 
	   color_m[boost::target(*ei,g_m)])
	  {
	    os << "[style=\"bold\"]";
	  }
	os << ";" << std::endl;
      }
    os << "}" << std::endl;
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
  mutable std::vector<int> conflicts_m;

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
	    ++conflicts_m[boost::source(*ei,g_m)];
	    ++conflicts_m[boost::target(*ei,g_m)];
	    ++cost_m;
	  }
      }
  }

};

class vcp_set : public mets::mana_move
{
public:

  vcp_set(int i, int c) : i_m(i), c_m(c) {}

  void set(int i, int c) 
  { i_m = i; c_m = c; }

  int node() const
  { return i_m; }

  int color() const
  { return c_m; }

  void apply(mets::feasible_solution& sol) const
  {
    vcp& v = static_cast<vcp&>(sol);
    v.color(i_m, c_m);
  }

  double evaluate(const mets::feasible_solution& sol) const
  {
    const vcp& v = static_cast<const vcp&>(sol);
    return v.evaluate(i_m, c_m);;
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
};

class vcp_neighborhood : public mets::move_manager
{
public:
  vcp_neighborhood(int ki) : ki_m(ki) 
  { }
  
  void refresh(mets::feasible_solution& other)
  { 
    for(iterator mi = moves_m.begin(); mi != moves_m.end(); ++mi)
      {
	delete *mi;
      }
    moves_m.clear();

    vcp& v = static_cast<vcp&>(other);
    for(unsigned int i(0); i!=v.size(); ++i)
      for(int c(0); c!=ki_m; ++c)
	{
	  if(v.conflicts(i) && v.color(i) != c)
	    moves_m.push_back(new vcp_set(i, c));
	}
  }
  int ki_m;
};

template<typename neighborhood_t>
struct logger : public mets::search_listener<neighborhood_t>
{
  explicit
  logger(std::ostream& o) 
    : mets::search_listener<neighborhood_t>(), 
      iteration(0), 
      os(o), 
      gen(time(NULL))
  { }
  
  void 
  update(mets::abstract_search<neighborhood_t>* as) 
  {
    vcp& v = static_cast<vcp&>(as->working());
    if(as->step() == mets::abstract_search<neighborhood_t>::MOVE_MADE)
      {
	iteration++;
	v.print(os);
	os << v.cost_function() << '\n';
	if(iteration % 20 == 0) os << std::flush;
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
  bool fast = false;

  if(argc == 4)
    fast = true;
  else if(argc != 3) 
    {
      clog << "vcp file.col colors" << endl;
      exit(1);
    }

  std::ifstream in(argv[1]);
  int num_colors = ::atoi(argv[2]);
  vector< pair<int, int> > edges;
  int n = 0, m = 0;
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

  assert(n);
  assert(m);
  assert( m == edges.size() );

  boost::mt19937 gen(time(NULL));

  g_colors = num_colors;

  vcp point(n, edges);
  point.randomize(g_colors, gen);

  // storage for the best known solution.
  vcp best(n, edges);
  mets::best_ever_solution best_recorder(best);
  
  boost::uniform_int<> tl_dist(1, n);
  boost::variate_generator<boost::mt19937&, boost::uniform_int<> >
    tlg(gen, tl_dist);

  vcp_neighborhood neigh(g_colors);

  ofstream flog((string(argv[1]) + ".log").c_str());
  logger<vcp_neighborhood> log(flog);

  if(!fast) for(int run=0; run!=10; ++run)
    {
      point.randomize(g_colors, gen);

      vcp major_store(n, edges);
      mets::best_ever_solution major_best(major_store);

      // simple tabu list (recency on moves)
      mets::simple_tabu_list tabu_list(tlg());
      
      // simple aspiration criteria
      mets::best_ever_criteria aspiration_criteria;
      
      mets::threshold_termination_criteria 
	threshold(0);
      
      // combine threshold with a max noimprove criterion
      mets::noimprove_termination_criteria 
	minor_it_criteria(&threshold, 20);
      
      while(!minor_it_criteria(major_best.best_seen()))
	{
	  vcp minor_store(n, edges);
	  mets::best_ever_solution minor_best(minor_store);
	  // combine threshold with a max noimprove criterion
	  mets::noimprove_termination_criteria noimprove(&threshold, 5000);
	  // random tabu list tenure
	  tabu_list.tenure(tlg());
	  // Create a tabu_search algorithm instance starting from "model",
	  // recording the best solution in "best", exploring the neighborhood
	  // using "neigh", using the tabu list "tabu_list", the best ever
	  // aspiration criteria "aspiration_criteria" and the combined
	  // termination criteria "threshold_noimprove".
	  mets::tabu_search<vcp_neighborhood> algorithm(point, 
							minor_best, 
							neigh, 
							tabu_list, 
							aspiration_criteria, 
							noimprove);
	  algorithm.attach(log);
	  std::clog << "New iteration with tenure: " 
		    << tabu_list.tenure();
	  algorithm.search();
	  std::clog << " -> " << minor_best.best_cost() << std::endl;
	  major_best.accept(minor_best.best_seen());
	  point.copy_from(major_best.best_seen());
	  point.perturbate(g_colors, n/4, gen);
	}
      
       best_recorder.accept(major_best.best_seen());
      clog << "Best of this run/so far: " 
	   << major_best.best_cost()  
	   << "/"  << best_recorder.best_cost() << endl;
      
      if(best_recorder.best_cost() == 0) break;
    }

  vcp_neighborhood neigh2(g_colors);
  mets::local_search<vcp_neighborhood> algo2(point, 
					     best_recorder, 
					     neigh2, 
					     false);
  
  algo2.attach(log);
  if(!fast) algo2.search();

  flog.close();

  clog << "Best solution: " << best_recorder.best_cost()  << endl;

  vcp& incumbent = ((vcp&)best_recorder.best_seen());

  ofstream sol((string(argv[1]) + ".log").c_str());
  incumbent.print(sol);
  sol.close();

  ofstream dot((string(argv[1]) + ".dot").c_str());
  incumbent.print_dot(g_colors, dot);
  dot.close();

  incumbent.display_conflicts(clog);

  return 0;
}
