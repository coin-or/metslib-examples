#pragma once
#include <typeinfo>
#include <metslib/mets.h>

#include "qap_model.hpp"

template<typename T>
class qap_neighborhood;

/// @brief A mets::mana_move that swaps two object
///
/// Each instance swaps two specific objects.
///
class qap_move : public mets::mana_move 
{
public:  
  template<typename T>
  friend class qap_neighborhood;

  /// @brief A move that swaps from and to.
  qap_move(int from, int to) 
    : p1(std::min(from,to)), p2(std::max(from,to)) 
  { }

  /// @brief Virtual method that applies the move on a point
  void
  apply(mets::feasible_solution& s)
  {
    qap_model& sol = reinterpret_cast<qap_model&>(s);
    sol.swap(p1, p2);
  }

  /// @brief Unapply the last move: in case of a swap the inverse move
  /// is just the same swap.
  void
  unapply(mets::feasible_solution& s)
  {
    this->apply(s);
  }

  /// @brief A method to clone self. Needed to insert the move in a
  /// tabu list.
  mana_move* 
  clone() const 
  { return new qap_move(p1, p2); }

  /// @brief An hash function used by the tabu list (the hash value is
  /// used to insert the move in an hash set).
  size_t
  hash() const
  { return (p1)<<16^(p2); }

  /// @brief Comparison operator used to tell if this move is equal to
  /// a move in the tabu list.
  bool 
  operator==(const mets::mana_move& o) const
  {
    try {
      const qap_move& other = dynamic_cast<const qap_move&>(o);
      return (this->p1 == other.p1 && this->p2 == other.p2);
    } catch (std::bad_cast& e) {
      return false;
    }
  }
  
protected:
  int p1;
  int p2;
};

