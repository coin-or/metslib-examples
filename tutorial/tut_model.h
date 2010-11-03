#pragma once
// METSlib tutorial source file - tut_model.cc                 -*- C++ -*-
//
// Copyright (C) 2009 Mirko Maischberger <mirko.maischberger@gmail.com>
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <vector>
#include <numeric>
#include <cassert>
#include <algorithm>
#include <metslib/mets.hh>

/// @brief The tutorial model is a simple model for the subset sum problem
class subsetsum : public mets::copyable_solution {
  /// @brief The binary variables 
  std::vector<bool> delta_m;
  /// @brief The values (parameters) of the problem
  std::vector<int> set_m;
  /// @brief The target sum
  int target_sum_m;
  /// @brief The actual cost
  int current_sum_m;
public:
  /// @brief Ctor.
  subsetsum(const std::vector<int>& set, int sum) 
    : delta_m(set.size(), false), 
      set_m(set.begin(), set.end()),
      target_sum_m(sum),
      current_sum_m(0)
  { }

  /// @brief The cost_function that we want minimized
  ///
  /// min set_m' delta_m
  /// s.t. set_m' delta_m <= target_sum_m
  ///
  mets::gol_type cost_function() const
  {
    // the method allows, but hardly penalizes a constraint violation
    int diff = target_sum_m - current_sum_m;
    if(diff < 0)
      return -100 * diff;
    else
      return diff;
  }

  /// @brief This method is needed by the algorithm to record the best
  /// solution.
  void copy_from(const mets::copyable& o)
  {
    const subsetsum& s = dynamic_cast<const subsetsum&>(o);
    delta_m = s.delta_m;
    set_m = s.set_m;
    target_sum_m = s.target_sum_m;
    current_sum_m = s.current_sum_m;
  }

  /// @brief The size of the problem
  size_t size() const
  { return delta_m.size(); }

  /// @brief Evaluates the cost of a change without actually doing it.
  mets::gol_type what_if(int i, bool val) const
  {
    int newcost = current_sum_m;
    if(delta_m[i] && !val)
      newcost -= set_m[i];
    else if(!delta_m[i] && val)
      newcost += set_m[i];
    // the method allows, but hardly penalizes a constraint violation
    int diff = target_sum_m - newcost;
    if(diff < 0)
      return -100 * diff;
    else
      return diff;
  }
  
  /// @brief Return actual delta[i] value
  bool delta(int i) const { return delta_m[i]; }

  /// @brief Set delta[i] to val and update the cost
  void delta(int i, bool val) 
  {
    if(delta_m[i] && !val)
      current_sum_m -= set_m[i];
    else if(!delta_m[i] && val)
      current_sum_m += set_m[i];
    delta_m[i] = val;
  }

  int element(int i) const { return set_m[i]; }

};

/// @brief Print the solution on a stream.
std::ostream& operator<<(std::ostream& o, const subsetsum& s)
{
  for(int ii(0); ii!=s.size(); ++ii)
    if(s.delta(ii))
      o << s.element(ii) << ' ';
  return o;
}
