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
#include <metslib/mets.h>

/// @brief The tutorial model is a simple model for the subset sum problem
class subsetsum : public mets::feasible_solution {
  std::vector<bool> delta_m;
  std::vector<int> set_m;
  int sum_m;
public:
  subsetsum(const std::vector<int>& set, int sum) 
    : delta_m(set.size()), 
      set_m(set.begin(), set.end()),
      sum_m(sum)
  { }

  /// @brief Pure virtual cost_function to be minimized
  ///
  /// min set_m' delta_m
  /// s.t. set_m' delta_m <= sum_m
  ///
  mets::gol_type cost_function() const
  {
    int diff = sum_m - std::inner_product
      (set_m.begin(), set_m.end(), delta_m.begin(), 0);
    
    // allow, but hardly penalizes a constraint violation
    if(diff < 0)
      return -100 * diff;
    else
      return diff;
  }

   /// @brief This method is needed by the algorithm to record the
  /// best solution.
  void copy_from(const mets::feasible_solution& o)
  {
    const subsetsum& s = dynamic_cast<const subsetsum&>(o);
    delta_m = s.delta_m;
    set_m = s.set_m;
  }

  /// @brief Utility methods
  size_t size() const
  { return delta_m.size(); }

  /// @brief Read/write operators
  bool delta(int i) const { return delta_m[i]; }
  void delta(int i, bool val) { delta_m[i] = val; }
  int element(int i) const { return set_m[i]; }
};

std::ostream& operator<<(std::ostream& o, const subsetsum& s)
{
  for(int ii(0); ii!=s.size(); ++ii)
    if(s.delta(ii))
      o << s.element(ii) << ' ';
  return o;
}
