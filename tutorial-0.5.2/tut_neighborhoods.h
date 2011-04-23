#pragma once
// METSlib tutorial source file - tut_neighborhoods.hh         -*- C++ -*-
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

#include <iostream>
#include "tut_moves.h"

class full_neighborhood  
{
public:
  std::vector<toggle*> moves_m;
  typedef std::vector<toggle*>::iterator iterator;
  iterator begin() { return moves_m.begin(); }
  iterator end() { return moves_m.end(); }
  
  full_neighborhood(int problem_size) 
  {
    for(unsigned int ii = 0; ii != problem_size; ++ii) 
      moves_m.push_back(new toggle(ii));
  }
  
  ~full_neighborhood()
  {
    // delete all moves
    for(iterator ii = begin(); ii != end(); ++ii)
      delete (*ii);
  }

  void refresh(mets::feasible_solution& s) 
  { 
    // This method can be used to adapt the neighborhood to the
    // current problem instance.  In our simple case there is no need
    // to update the neighborhood since its moves does not depend on
    // the current solution considered.
  }

};
