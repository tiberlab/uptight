/*
 * This file is part of uptight.
 *
 * uptight is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * uptight is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with uptight. If not, see <https://www.gnu.org/licenses/>.
 */
// $Id: SolveFailedException.h 126 2006-11-24 08:10:44Z maufder $


#ifndef _SOLVEFAILEDEXCEPTION_H_
#define _SOLVEFAILEDEXCEPTION_H_

#include <stdexcept>
#include <string>

//! An exception class for failed solve
class SolveFailedException : public std::runtime_error
{

  public:
    SolveFailedException(const char* msg)
      : std::runtime_error(msg) {};

    SolveFailedException(const std::string& msg)
      : std::runtime_error(msg) {};


  private:

};



#endif // _SOLVEFAILEDEXCEPTION_H_
