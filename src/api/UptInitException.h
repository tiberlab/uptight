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
#ifndef _UPTINITEXCEPTION_H_
#define _UPTINITEXCEPTION_H_

#include "InitFailedException.h"
#include "exception_codes.h"

//! An exception class for the solver interfaces
class ETBInitException : public InitFailedException 
{

  public:
     ETBInitException(const int errorcode)
	     : InitFailedException(""), _error(errorcode){};

     virtual const char* what() const throw();  

  private:
     int _error;
};


#endif // _UPTINITEXCEPTION_H_
