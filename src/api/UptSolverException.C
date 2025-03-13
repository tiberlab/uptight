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
#include "UptSolverException.h"

const char* ETBSolverException::what(void) const throw()
{
   switch(_error)
   {
   case _ERR_ALLOC_ERR:
     return "allocation error"; 
	    
   case _ERR_ATM_MISMTCH:
     return "number of atom mismatch"; 	    

   case _ERR_EGV_VOID:
     return "empty eigenvectors"; 	    

   case _ERR_EGV_MISMTCH: 
     return "eigenvector mismatch"; 	    

   case _ERR_MAT_UNALLOC: 
     return "matrix not allocated"; 	    

   case _ERR_MAT_NOTSQRE: 
     return "matrix not square"; 	    

   case _ERR_NN_LIST:  
     return "nearest neighbour list error"; 	    

   case _ERR_HAM_UNMAT:   
     return "unrecognized material"; 	    

   case _ERR_HAM_UNPAIR:  
     return "unrecognized pair";

   case _ERR_HAM_ZEROLN:
     return "unsupported matrix with whole line of zero";  

   case _ERR_REF_ATMLIST:
     return "uncorrect reference atom list: check code";

   case _ERR_REF_NAME:
     return "symbol not found in reference table";  

   case _ERR_ALLOY_TBBLK:
     return "TB block mismatch in alloy";

   case _ERR_DB_PAIR:
     return "wrong pair - check database";
     
   case _ERR_INPUT:
     return "uncorrect input or database";

   case _ERR_DB_NOINT:
     return "integer not found";

   case _ERR_DB_NOBOOL:
     return "logical flag not found";

   case _ERR_FILE_OPEN:
     return "file not found";

   case _ERR_OUTPUT:
     return "output error";

   case _ERR_LANCZ_DIAG:
     return "error in lanczos solver";

   case _ERR_FEAST_DIAG:
     return "error in feast solver";

   case _ERR_LAPACK_DIAG: 
     return "error in lapack solver";

   default:
     return "general internal error";
   }
}


extern "C" {
  void throw_solve_exception_(const int& errcode)
  {
     throw ETBSolverException(errcode);
  }
}
