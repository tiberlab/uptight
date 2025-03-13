! This file is part of uptight.
!
! uptight is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! uptight is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with uptight. If not, see <https://www.gnu.org/licenses/>.
!
SUBROUTINE throw_init_exception(CODE)
  USE exceptions
  implicit none     
  integer(4) :: CODE

  SELECT CASE(CODE)
  
   case (ERR_ALLOC_ERR)
     write(*,*) "allocation error";
   case (ERR_ATM_MISMTCH)
     write(*,*) "number of atom mismatch";  
   case (ERR_EGV_VOID)
     write(*,*) "empty eigenvectors";     
   case (ERR_EGV_MISMTCH) 
     write(*,*) "eigenvector mismatch";    
   case (ERR_MAT_UNALLOC) 
     write(*,*) "matrix not allocated";    
   case (ERR_MAT_NOTSQRE) 
     write(*,*) "matrix not square";     
   case (ERR_NN_LIST)  
     write(*,*) "nearest neighbour list error";   
   case (ERR_SCHE_INCON) 
     write(*,*) "ETB schemes of .etb files not the same";     
   case (ERR_HAM_UNSCHE) 
     write(*,*) "unrecognized ETB scheme";           
   case (ERR_HAM_UNMAT)   
     write(*,*) "unrecognized material";     
   case (ERR_HAM_UNPAIR)  
     write(*,*) "unrecognized pair";
   case (ERR_HAM_ZEROLN)
     write(*,*) "unsupported matrix with whole line of zero";  
   case (ERR_REF_ATMLIST)
     write(*,*) "uncorrect reference atom list check code";
   case (ERR_REF_NAME)
     write(*,*) "symbol not found in reference table";  
   case (ERR_ALLOY_TBBLK)
     write(*,*) "TB block mismatch in alloy";
   case (ERR_DB_PAIR)
     write(*,*) "wrong pair - check database";
   case (ERR_INPUT)
     write(*,*) "incorrect input or database";
   case (ERR_DB_NOINT)
     write(*,*) "integer not found";
   case (ERR_DB_NOBOOL)
     write(*,*) "logical flag not found";
   case (ERR_FILE_OPEN)
     write(*,*) "file not found";
   case (ERR_OUTPUT)
     write(*,*) "output error";
   case default
     write(*,*) "general internal error";
   end select
  
   STOP

END SUBROUTINE       
 
SUBROUTINE throw_solve_exception(CODE)
  USE exceptions
  implicit none     
  integer(4) :: CODE

  SELECT CASE(CODE)
  
   case (ERR_ALLOC_ERR)
     write(*,*) "allocation error";
   case (ERR_ATM_MISMTCH)
     write(*,*) "number of atom mismatch";   
   case (ERR_EGV_VOID)
     write(*,*) "empty eigenvectors"; 
   case (ERR_EGV_MISMTCH) 
     write(*,*) "eigenvector mismatch";    
   case (ERR_MAT_UNALLOC) 
     write(*,*) "matrix not allocated";    
   case (ERR_MAT_NOTSQRE) 
     write(*,*) "matrix not square";
   case (ERR_NN_LIST)  
     write(*,*) "nearest neighbour list error";    
   case (ERR_BOND_NUM)
     write(*,*) "more 1st and/or 2nd neighbors than expected";
   case (ERR_HAM_UNMAT)   
     write(*,*) "unrecognized material";    
   case (ERR_HAM_UNPAIR)  
     write(*,*) "unrecognized pair";
   case (ERR_HAM_ZEROLN)
     write(*,*) "unsupported matrix with whole line of zero";  
   case (ERR_REF_ATMLIST)
     write(*,*) "uncorrect reference atom list check code";
   case (ERR_REF_NAME)
     write(*,*) "symbol not found in reference table";  
   case (ERR_ALLOY_TBBLK)
     write(*,*) "TB block mismatch in alloy";
   case (ERR_DB_PAIR)
     write(*,*) "wrong pair - check database";
   case (ERR_INPUT)
     write(*,*) "uncorrect input or database";
   case (ERR_DB_NOINT)
     write(*,*) "integer not found";
   case (ERR_DB_NOBOOL)
     write(*,*) "logical flag not found";
   case (ERR_FILE_OPEN)
     write(*,*) "file not found";
   case (ERR_OUTPUT)
     write(*,*) "output error";
   case default
     write(*,*) "general internal error";
   end select
  
   STOP

END SUBROUTINE        
