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
MODULE exceptions
#include "exception_codes.h"
   implicit none
     

INTEGER(4),PARAMETER:: ERR_ALLOC_ERR  =  _ERR_ALLOC_ERR 
                                                        
INTEGER(4),PARAMETER:: ERR_ATM_MISMTCH = _ERR_ATM_MISMTCH
INTEGER(4),PARAMETER:: ERR_EGV_VOID    = _ERR_EGV_VOID  
INTEGER(4),PARAMETER:: ERR_EGV_MISMTCH = _ERR_EGV_MISMTCH
INTEGER(4),PARAMETER:: ERR_MAT_UNALLOC = _ERR_MAT_UNALLOC
INTEGER(4),PARAMETER:: ERR_MAT_NOTSQRE = _ERR_MAT_NOTSQRE
INTEGER(4),PARAMETER:: ERR_MAT_NOTINV  = _ERR_MAT_NOTINV 
                                  
INTEGER(4),PARAMETER:: ERR_NN_LIST     = _ERR_NN_LIST
INTEGER(4),PARAMETER:: ERR_BOND_NUM    = _ERR_BOND_NUM
INTEGER(4),PARAMETER:: ERR_NN_NOTSYM   = _ERR_NN_NOTSYM
INTEGER(4),PARAMETER:: ERR_NN_DGBOND   = _ERR_NN_DGBOND
                                                        
INTEGER(4),PARAMETER:: ERR_SCHE_INCON  = _ERR_SCHE_INCON
INTEGER(4),PARAMETER:: ERR_HAM_UNSCHE  = _ERR_HAM_UNSCHE
INTEGER(4),PARAMETER:: ERR_HAM_UNMAT   = _ERR_HAM_UNMAT  
INTEGER(4),PARAMETER:: ERR_HAM_UNPAIR  = _ERR_HAM_UNPAIR 
INTEGER(4),PARAMETER:: ERR_HAM_ZEROLN  = _ERR_HAM_ZEROLN 

INTEGER(4),PARAMETER:: ERR_REF_ATMLIST = _ERR_REF_ATMLIST
INTEGER(4),PARAMETER:: ERR_REF_NAME    = _ERR_REF_NAME   

INTEGER(4),PARAMETER:: ERR_ALLOY_TBBLK = _ERR_ALLOY_TBBLK

INTEGER(4),PARAMETER:: ERR_INPUT       = _ERR_INPUT
INTEGER(4),PARAMETER:: ERR_OUTPUT      = _ERR_OUTPUT
INTEGER(4),PARAMETER:: ERR_FILE_OPEN   = _ERR_FILE_OPEN
INTEGER(4),PARAMETER:: ERR_DB_NOINT    = _ERR_DB_NOINT 
INTEGER(4),PARAMETER:: ERR_DB_NOBOOL   = _ERR_DB_NOBOOL
INTEGER(4),PARAMETER:: ERR_DB_ATOM     = _ERR_DB_ATOM   
INTEGER(4),PARAMETER:: ERR_DB_PAIR     = _ERR_DB_PAIR  

INTEGER(4),PARAMETER:: ERR_LANCZ_DIAG  = _ERR_LANCZ_DIAG 
INTEGER(4),PARAMETER:: ERR_JD_DIAG     = _ERR_JD_DIAG 
INTEGER(4),PARAMETER:: ERR_LAPACK_DIAG = _ERR_LAPACK_DIAG 
INTEGER(4),PARAMETER:: ERR_FEAST_DIAG  = _ERR_FEAST_DIAG 

INTEGER(4),PARAMETER:: ERR_GENERAL     = _ERR_GENERAL
INTEGER(4),PARAMETER:: ERR_IO_ERROR    = _ERR_IO_ERROR

END MODULE exceptions
