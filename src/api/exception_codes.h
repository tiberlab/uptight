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
#ifndef _EXCEPTION_CODES_H
#define _EXCEPTION_CODES_H

#define _ERR_ALLOC_ERR   100

#define _ERR_ATM_MISMTCH 201
#define _ERR_EGV_VOID    202
#define _ERR_EGV_MISMTCH 203
#define _ERR_MAT_UNALLOC 204
#define _ERR_MAT_NOTSQRE 205
#define _ERR_MAT_NOTINV  206

#define _ERR_NN_LIST     250
#define _ERR_BOND_NUM    251
#define _ERR_NN_NOTSYM   260
#define _ERR_NN_DGBOND   261

#define _ERR_HAM_UNMAT   301
#define _ERR_HAM_UNPAIR  304
#define _ERR_HAM_ZEROLN  305
#define _ERR_HAM_UNSCHE  306 
#define _ERR_SCHE_INCON  316 
#define _ERR_ALLOY_TBBLK 326

#define _ERR_REF_ATMLIST 401
#define _ERR_REF_NAME    402

#define _ERR_FILE_OPEN   500
#define _ERR_INPUT       501
#define _ERR_OUTPUT      502
#define _ERR_DB_NOINT    503
#define _ERR_DB_NOBOOL   504
#define _ERR_IO_ERROR    510

#define _ERR_DB_ATOM     601
#define _ERR_DB_PAIR     602

#define _ERR_LANCZ_DIAG  701
#define _ERR_LAPACK_DIAG 721
#define _ERR_FEAST_DIAG  731
#define _ERR_JD_DIAG     741

#define _ERR_GENERAL     666



#endif

