/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! FILE created  by Huiying Xu 2007
!!      modified by Eric Polizzi 2009 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/


#include "f90_noruntime.h"



double wdsin_(double *x){
return sin(*x);
}

double wdcos_(double *x){
return cos(*x);
}


float wssin_(float *x){
return sin(*x);
}

float wscos_(float *x){
return cos(*x);
}



/*void wstop_ () {

	exit(0);

}*/


void wdeallocate_all_type(array_all_dim *p) {

/*
 * This function deallocates the space of an array.
 * It works for all dimensions of arrays.
 */ 

        long int i, j;     
/*
        if (p->base != (void *) 0) free(p->base);   
*/
        free(p->base);   
        p->base = (void *) 0;   

	return;
   
        j = p->ndim;  
  
        p->size = 0;   
        p->offset = (void *) 0;   
        p->bitsets = 0;   
        p->ndim = 0;   
        p->reserved = 0;   
  
        for (i = 0; i < j; i++)  
                memset(p->dim_info+i, 0, sizeof(dimension_info));   
}


void wdeallocate_1i_(array_dim1 *p) {  

/*
 * This function deallocates the space of one-dimensional integer array.
 */

	array_all_dim *p1 = (array_all_dim *) p; 

	wdeallocate_all_type(p1);

} 
 
 
void wdeallocate_2i_(array_dim2 *p) {  

/*
 * This function deallocates the space of two-dimensional integer array.
 */

	array_all_dim *p1 = (array_all_dim *) p;
 
	wdeallocate_all_type(p1);    

}  


void wdeallocate_3i_(array_dim3 *p) { 

/*
 * This function deallocates the space of three-dimensional integer array.
 */

	array_all_dim *p1 = (array_all_dim *) p;  
  
        wdeallocate_all_type(p1);     
 
}   


void wdeallocate_1d_(array_dim1 *p) {

/* 
 * This function deallocates the space of one-dimensional double array.
 */

	array_all_dim *p1 = (array_all_dim *) p;
 
	wdeallocate_all_type(p1);
	
}

void wdeallocate_2d_(array_dim2 *p) {    
   
/*
 * This function deallocates the space of two-dimensional double array.
 */
        
	array_all_dim *p1 = (array_all_dim *) p;

        wdeallocate_all_type(p1);   
}   

void wdeallocate_3d_(array_dim3 *p) {

/*
 * This function deallocates the space of three-dimensional double array.
 */

	array_all_dim *p1 = (array_all_dim *) p;     
            
        wdeallocate_all_type(p1);    
}    



void wdeallocate_1z_(array_dim1 *p) {

/* 
 * This function deallocates the space of one-dimensional double complex array.
 */

	array_all_dim *p1 = (array_all_dim *) p;
 
	wdeallocate_all_type(p1);
	
}

void wdeallocate_2z_(array_dim2 *p) {    
   
/*
 * This function deallocates the space of two-dimensional double complex array.
 */
        
	array_all_dim *p1 = (array_all_dim *) p;

        wdeallocate_all_type(p1);   
}   

void wdeallocate_3z_(array_dim3 *p) {

/*
 * This function deallocates the space of three-dimensional double complex array.
 */

	array_all_dim *p1 = (array_all_dim *) p;     
            
        wdeallocate_all_type(p1);    
}    



void wdeallocate_1s_(array_dim1 *p) {

/* 
 * This function deallocates the space of one-dimensional real array.
 */

	array_all_dim *p1 = (array_all_dim *) p;
 
	wdeallocate_all_type(p1);
	
}

void wdeallocate_2s_(array_dim2 *p) {    
   
/*
 * This function deallocates the space of two-dimensional real array.
 */
        
	array_all_dim *p1 = (array_all_dim *) p;

        wdeallocate_all_type(p1);   
}   

void wdeallocate_3s_(array_dim3 *p) {

/*
 * This function deallocates the space of three-dimensional real array.
 */

	array_all_dim *p1 = (array_all_dim *) p;     
            
        wdeallocate_all_type(p1);    
}    



void wdeallocate_1c_(array_dim1 *p) {

/* 
 * This function deallocates the space of one-dimensional complex array.
 */

	array_all_dim *p1 = (array_all_dim *) p;
 
	wdeallocate_all_type(p1);
	
}

void wdeallocate_2c_(array_dim2 *p) {    
   
/*
 * This function deallocates the space of two-dimensional complex array.
 */
        
	array_all_dim *p1 = (array_all_dim *) p;

        wdeallocate_all_type(p1);   
}   

void wdeallocate_3c_(array_dim3 *p) {

/*
 * This function deallocates the space of three-dimensional complex array.
 */

	array_all_dim *p1 = (array_all_dim *) p;     
            
        wdeallocate_all_type(p1);    
}    





void wallocate_all_type(array_all_dim *p, int * firstNo, int *secondNo, int *thirdNo, int type, int *alloc_info) {

/*
 * This function allocates the space for an array.
 * It works for all dimensions of arrays.
 */

        dimension_info dim;   
        void *base = NULL;   
        /* int element_size; */
	size_t element_size;

	switch(type) {
	    case F90_INTEGER:
		element_size = sizeof(int);
		break;
	    case F90_REAL:
		element_size = sizeof(float);
		break;
            case F90_COMPLEX:
		element_size = 2*sizeof(float);
		break;
	    case F90_DOUBLE_COMPLEX:
		element_size = 2*sizeof(double);
		break;		
	    default:   /*case F90_DOUBLE_PRECISION: */
		element_size = sizeof(double);
	}


	if ((*thirdNo) != 0) {
		base = (void *) malloc(element_size * (*firstNo) * (*secondNo) * (*thirdNo));
		if (base == NULL) {
			*alloc_info = ALLOC_FAIL;
			return;
		}
		else
			*alloc_info = ALLOC_SUCCESS;
		p->ndim = 3;
	}
	else if ((*secondNo) != 0) {
		base = (void *) malloc(element_size * (*firstNo) * (*secondNo));
                if (base == NULL) {
                        *alloc_info = ALLOC_FAIL;
			return;
		}
                else
                        *alloc_info = ALLOC_SUCCESS; 
                p->ndim = 2; 
        }
	else { 
                base = (void *) malloc((*firstNo) * element_size);
                if (base == NULL) {
                        *alloc_info = ALLOC_FAIL;
			return;
		}
                else
                        *alloc_info = ALLOC_SUCCESS; 
                p->ndim = 1; 
	}


        p->base = base;   
        p->size = (long int) element_size;    
        p->offset = 0;   
        p->bitsets = 5;   
        p->reserved = 0;  

        dim.nelts = (long int)(*firstNo);      
        dim.stride = (long int) element_size;      
        dim.lower_bound = 1;     
        memcpy(p->dim_info, &dim, sizeof(dimension_info));  	

	if ((*secondNo) != 0) {
		dim.nelts = (long int)(*secondNo);       
        	dim.stride = (long int)((*firstNo) * element_size);       
        	dim.lower_bound = 1;    
		memcpy(p->dim_info+1, &dim, sizeof(dimension_info));
	}

        if ((*thirdNo) != 0) { 
                dim.nelts = (long int) (*thirdNo);        
                dim.stride = (long int) ((*firstNo) * (*secondNo) * element_size);        
                dim.lower_bound = 1;     
                memcpy(p->dim_info+2, &dim, sizeof(dimension_info)); 
        } 
}




void wallocate_1i_(array_dim1 *p, int * element_No, int *alloc_info) {   

/* 
 * This function allocates the space for one-dimensional integer array.
 */

	int i = 0, j = 0;
	array_all_dim *p1 = (array_all_dim *)p;

	wallocate_all_type(p1, element_No, &i, &j, F90_INTEGER, alloc_info);

}  


 
void wallocate_2i_(array_dim2 *p, int * rowNo, int *columNo, int *alloc_info) {  
	
/*
 * This function allocates the space for two-dimensional integer array.
 */

	int j = 0;
	array_all_dim *p1 = (array_all_dim *) p;
 
	wallocate_all_type(p1, rowNo, columNo, &j, F90_INTEGER, alloc_info);

} 


void wallocate_3i_(array_dim3 *p, int * firstNo, int *secondNo, int *thirdNo, int *alloc_info) {   

/* 
 * This function allocates the space for three-dimensional integer array.
 */   
      
	array_all_dim *p1 = (array_all_dim *) p;
        wallocate_all_type(p1, firstNo, secondNo, thirdNo, F90_INTEGER, alloc_info); 
}  




void wallocate_1d_(array_dim1 *p, int * element_No, int *alloc_info) {  

/*
 * This function allocates the space for one-dimensional double array.
 */

	int i = 0, j = 0;
	array_all_dim *p1 = (array_all_dim *) p;
 
	wallocate_all_type(p1, element_No, &i, &j, F90_DOUBLE_PRECISION, alloc_info);

} 



void wallocate_2d_(array_dim2 *p, int * rowNo, int *columNo, int *alloc_info) {

/*
 * This function allocates the space for two-dimensional double array.
 */

	int j = 0; 
	array_all_dim *p1 = (array_all_dim *) p;
	
	wallocate_all_type(p1, rowNo, columNo, &j, F90_DOUBLE_PRECISION, alloc_info);

}


void wallocate_3d_(array_dim3 *p, int * firstNo, int *secondNo, int *thirdNo, int *alloc_info) { 

/*
 * This function allocates the space for three-dimensional double array.
 */ 

	array_all_dim *p1 = (array_all_dim *) p;
        wallocate_all_type(p1, firstNo, secondNo, thirdNo, F90_DOUBLE_PRECISION, alloc_info); 
 
} 


void wallocate_1z_(array_dim1 *p, int * element_No, int *alloc_info) {  

/*
 * This function allocates the space for one-dimensional double array.
 */

	int i = 0, j = 0;
	array_all_dim *p1 = (array_all_dim *) p;
 
	wallocate_all_type(p1, element_No, &i, &j, F90_DOUBLE_COMPLEX, alloc_info);

} 



void wallocate_2z_(array_dim2 *p, int * rowNo, int *columNo, int *alloc_info) {

/*
 * This function allocates the space for two-dimensional double array.
 */

	int j = 0; 
	array_all_dim *p1 = (array_all_dim *) p;
	
	wallocate_all_type(p1, rowNo, columNo, &j, F90_DOUBLE_COMPLEX, alloc_info);

}


void wallocate_3z_(array_dim3 *p, int * firstNo, int *secondNo, int *thirdNo, int *alloc_info) { 

/*
 * This function allocates the space for three-dimensional double array.
 */ 

	array_all_dim *p1 = (array_all_dim *) p;
        wallocate_all_type(p1, firstNo, secondNo, thirdNo, F90_DOUBLE_COMPLEX, alloc_info); 
 
} 




void wallocate_1s_(array_dim1 *p, int * element_No, int *alloc_info) {  

/*
 * This function allocates the space for one-dimensional real array.
 */

	int i = 0, j = 0;
	array_all_dim *p1 = (array_all_dim *) p;
 
	wallocate_all_type(p1, element_No, &i, &j, F90_REAL, alloc_info);

} 



void wallocate_2s_(array_dim2 *p, int * rowNo, int *columNo, int *alloc_info) {

/*
 * This function allocates the space for two-dimensional real array.
 */

	int j = 0; 
	array_all_dim *p1 = (array_all_dim *) p;
	
	wallocate_all_type(p1, rowNo, columNo, &j, F90_REAL, alloc_info);

}


void wallocate_3s_(array_dim3 *p, int * firstNo, int *secondNo, int *thirdNo, int *alloc_info) { 

/*
 * This function allocates the space for three-dimensional real array.
 */ 

	array_all_dim *p1 = (array_all_dim *) p;
        wallocate_all_type(p1, firstNo, secondNo, thirdNo, F90_REAL, alloc_info); 
 
} 


void wallocate_1c_(array_dim1 *p, int * element_No, int *alloc_info) {  

/*
 * This function allocates the space for one-dimensional complex array.
 */

	int i = 0, j = 0;
	array_all_dim *p1 = (array_all_dim *) p;
 
	wallocate_all_type(p1, element_No, &i, &j, F90_COMPLEX, alloc_info);

} 



void wallocate_2c_(array_dim2 *p, int * rowNo, int *columNo, int *alloc_info) {

/*
 * This function allocates the space for two-dimensional complex array.
 */

	int j = 0; 
	array_all_dim *p1 = (array_all_dim *) p;
	
	wallocate_all_type(p1, rowNo, columNo, &j, F90_COMPLEX, alloc_info);

}


void wallocate_3c_(array_dim3 *p, int * firstNo, int *secondNo, int *thirdNo, int *alloc_info) { 

/*
 * This function allocates the space for three-dimensional complex array.
 */ 

	array_all_dim *p1 = (array_all_dim *) p;
        wallocate_all_type(p1, firstNo, secondNo, thirdNo, F90_COMPLEX, alloc_info); 
 
} 



/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/



void wopen_(unsigned long int *file_id, char *filename, char *status) {

/* In Fortran, once a file exists, we cannot open it as 'new'. 
 * In C, we can choose the appropriate model so that we don't care 
 * if the file already exists.
 */  

	FILE *fp = fopen(filename, status); 
	*file_id = (unsigned long int) fp;
}


void wread_line_(unsigned long int *file_id, void *buffer, int *type) { 
 
/* 
 * This function writes the char (or string), integer, real, and double precision into the file.
 * When *type = 0 (F90_CHAR), the function reads the first character in the line 
 * When *type = 1 (F90_STRING), the function reads the whole line 
 * When *type = 2 F90_INTEGER, the function reads the first integer in the line 
 * When *type = -1 (F90_LOGICAL), the function read the first integer (logical) in the line
 * When *type = 3 (F90_REAL), the function reads the first float in the line 
 * WHen *type = 4 (F90_DOUBLE_PRECISION), the function reads the first double in the line 
 */ 

	int i;
	char *temp = NULL;
	size_t len = 0;
	getline(&temp, &len, (FILE *) *file_id);

        switch (*type) { 
	    case F90_LOGICAL:
		if (temp[0] == 'T' || temp[0] == 't' || strncmp(temp, "-1", 2) == 0 || strncmp(temp, ".true.", 6) == 0)
			* (int *) buffer = -1;
		else
			* (int *) buffer = 0;
	    case F90_CHAR:
		* (char *) buffer = temp[0];
		break;
            case F90_STRING: 
		for (i = 0; i < 500; i++)  /* for (i = 0; i < len; i++) */
		{
			if (temp[i] != '!' && temp[i] != '\0' && temp[i] != ' ' && temp[i] != '\t') {
				* (char *) ((unsigned long int)buffer + i*sizeof(char)) = temp[i];
			}
			else {
				* (char *) ((unsigned long int)buffer + i*sizeof(char)) = '\0';
				break;
			}
		}
                break; 
            case F90_INTEGER:
		* (int *) buffer = atoi(temp); 
                break; 
            case F90_REAL:
		* (float *) buffer = atof(temp); 
                break; 
            case F90_DOUBLE_PRECISION: 
		* (double *) buffer = atof(temp);
		break;

	    /* We might add more choices here, depending on our need */

	    case F90_INTEGER_STAR_8:
		* (unsigned long int *) buffer = strtoul(temp, NULL, 16); 
		break;
	    default:
		return;  
        } 
} 



void wwrite_(unsigned long int *file_id, void *buffer, int *type) {

/*
 * This function writes the char (or string), integer, real, and double precision into the file.
 * When *type = -2 (F90_WRITE_NOTATION), the function writes '\n', '\t', ... into the file
 * When *type = 1 (F90_STRING), the function writes a character or string into the file
 * When *type = 2 (F90_INTEGER), the function writes an integer into the file
 * When *type = 3 (F90_REAL), the function writes a real into the file
 * WHen *type = 4 (F90_DOUBLE_PRECISION), the function writes double precision into the file 
 */

	int flag = 0;

        if (*file_id == (unsigned long int) 6) flag = -1;

	switch (*type) {
	    case F90_LOGICAL:
		if (*(int *)buffer == -1) {
			if (flag == 0) fprintf((FILE *) *file_id, "T");
			else printf("T");
		}
		else {
			if (flag == 0) fprintf((FILE *) *file_id, "F"); 
                        else printf("F");
		}
		break;	
	    case F90_WRITE_NOTATION:

		if (strncmp((char *) buffer, "\\n", 2) == 0) {
			if (flag == 0) fprintf((FILE *) *file_id, "\n");
			else printf("\n");
		}
		else {  /* \\t */
			if (flag == 0) fprintf((FILE *) *file_id, "\t");
			else printf("\t");
		}


		/* we might add more choicse here, depending on our need */

		break;
	    case F90_CHAR:
		if (flag == 0) fprintf((FILE *) *file_id, "%c", * (char *) buffer);
		else printf("%c", * (char *) buffer);
		break;
	    case F90_STRING:	/*    case F90_WRITE_NOTATION:  */
		if (flag == 0) fprintf((FILE *) *file_id, "%s", (char *) buffer);
		else printf("%s", (char *) buffer);
		break;
	    case F90_INTEGER:
		if (flag == 0) fprintf((FILE *) *file_id, "%d", * (int *)buffer);
		else printf("%d", * (int *)buffer);
		break;
	    case F90_REAL:
		if (flag == 0) fprintf((FILE *) *file_id, "%f", * (float *)buffer);
		else printf("%f", * (float *)buffer);
		break;
	    case F90_DOUBLE_PRECISION:
		if (flag == 0) fprintf((FILE *) *file_id, "%e", * (double *)buffer);
		else printf("%.15e", * (double *)buffer);
		/*
		printf("The double precision being written is %e\n", * (double *) buffer);
		*/
		break;
	    case F90_INTEGER_STAR_8:
		if (flag == 0) fprintf((FILE *) *file_id, "%lu", * (unsigned long int *) buffer);
		else printf("%lu", * (unsigned long int *) buffer);
		break;
	    default:
		return;
	}
}



void wclose_(unsigned long int *file_id) {

/* 
 * This function works for close an opened file
 */

	fclose((FILE *) *file_id);
}
	

