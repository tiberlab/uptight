#include <stdio.h>
#include <stdlib.h>


#include <cuda.h>
#include "cusparse_v2.h"
#include "cublas_v2.h"
#include <cuda_runtime.h>
#include <curand.h>

#include "mkl.h"
//include mkl_lapacke.h // this does not work with mkl 10.2

#include <unistd.h>
#include <sys/unistd.h>
#include <ctype.h>

#include <time.h>
#include <sys/time.h>

#define cuDouble double

#define CLEANUP(s) \
do { \
printf ("%s\n", s); \
if (p_current) free(p_current);\
if (p_previous) free(p_previous);\
if (p_next) free(p_next);\
if (p_aux) free(p_aux);\
if (p_aux2) free(p_aux2);\
if (eigen_vec) free(eigen_vec);\
if (eigen_val) free(eigen_val);\
if (p_current_device) cudaFree(p_current_device);\
if (p_previous_device) cudaFree(p_previous_device);\
if (p_next_device) cudaFree(p_next_device);\
if (p_aux_device) cudaFree(p_aux_device);\
if (p_aux2_device) cudaFree(p_aux2_device);\
if (eigen_val_device) cudaFree(eigen_val_device);\
if (eigen_vec_device) cudaFree(eigen_vec_device);\
if (eigen_seed_device) cudaFree(eigen_seed_device);\
if (norm_eigen_seed_device) cudaFree(norm_eigen_seed_device);\
if (alpha_device) cudaFree(alpha_device);\
if (beta_device) cudaFree(beta_device);\
if (temp1_device) cudaFree(temp1_device);\
if (temp_1_device) cudaFree(temp_1_device);\
if (temp_2_device) cudaFree(temp_2_device);\
fflush (stdout); \
} while (0);

#define THREADS_PER_VECTOR  32
#define VECTORS_PER_BLOCK  (192/THREADS_PER_VECTOR)

#define WARP_SIZE  32 
#define BLOCK_SIZE 1024 
#define BLOCK_SIZE_MUL 256 

extern "C" {
  void *__gxx_personality_v0;

  void fast_lanczos_ev_cusparse_real_(int* Size_Mat, double* valptr, int* rowptr, int* colptr, char* sparse_fmt, 
      int* Min_step, int* Long_step, int* Max_step, double* eigen_seed, int*  N_ham, double* Fast_tol, 
      double* Long_tol, double* Ort_tol_in, int* COUNTER_inout, int* res_flag, double* eigen_v, int* NR_eigv, 
      double* Energy, double* DeltaE);
}

void project_out(cublasHandle_t cubl_handle, cuDouble* vec, cuDouble* eigen_v_device, int nr_eigv, 
                 int n_ham, int numBlocks, int threadsPerBlock);

// v = 0
__global__ void vct_init_zero_kernel(cuDouble* vct, int n_ham);
// v2 = v1
__global__ void cpy_vct1_to_vct2_kernel(cuDouble* vct_scr, cuDouble* vct_des, int n_ham);


// ---------------------------------------------------------------------------------------------------------------------
// BLAS KERNELS
// ---------------------------------------------------------------------------------------------------------------------
// v1 = v1 + v2
__global__ void vct_pls_vct_kernel(cuDouble* vct1, cuDouble* vct2, int n_ham);

// v1 = s* v2
__global__ void scl_mul_vct_kernel(cuDouble* vct1, cuDouble* vct2, int n_ham, double scalar);

// v1 = v1 + s* v2
__global__ void vct_pls_scl_mul_vct_kernel(cuDouble* vct1, cuDouble* vct2, double scalar, int n_ham);


__global__ void project_out_kernel(cuDouble* vec, cuDouble* eigen_vector, cuDouble* alpha, int n_ham);



void fast_lanczos_ev_cusparse_real(int* Size_Mat, cuDouble* valptr, int* rowptr, int* colptr, char* sparse_fmt, 
                                   int* Min_step, int* Long_step, int* Max_step, cuDouble* eigen_seed, int* N_ham, 
                                   double* Fast_tol, double* Long_tol, double* Ort_tol_in, int* COUNTER_inout, int* res_flag, 
                                   cuDouble* eigen_v, int* NR_eigv, double* Energy, double* DeltaE) 
{

  // internal variables copies
  int min_step = *Min_step;
  int long_step = *Long_step;
  int max_step = *Max_step;
  int n_ham = *N_ham;
  int nr_eigv = *NR_eigv;
  double fast_tol = *Fast_tol;
  double long_tol = *Long_tol;
  double ort_tol_in = *Ort_tol_in;
  int COUNTER = *COUNTER_inout;
  int size_mat = *Size_Mat;
  
  // output variables
  double energy, deltaE;

  //----------Local Variables-------------------------------------------
  
  const cuDouble  scalar_1 = 1.0000;
  const cuDouble  scalar_2 = 0.0000;

  cuDouble* alpha = (cuDouble *) malloc(sizeof(cuDouble)); 
  cuDouble* beta = (cuDouble *) malloc(sizeof(cuDouble)); 
  cuDouble* alpha_device, *beta_device;

  // Tolerance test variables
  double tol, errors;
  double test_E;
  double ort_tol = ort_tol_in;
  double ort_chk = 0;
  int test_fast, test;
  // Counters
  int counter, j, nstep;
  
  cuDouble* temp1 = (cuDouble *) malloc(sizeof(cuDouble));
  cuDouble* temp1_device;
  cuDouble* temp_1, * temp_2, * temp_1_device, * temp_2_device;
  curandGenerator_t gen;

  // Matrix on device
  cuDouble* valptr_device;
  int* colptr_device;
  int* rowptr_device;
  cuDouble* p_current;
  cuDouble* p_previous;
  cuDouble* p_next;
  cuDouble* p_aux;
  cuDouble* p_aux2;
  cuDouble* p_current_device;
  cuDouble* p_previous_device;
  cuDouble* p_next_device;
  cuDouble* p_aux_device;
  cuDouble* p_aux2_device;

  cuDouble* eigen_seed_device;
  cuDouble* eigen_v_device;
  double* norm_eigen_seed = (double*) malloc(sizeof(double));
  double* norm_eigen_seed_device;

  p_current = (cuDouble*) malloc(n_ham*sizeof(cuDouble));
  memset(p_current, 0, n_ham*sizeof(cuDouble));
  p_previous = (cuDouble*) malloc(n_ham*sizeof(cuDouble));
  memset(p_previous, 0, n_ham*sizeof(cuDouble));
  p_next = (cuDouble*) malloc(n_ham*sizeof(cuDouble));
  memset(p_next, 0, n_ham*sizeof(cuDouble));
  p_aux = (cuDouble*) malloc(n_ham*sizeof(cuDouble));
  memset(p_aux, 0, n_ham*sizeof(cuDouble));
  p_aux2 = (cuDouble*) malloc(n_ham*sizeof(cuDouble));
  memset(p_aux2, 0, n_ham*sizeof(cuDouble));

  // LAPACK diagonalizer variables-----------------------------------------
  double* eigen_val;
  double* eigen_val_device;
  double* eigen_vec;
  double* eigen_vec_device;
  double* ad1;
  double* ad2;
  double* energy_test;
  double* work;
  double temp;
  int ldz = max_step;
  char JOBZ = 'V';
  char RANGE = 'I';
  MKL_INT IL = 1;
  MKL_INT IU = 1;
  MKL_INT LZ=  ldz;
  double VL, VU;
  MKL_INT* ifail;
  MKL_INT* iwork;
  MKL_INT M_out;
  MKL_INT error_code=0;
  double ABSTOL = 1e-12;
  eigen_val = (double *) malloc(2*sizeof(double));
  memset(eigen_val,0,2*sizeof(double));
  eigen_vec = (double *) malloc(ldz*2*sizeof(double));
  memset(eigen_vec,0,ldz*2*sizeof(double));
  work = (double *) malloc( 5*ldz*sizeof(double));
  memset(work,0,ldz*5*sizeof(double));
  iwork = (MKL_INT *) malloc( 5*ldz*sizeof(MKL_INT));
  memset(iwork,0,ldz*5*sizeof(MKL_INT));
  ifail = (MKL_INT *) malloc(ldz*sizeof(MKL_INT));
  memset(ifail,0,ldz*sizeof(MKL_INT));
  ad1 = (double *) malloc(ldz*sizeof(double));
  memset(ad1,0,ldz*sizeof(double));
  ad2 = (double *) malloc(ldz*sizeof(double));
  memset(ad2,0,ldz*sizeof(double));
  energy_test = (double *) malloc(10*sizeof(double));
  memset(energy_test,0,10*sizeof(double));
  //printf("Memory required for lanczos vectors %d bytes\n", n_ham*5*16);
  //printf("Memory required for lanczos matrix %d bytes hello\n", ldz*13*8);
  //------------------------------------------------------------------------

  //$$$$$$$$$$$$$$$$$$$$$$$$$$$ cusparse
  int threadsPerBlock = 1024;
  int numBlocks=(n_ham/threadsPerBlock)+1;
  //int numBlocks=8;

  //int * Min_step, int * Long_step, int * Max_step, double * Fast_tol, double * Long_tol, double * Ort_tol_in, int * COUNTER_inout
  //printf("Min_step = %d\n", min_step);
  //printf("Long_step = %d\n", long_step);
  //printf("Max_step = %d\n", max_step);
  //printf("Fast_tol = %lf\n", fast_tol);
  //printf("Long_tol = %lf\n", long_tol);
  //printf("Ort_tol_in = %lf\n", ort_tol_in);
  //printf("Counter_inout = %d\n", COUNTER);

  cusparseStatus_t cusp_status;
  cusparseHandle_t cusp_handle=0;
  cusparseMatDescr_t cusp_descra=0;
  cublasStatus_t cubl_status;
  cublasHandle_t cubl_handle=0;
  cudaDeviceReset();
  cudaThreadExit();
  cudaSetDevice(0);
  /* create a sparse and dense vector */

  temp_1 = (cuDouble *) malloc(sizeof(cuDouble));
  temp_2 = (cuDouble *) malloc(sizeof(cuDouble));

  cudaError_t cudaStat[23];
  for (counter=1; counter < 23; counter++) cudaStat[counter] = cudaSuccess;

  cudaStat[1]  = cudaMalloc((void**) &colptr_device, size_mat*sizeof(int));
  cudaStat[2]  = cudaMalloc((void**) &rowptr_device, (n_ham+1)*sizeof(int));
  cudaStat[3]  = cudaMalloc((void**) &valptr_device, size_mat*sizeof(cuDouble));
  cudaStat[4]  = cudaMalloc((void**) &p_current_device, n_ham*sizeof(cuDouble));
  cudaStat[5]  = cudaMalloc((void**) &p_previous_device, n_ham*sizeof(cuDouble));
  cudaStat[6]  = cudaMalloc((void**) &p_next_device, n_ham*sizeof(cuDouble));
  cudaStat[7]  = cudaMalloc((void**) &p_aux_device, n_ham*sizeof(cuDouble));
  cudaStat[8]  = cudaMalloc((void**) &p_aux2_device, n_ham*sizeof(cuDouble));
  cudaStat[9] = cudaMalloc((void**) &eigen_seed_device, n_ham*sizeof(cuDouble));
  cudaStat[10] = cudaMalloc((void**) &eigen_v_device, n_ham*nr_eigv*sizeof(cuDouble));
  cudaStat[14] = cudaMalloc((void**) &eigen_vec_device, ldz*2*sizeof(double));
  cudaStat[15] = cudaMalloc((void**) &norm_eigen_seed_device, sizeof(double));
  cudaStat[16] = cudaMalloc((void**) &alpha_device, sizeof(cuDouble));
  cudaStat[17] = cudaMalloc((void**) &beta_device, sizeof(cuDouble));
  cudaStat[18] = cudaMalloc((void**) &temp1_device, sizeof(cuDouble));
  cudaStat[20] = cudaMalloc((void**) &temp_1_device, sizeof(cuDouble));
  cudaStat[21] = cudaMalloc((void**) &temp_2_device, sizeof(cuDouble));

  for (counter=1; counter < 23; counter++)
    if (cudaStat[counter] != cudaSuccess){CLEANUP("Device malloc failed\n")};

  test = 0; // FALSE
  nstep = min_step;

  for(counter=0; counter < n_ham; counter++)
  {
    eigen_seed[counter] = eigen_seed[counter]+(rand()%100);
  }

  cudaStat[1] = cudaMemcpy(colptr_device, colptr, (size_t)(size_mat*sizeof(int)), cudaMemcpyHostToDevice);
  cudaStat[2] = cudaMemcpy(rowptr_device, rowptr, (size_t)((n_ham+1)*sizeof(int)), cudaMemcpyHostToDevice);
  cudaStat[3] = cudaMemcpy(valptr_device, valptr, (size_t)(size_mat*sizeof(cuDouble)), cudaMemcpyHostToDevice);
  cudaStat[4] = cudaMemcpy(eigen_seed_device, eigen_seed, (size_t)(n_ham*sizeof(cuDouble)), cudaMemcpyHostToDevice);
  cudaStat[5] = cudaMemcpy(eigen_v_device, eigen_v, (size_t)(n_ham*nr_eigv*sizeof(cuDouble)), cudaMemcpyHostToDevice);

  for (counter=1; counter <= 5; counter++)
    if (cudaStat[counter] != cudaSuccess){ CLEANUP("Device malloc failed\n"); }

  // initalization of CUBLAS library
  cubl_status = cublasCreate(&cubl_handle);
  if (cubl_status != CUBLAS_STATUS_SUCCESS){ CLEANUP("CUBLAS initilization failed \n");}

  cubl_status = cublasSetPointerMode(cubl_handle, CUBLAS_POINTER_MODE_DEVICE);
  if (cubl_status != CUBLAS_STATUS_SUCCESS){ CLEANUP("CUBLAS setting device pointer mode failed\n");}

  /* initialize cusparse library */
  cusp_status= cusparseCreate(&cusp_handle);
  if (cusp_status != CUSPARSE_STATUS_SUCCESS){ CLEANUP("CUSPARSE Library initialization failed\n");}

  /* create and setup matrix descriptor */
  cusp_status= cusparseCreateMatDescr(&cusp_descra);
  if (cusp_status != CUSPARSE_STATUS_SUCCESS){ CLEANUP("Matrix descriptor initialization failed 1");}

  cusp_status=cusparseSetMatType(cusp_descra, CUSPARSE_MATRIX_TYPE_GENERAL);
  if (cusp_status != CUSPARSE_STATUS_SUCCESS){ CLEANUP("Matrix descriptor initialization failed 2");}
  
  cusp_status=cusparseSetMatIndexBase(cusp_descra, CUSPARSE_INDEX_BASE_ONE);
  if (cusp_status != CUSPARSE_STATUS_SUCCESS){ CLEANUP("Matrix descriptor initialization failed 3");}

  cubl_status = cublasDznrm2(cubl_handle, n_ham, eigen_seed_device, 1, norm_eigen_seed_device);

  if(cubl_status != CUBLAS_STATUS_SUCCESS){ CLEANUP("CUBLAS initilization failed \n");}

  cudaStat[1] = cudaMemcpy(norm_eigen_seed, norm_eigen_seed_device,(size_t)(sizeof(double)), cudaMemcpyDeviceToHost);
  if (cudaStat[1] != cudaSuccess){CLEANUP("copy to host failed 1\n"); }

  scl_mul_vct_kernel<<<numBlocks, threadsPerBlock>>>(eigen_seed_device, eigen_seed_device, n_ham, 1.0/(*norm_eigen_seed));

  if (nr_eigv > 1) project_out(cubl_handle, eigen_seed_device, eigen_v_device, nr_eigv-1, n_ham, 
                                                                                         numBlocks, threadsPerBlock);

  // LOOP ////////////////////////////////////////////////////////////////////////////////////////////////////////////
       //"--------+-----+-----+------------------+------------------+-------------------+-------------------+"
  printf("         niters nsteps       test_E             energy              error            ort_chk\n"); 


  while(!test)
  {
 
    alpha[0] = 0.0;
    beta[0] = 0.0;

    vct_init_zero_kernel<<<numBlocks, threadsPerBlock>>>( p_previous_device, n_ham);
    vct_init_zero_kernel<<<numBlocks, threadsPerBlock>>>( p_next_device, n_ham);
    vct_init_zero_kernel<<<numBlocks, threadsPerBlock>>>( p_aux_device, n_ham);
    cpy_vct1_to_vct2_kernel<<<numBlocks, threadsPerBlock>>>( eigen_seed_device, p_current_device, n_ham);

    test_fast = 0; 

    // Lanczos Loop	
    for(j = 0; j < nstep; j++)
    {

      // update | J + 1 > to  | J + 1 > + A * A | j >
      //-------------------------------------------------------------------
      cusp_status= cusparseDcsrmv(cusp_handle,CUSPARSE_OPERATION_NON_TRANSPOSE, n_ham, n_ham, size_mat, 
                                    &scalar_1, cusp_descra, valptr_device, rowptr_device, colptr_device, 
                                    p_current_device, &scalar_2, p_aux_device);

      if (cusp_status != CUSPARSE_STATUS_SUCCESS){ CLEANUP("Matrixvector multiplication failed 1 \n"); }

      cusp_status= cusparseDcsrmv(cusp_handle,CUSPARSE_OPERATION_NON_TRANSPOSE, n_ham, n_ham, size_mat, 
                                  &scalar_1, cusp_descra, valptr_device, rowptr_device, colptr_device, p_aux_device, 
                                  &scalar_2, p_aux2_device);

      if (cusp_status != CUSPARSE_STATUS_SUCCESS){ CLEANUP("Matrixvector multiplication failed 2"); }


      vct_pls_vct_kernel<<<numBlocks, threadsPerBlock>>>( p_next_device, p_aux2_device, n_ham);
 
      // alpha = < j | A * A | j >
      //-------------------------------------------------------------------
      cubl_status = cublasDdotc(cubl_handle, n_ham, p_aux_device, 1, p_aux_device, 1, alpha_device);
      if(cubl_status != CUBLAS_STATUS_SUCCESS){ CLEANUP("CUBLAS dot product failed hello1  %d\n"); }

      cudaStat[1] = cudaMemcpy(alpha, alpha_device,(size_t)(sizeof(cuDouble)), cudaMemcpyDeviceToHost);
      if (cudaStat[1] != cudaSuccess){ CLEANUP("copy to host failed 4\n"); }

      //  update | J + 1 > to | J + 1 > - alpha | j >
      //-------------------------------------------------------------------    
      vct_pls_scl_mul_vct_kernel<<<numBlocks, threadsPerBlock>>>( p_next_device, p_current_device, -alpha[0], n_ham);
 
      //  project out existing eigenvectors
      //-------------------------------------------------------------------    
      if (nr_eigv > 1) project_out(cubl_handle, p_next_device, eigen_v_device, nr_eigv-1, n_ham, 
                                                                                         numBlocks, threadsPerBlock);
      // beta = SQRT( < J + 1 | J + 1 > )
      //------------------------------------------------------------------- 
      cubl_status = cublasDdotc(cubl_handle, n_ham, p_next_device, 1, p_next_device, 1, beta_device);
      if(cubl_status != CUBLAS_STATUS_SUCCESS){ CLEANUP("CUBLAS dot product failed hello2\n");}

      cudaStat[1] = cudaMemcpy(beta, beta_device,(size_t)(sizeof(cuDouble)), cudaMemcpyDeviceToHost);
      if (cudaStat[1] != cudaSuccess){ CLEANUP("copy to host failed 5\n"); }

      beta[0] = sqrt(beta[0]);
      //------------------------------------------------------------------- 
      // orthogonality check
      if (j == nstep-1 && nstep > 2 && nstep < max_step-10)
      {
        cubl_status = cublasDdotc(cubl_handle, n_ham, p_next_device, 1, eigen_seed_device, 1, temp1_device);
        if(cubl_status != CUBLAS_STATUS_SUCCESS){CLEANUP("CUBLAS dot product failed 3\n");}

        cudaStat[1] = cudaMemcpy(temp1, temp1_device,(size_t)(sizeof(cuDouble)), cudaMemcpyDeviceToHost);
        if (cudaStat[1] != cudaSuccess){ CLEANUP("copy to host failed 6\n");}

        ort_chk = fabs(temp1[0]/beta[0]);

        if(ort_chk < ort_tol)  nstep = nstep + 10;
      }
      //------------------------------------------------------------------- 
      ad2[j] = alpha[0];
      ad1[j] = beta[0];

      //------------------------------------------------------------------- 
      // update |j - 1> to | j >
      cpy_vct1_to_vct2_kernel<<<numBlocks, threadsPerBlock>>>( p_current_device, p_previous_device, n_ham);
      // update |j> to |j + 1>/beta
      scl_mul_vct_kernel<<<numBlocks, threadsPerBlock>>>( p_current_device, p_next_device, n_ham, 1.0/beta[0]);
      // initialize |j + 1> = -beta |j> (= -beta |j-1>)
      scl_mul_vct_kernel<<<numBlocks, threadsPerBlock>>>( p_next_device, p_previous_device, n_ham, -beta[0]);

    } // end of lanczos loop

    //**********************************************************************************************************************
    COUNTER = COUNTER + nstep;

    //clock_gettime(CLOCK_REALTIME, &tstart);

    //diagonalize the T matrix
    //NOTE: mkl 10.2 does not define LAPACKE so using standard version
    //error_code = LAPACKE_dstevx( LAPACK_COL_MAJOR, 'V', 'I', nstep, ad2, ad1, VL, VU, IL, IU, ABSTOL, 
    //                          &M_out, eigen_val, eigen_vec, LZ, ifail );

    dstevx( &JOBZ, &RANGE, &nstep, ad2, ad1, &VL, &VU, &IL, &IU, &ABSTOL, 
           &M_out,eigen_val,eigen_vec, &LZ, work, iwork, ifail,&error_code);

    //printf("\n no of eigen values found is %d, %d\n", M_out, error_code);
    //printf("\n eigen values found is %f\n", eigen_val[0]);

    cudaStat[1] = cudaMemcpy(eigen_vec_device, eigen_vec, (size_t)(2*ldz*sizeof(double)), cudaMemcpyHostToDevice);
    if (cudaStat[1] != cudaSuccess){ CLEANUP("copy to device failed \n"); }

    energy = sqrt(fabs(eigen_val[0]));

    double mean_en = 0.0;
    for(counter=0; counter < 9; counter++) 
    {
      energy_test[counter]= energy_test[counter+1];
      mean_en += energy_test[counter];
    }
    energy_test[9]=energy;
    mean_en += energy;  mean_en /= 10.0;

    temp=0.0;
    for(counter=0; counter < 10; counter++)
      temp += pow((energy_test[counter] - mean_en),2);

    temp=temp/10.0;
    double dev_en = sqrt(temp);

    if(nstep == min_step)
    {
      if((dev_en/mean_en) < fast_tol) test_fast = 1; 
    }

    //********************************************************************************************************************
    // SECOND LOOP: project Ritz eigenvector. 
    //********************************************************************************************************************
    // RESTART ITERATION FROM  eigen_seed 
    // |j> = |seed>
    cpy_vct1_to_vct2_kernel<<<numBlocks, threadsPerBlock>>>( eigen_seed_device, p_current_device, n_ham);

    vct_init_zero_kernel<<<numBlocks, threadsPerBlock>>>( p_previous_device, n_ham);
    vct_init_zero_kernel<<<numBlocks, threadsPerBlock>>>( p_next_device, n_ham);
    vct_init_zero_kernel<<<numBlocks, threadsPerBlock>>>( p_aux_device, n_ham);
    vct_init_zero_kernel<<<numBlocks, threadsPerBlock>>>( eigen_seed_device, n_ham);

    alpha[0] = 0.0;
    beta[0] = 0.0;
    
    for( j= 0; j < nstep; j++)
    { 
      // |aux> = A |j> 
      cusp_status = cusparseDcsrmv(cusp_handle,CUSPARSE_OPERATION_NON_TRANSPOSE, n_ham, n_ham, size_mat, 
                                   &scalar_1, cusp_descra, valptr_device, rowptr_device, colptr_device, p_current_device, 
                                   &scalar_2, p_aux_device);

      if (cusp_status != CUSPARSE_STATUS_SUCCESS){ CLEANUP("Matrixvector multiplication failed"); }

      // |aux2> = A |aux> 
      cusp_status= cusparseDcsrmv(cusp_handle,CUSPARSE_OPERATION_NON_TRANSPOSE, n_ham, n_ham, size_mat, 
                                  &scalar_1, cusp_descra, valptr_device, rowptr_device, colptr_device, p_aux_device, 
                                  &scalar_2, p_aux2_device);

      if (cusp_status != CUSPARSE_STATUS_SUCCESS){ CLEANUP("Matrixvector multiplication failed"); }

      // |j+1> = |j+1> + |aux2> = |j+1> + A^2 |j> 
      vct_pls_vct_kernel<<<numBlocks, threadsPerBlock>>>( p_next_device, p_aux2_device, n_ham);

      //  alpha = <j| A A |j> = <j | j+1>  
      cubl_status = cublasDdotc(cubl_handle, n_ham, p_aux_device, 1, p_aux_device, 1, alpha_device);
      if(cubl_status != CUBLAS_STATUS_SUCCESS){ CLEANUP("CUBLAS dot product failed hello4\n"); }

      cudaStat[1] = cudaMemcpy(alpha, alpha_device,(size_t)(sizeof(cuDouble)), cudaMemcpyDeviceToHost);
      if (cudaStat[1] != cudaSuccess){ CLEANUP("copy to host failed 7\n");}

      // |j+1> = |j+1> - alpha |j>
      vct_pls_scl_mul_vct_kernel<<<numBlocks, threadsPerBlock>>>( p_next_device, p_current_device, -alpha[0], n_ham);

      if (nr_eigv > 1) project_out(cubl_handle, p_next_device, eigen_v_device, nr_eigv-1, n_ham, 
                                                                                         numBlocks, threadsPerBlock);

      // beta = <j+1|j+1> 
      cubl_status = cublasDdotc(cubl_handle, n_ham, p_next_device, 1, p_next_device, 1, beta_device);
      if (cubl_status != CUBLAS_STATUS_SUCCESS){ CLEANUP("CUBLAS dot product failed hello5\n");}

      cudaStat[1] = cudaMemcpy(beta, beta_device,(size_t)(sizeof(cuDouble)), cudaMemcpyDeviceToHost);
      if (cudaStat[1] != cudaSuccess){ CLEANUP("copy to host failed 8\n"); }

      beta[0] = sqrt(beta[0]);

      // Ritz Projection: |u> = |u> + Vj*|j> 
      //vct_1_add_vct_2__mul_vct_3_kernel<<<numBlocks, threadsPerBlock>>>( eigen_seed_device, eigen_vec_device, 
      //                                                                    p_current_device, n_ham, j);

      vct_pls_scl_mul_vct_kernel<<<numBlocks, threadsPerBlock>>>( eigen_seed_device,  p_current_device, 
                                                                  eigen_vec[j], n_ham);
      // |j-1>  <- |j> 
      cpy_vct1_to_vct2_kernel<<<numBlocks, threadsPerBlock>>>( p_current_device, p_previous_device, n_ham);

      // |j>  <- |j+1>/beta 
      scl_mul_vct_kernel<<<numBlocks, threadsPerBlock>>>( p_current_device, p_next_device, n_ham, 1.0/beta[0]);

      // Initialize |j+1> = - beta |j>  ( = -beta |j-1> )
      scl_mul_vct_kernel<<<numBlocks, threadsPerBlock>>>( p_next_device, p_previous_device, n_ham, -beta[0]);

    } // end second loop


    if (nr_eigv > 1) project_out(cubl_handle, eigen_seed_device, eigen_v_device, nr_eigv-1, n_ham, 
                                                                                         numBlocks, threadsPerBlock);
    // temp_1 = <u|u> 
    cubl_status = cublasDznrm2(cubl_handle, n_ham, eigen_seed_device, 1, norm_eigen_seed_device);

    cudaStat[1] = cudaMemcpy(norm_eigen_seed, norm_eigen_seed_device,(size_t)(sizeof(double)), cudaMemcpyDeviceToHost);
    if (cudaStat[1] != cudaSuccess){ CLEANUP("copy to host failed 9\n"); }

    // normalize |u>  
    scl_mul_vct_kernel<<<numBlocks, threadsPerBlock>>>(eigen_seed_device, eigen_seed_device, n_ham, 1.0/(*norm_eigen_seed));

    cudaStat[1] = cudaMemcpy(eigen_seed, eigen_seed_device, (size_t)(n_ham*sizeof(cuDouble)), cudaMemcpyDeviceToHost);
    if (cudaStat[1] != cudaSuccess){ CLEANUP("copy to device failed 2\n");}

    //************************************************************************************************************************
 
    // CALCULATION of <u| H |u>    
    //************************************************************************************************************************
    // |v> = H |u> 

    cusp_status= cusparseDcsrmv(cusp_handle,CUSPARSE_OPERATION_NON_TRANSPOSE, n_ham, n_ham, size_mat,
                                &scalar_1, cusp_descra, valptr_device, rowptr_device, colptr_device, eigen_seed_device, 
                                &scalar_2, p_aux_device);

    // temp_2 = <u|v>
    cubl_status = cublasDdotc(cubl_handle, n_ham, eigen_seed_device, 1, p_aux_device, 1, temp_2_device);
    cudaStat[1] = cudaMemcpy(temp_2, temp_2_device,(size_t)(sizeof(cuDouble)), cudaMemcpyDeviceToHost);
    if (cudaStat[1] != cudaSuccess){ CLEANUP("copy to host failed 10\n");}

    test_E = temp_2[0];

    deltaE = fabs(energy - (fabs(test_E)));

    if(nstep == min_step)
    {
       printf("cuspFast %5d %5d %18.8g %18.8g %18.8g %18.8g\n", 
               COUNTER, nstep, test_E, energy, deltaE, ort_chk);
    }
    else
    {
       printf("cuspStan %5d %5d %18.8g %18.8g %18.8g %18.8g\n", 
               COUNTER, nstep, test_E, energy, deltaE, ort_chk);
    }
  
    if(deltaE < long_tol)
    {
      test = 1; // true
      printf("Eig value = %f, Error = %f\n", energy, deltaE);
      *Energy = test_E;
      *DeltaE = deltaE;
      *COUNTER_inout = COUNTER;
    }
    else
    {
      test = 0; // false
    }

    if(nstep == min_step)
    {
      if(test_fast)
      {
        test_fast = 0; // false
        nstep = long_step;
      }
      else
      {
        nstep = min_step;
      }
    }

  } // while (!test)


}

// |v> = |v> - Sum_k |e_k><e_k|v>
void project_out(cublasHandle_t cubl_handle, cuDouble* vec, cuDouble* eigen_v_device, int nr_eigv, 
                 int n_ham, int numBlocks, int threadsPerBlock)
{
  
  cublasStatus_t cubl_status;
  
  //cuDouble* val = (cuDouble *) malloc(sizeof(cuDouble));
 
  cuDouble* alpha;  
  cudaMalloc((void**) &alpha, sizeof(cuDouble));

  for (unsigned int k=0; k < nr_eigv; k++)
  {
     cubl_status = cublasDdotc(cubl_handle, n_ham,  &eigen_v_device[k*n_ham], 1, vec, 1, alpha);
     project_out_kernel<<<numBlocks, threadsPerBlock>>>(vec, &eigen_v_device[k*n_ham], alpha, n_ham);
  }

  cudaFree(alpha);
}


__global__ void scl_mul_vct_kernel(cuDouble* v1, cuDouble* v2, int vct_size, double slr)
{
  int ID;
  //cuDouble * dptr;
  ID = blockDim.x * blockIdx.x + threadIdx.x;
  //dptr = (cuDouble *)((char *)V+j*pitch);
  //for(m=0; m < Vct_Size; m++){
  if(ID < vct_size)
  {
    v1[ID] = slr * v2[ID];
  }
}



__global__ void vct_init_zero_kernel(cuDouble* vct, int n_ham)
{
  int ID;
  ID = blockDim.x * blockIdx.x + threadIdx.x;
  if(ID < n_ham)
  {
    vct[ID] = 0.00000000000000;
  }
}


__global__ void cpy_vct1_to_vct2_kernel(cuDouble* vct_scr, cuDouble* vct_des, int n_ham)
{
  int ID;
  ID = blockDim.x * blockIdx.x + threadIdx.x;
  if(ID < n_ham)
  {
    vct_des[ID] = vct_scr[ID];
  }
}


__global__ void vct_pls_vct_kernel(cuDouble* vct_1, cuDouble* vct_2, int n_ham)
{
  int ID;
  ID = blockDim.x * blockIdx.x + threadIdx.x;
  if(ID < n_ham)
  {
    vct_1[ID] = vct_1[ID] + vct_2[ID];
  }
}



__global__ void vct_pls_scl_mul_vct_kernel(cuDouble* vct1, cuDouble* vct2, double scalar, int n_ham)
{
  int ID;
  ID = blockDim.x * blockIdx.x + threadIdx.x;
  if(ID < n_ham)
  {
    vct1[ID] =  vct1[ID] + scalar * vct2[ID];
  }
}


__global__ void vct_pls_scl_mul_vct_kernel(cuDouble* vct1, cuDouble* vct2, cuDouble scalar, int n_ham)
{
  int ID;
  ID = blockDim.x * blockIdx.x + threadIdx.x;
  if(ID < n_ham)
  {
    vct1[ID] =  vct1[ID] + scalar * vct2[ID];
  }
}



__global__ void project_out_kernel(cuDouble* vec, cuDouble* eigen_vector, cuDouble* alpha, int n_ham)
{

  int  ID;
  ID = blockDim.x * blockIdx.x + threadIdx.x;

  if (ID < n_ham)
  {
    vec[ID] = vec[ID] - alpha[0] * eigen_vector[ID];
  }
}




