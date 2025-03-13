#include <cstdio>
#include <type_traits>
//#include <stdlib.h>


#include <cuda.h>
#include "cusparse_v2.h"
#include "cublas_v2.h"
#include <cuda_runtime.h>
#include <curand.h>

#include "mkl.h"
//include mkl_lapacke.h // this does not work with mkl 10.2

//#include <unistd.h>
//#include <sys/unistd.h>
//#include <ctype.h>

//#include <time.h>
//#include <sys/time.h>


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
  void fast_lanczos_ev_mod_(int* Size_Mat, cuDoubleComplex* valptr, int* rowptr, int* colptr, char* sparse_fmt, 
      int* Min_step, int* Long_step, int* Max_step, cuDoubleComplex* eigen_seed, int* N_ham, double* Fast_tol, 
      double* Long_tol, double* Ort_tol_in, int* COUNTER_inout, int* res_flag, cuDoubleComplex* eigen_v, 
      int* NR_eigv, double* Energy, double* DeltaE);

  void fast_lanczos_ev_cusparse_(int* Size_Mat, cuDoubleComplex* valptr, int* rowptr, int* colptr, char* sparse_fmt, 
      int* Min_step, int* Long_step, int* Max_step, cuDoubleComplex* eigen_seed, int*  N_ham, double* Fast_tol, 
      double* Long_tol, double* Ort_tol_in, int* COUNTER_inout, int* res_flag, cuDoubleComplex* eigen_v, int* NR_eigv, 
      double* Energy, double* DeltaE);

  void fast_lanczos_ev_cusparse_split_(int* Size_Mat_real, int* Size_Mat_img, double* valptr_real, int* rowptr_real, 
      int* colptr_real, double* valptr_img, int* rowptr_img, int* colptr_img, char* sparse_fmt, int* Min_step, 
      int* Long_step, int* Max_step, cuDoubleComplex* eigen_seed, int* N_ham, double* Fast_tol, 
      double* Long_tol, double* Ort_tol_in, int* COUNTER_inout, int* res_flag, cuDoubleComplex* eigen_v, 
      int* NR_eigv, double* Energy, double* DeltaE);

}

void spmv(int* rowptr_device_real, int*  colptr_device_real, double* valptr_device_real, int size_mat_real, 
         int* rowptr_device_img, int*  colptr_device_img, double* valptr_device_img, int size_mat_img,
         cuDoubleComplex* in, cuDoubleComplex* out, int n_ham,
         cusparseHandle_t&, cudaStream_t& stream1, cudaStream_t& stream2, cudaStream_t& stream3, cudaStream_t& stream4);

void project_out(cublasHandle_t cubl_handle, cuDoubleComplex* vec, cuDoubleComplex* eigen_v_device, int nr_eigv, 
                 int n_ham, int numBlocks, int threadsPerBlock);

// Set the GPU id to be used
/*
 * Can take the ID from the argument, or from the environment if gpu == -1
 */
static void set_gpu_id(int gpu = -1);



// v = 0
__global__ void vct_init_zero_kernel(cuDoubleComplex* vct, int n_ham);
// v2 = v1
__global__ void cpy_vct1_to_vct2_kernel(cuDoubleComplex* vct_scr, cuDoubleComplex* vct_des, int n_ham);
// vr = v.x; vi = v.y
__global__ void split_vct_kernel( cuDoubleComplex* scr_vct, double* dest_real, double* dest_img, int n_ham);
// v = vr + i vi
__global__ void form_complex_vct_kernel( double* real_a, double*  real_b, double* img_a, double* img_b, 
                                         cuDoubleComplex* result, int n_ham);

// ---------------------------------------------------------------------------------------------------------------------
// BLAS KERNELS
// ---------------------------------------------------------------------------------------------------------------------
// v1 = v1 + v2
__global__ void vct_pls_vct_kernel(cuDoubleComplex* vct1, cuDoubleComplex* vct2, int n_ham);

// v1 = s* v2
__global__ void scl_mul_vct_kernel(cuDoubleComplex* vct1, cuDoubleComplex* vct2, int n_ham, double scalar);
__global__ void scl_mul_vct_kernel(cuDoubleComplex* vct1, cuDoubleComplex* vct2, int n_ham, cuDoubleComplex scalar);


// v1 = v1 + s* v2
__global__ void vct_pls_scl_mul_vct_kernel(cuDoubleComplex* vct1, cuDoubleComplex* vct2, double scalar, int n_ham);
__global__ void vct_pls_scl_mul_vct_kernel(cuDoubleComplex* vct1, cuDoubleComplex* vct2, cuDoubleComplex scalar, int n_ham);

// ---------------------------------------------------------------------------------------------------------------------
// spMV KERNELS
// A v
// Implements naive kernel
__global__ void spmv_csr_simple_kernel(int num_rows, const int* ptr, const int* indices, const double* data, 
                                       const double* x, double* y); 

// Implement Mike Giles Kernel
__global__ void spmv_csr_vector_kernel(int num_rows, const int* ptr, const int* indices, const double* data, 
                                       const double* x, double* y, 
                                       int repeat, int coop);
 
// Mike Giles for Hybrid Real * cmplx operation 
__global__ void spmv_csr_hybrid_kernel(int num_rows, const int* ptr, const int* indices, const double* data, 
                                       const cuDoubleComplex* x, cuDoubleComplex*  y, int repeat, int coop);
__global__ void spmv_csr_hybrid_kernel(int num_rows, const int* rowPtrs, const int* colIdxs, const double* values, 
                                       const cuDoubleComplex* x, cuDoubleComplex* y);

__global__ void project_out_kernel(cuDoubleComplex* vec, cuDoubleComplex* eigen_vector, cuDoubleComplex* alpha, int n_ham);





void set_gpu_id(int gpu)
{

  if (gpu < 0)
  {
    // read it from the environment
    char* gpu_env = getenv("CUDA_SELECT_GPU");
    if (gpu_env != NULL)
      gpu = atoi(gpu_env);
    printf("gpu_env = %s, gpu = %d\n", gpu_env, gpu);
  }

  cudaSetDevice(gpu);
  // TODO add error checking
}





void fast_lanczos_ev_cusparse_split_(int* Size_Mat_real, int* Size_Mat_img, double* valptr_real, int* rowptr_real, int* colptr_real, 
                                        double* valptr_img, int* rowptr_img, int* colptr_img,  char* sparse_fmt, 
                                        int* Min_step, int* Long_step, int* Max_step, cuDoubleComplex* eigen_seed, int* N_ham, 
                                        double* Fast_tol, double* Long_tol, double* Ort_tol_in, int* COUNTER_inout, int* res_flag, 
                                        cuDoubleComplex*  eigen_v, int* NR_eigv, double* Energy, double* DeltaE) 
{

  // input copies
  int min_step = *Min_step;
  int long_step = *Long_step;
  int max_step = *Max_step;
  int n_ham = *N_ham;
  int nr_eigv = *NR_eigv;
  double fast_tol = *Fast_tol;
  double long_tol = *Long_tol;
  double ort_tol_in = *Ort_tol_in;
  int COUNTER = *COUNTER_inout;
  int size_mat_real = *Size_Mat_real;
  int size_mat_img = *Size_Mat_img;

  // outputs
  double energy, deltaE;

  //----------Local Variables-------------------------------------------

  const cuDoubleComplex jcmpx=make_cuDoubleComplex(0.0,1.0);   

  // counter variables
  int counter, j, nstep, info, ldz;
  // lanczos variables
  cuDoubleComplex* alpha, * beta;
  cuDoubleComplex* alpha_device, * beta_device;

  // Matrix real and imag parts
  double* valptr_device_real;
  int* colptr_device_real;
  int* rowptr_device_real;
  double* valptr_device_img;
  int* colptr_device_img;
  int* rowptr_device_img;

  cuDoubleComplex* p_current;
  cuDoubleComplex* p_previous;
  cuDoubleComplex* p_next;
  cuDoubleComplex* p_aux;
  cuDoubleComplex* p_aux2;

  double test_E; 
  double tol, errors;
  int test_fast, test;
  double* energy_test;
  double ort_tol = ort_tol_in;
  double ort_chk = 0;

  double* norm_eigen_seed = (double *) malloc(sizeof(double));
  double* norm_eigen_seed_device;

  double temp;
  cuDoubleComplex* temp1_device;
  cuDoubleComplex* temp_1_device, * temp_2_device;
  curandGenerator_t gen;

  cuDoubleComplex* temp1 = (cuDoubleComplex*) malloc(sizeof(cuDoubleComplex));
  cuDoubleComplex* temp_1 = (cuDoubleComplex *) malloc(sizeof(cuDoubleComplex));
  cuDoubleComplex* temp_2 = (cuDoubleComplex *) malloc(sizeof(cuDoubleComplex));

  cuDoubleComplex* p_current_device;
  cuDoubleComplex* p_previous_device;
  cuDoubleComplex* p_next_device;
  cuDoubleComplex* p_aux_device;
  cuDoubleComplex* p_aux1_device;
  cuDoubleComplex* p_aux2_device;
  cuDoubleComplex* eigen_seed_device;
  cuDoubleComplex* eigen_v_device; 

  cudaStream_t stream1, stream2, stream3, stream4, stream0, stream5, stream6, stream7;

  p_current = (cuDoubleComplex*) malloc(n_ham*sizeof(cuDoubleComplex));
  memset(p_current,0,n_ham*sizeof(cuDoubleComplex));
  p_previous = (cuDoubleComplex*) malloc(n_ham*sizeof(cuDoubleComplex));
  memset(p_previous,0,n_ham*sizeof(cuDoubleComplex));
  p_next = (cuDoubleComplex*) malloc(n_ham*sizeof(cuDoubleComplex));
  memset(p_next,0,n_ham*sizeof(cuDoubleComplex));
  p_aux = (cuDoubleComplex*) malloc(n_ham*sizeof(cuDoubleComplex));
  memset(p_aux,0,n_ham*sizeof(cuDoubleComplex));
  p_aux2 = (cuDoubleComplex*) malloc(n_ham*sizeof(cuDoubleComplex));
  memset(p_aux2,0,n_ham*sizeof(cuDoubleComplex));

  // LAPACK variables
  double* eigen_val;
  double* eigen_vec;
  double* eigen_val_device;
  double* eigen_vec_device;
  ldz = max_step;
  char JOBZ = 'V';
  char RANGE = 'I';
  MKL_INT IL = 1;
  MKL_INT IU = 1;
  double VL, VU;
  MKL_INT LZ = ldz;
  MKL_INT M_out;
  MKL_INT* ifail;
  MKL_INT error_code=0;
  MKL_INT* iwork;
  double* work;
  double* ad1;
  double* ad2;
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

  //$$$$$$$$$$$$$$$$$$$$$$$$$$$ cusparse
  int threadsPerBlock = BLOCK_SIZE;                         //blockSize
  int numBlocks=(n_ham/threadsPerBlock)+1;                //gridSize

  int coop= 8  ;    
  int repeat= 2 ;   // repeat = nrows * coop/ BLOCK_SIZE_MUL 
  int numBlocksMul =(n_ham*coop-1)/(repeat*BLOCK_SIZE_MUL)+1;  //gridSize
 
  
  printf(" (CUDA) fast_lanczos_ev_cusparse_real, threadsPerBlock %d \n", threadsPerBlock);
  printf(" (CUDA) gridSize  %d \n",numBlocks);
  printf(" (CUDA) fast_lanczos_ev_cusparse_real, threadsPerBlock %d \n", BLOCK_SIZE_MUL);
  printf(" (CUDA) gridSize for Mxv %d \n",numBlocksMul);

  cusparseStatus_t cusp_status;
  cusparseHandle_t cusp_handle=0;
  cusparseMatDescr_t cusp_descra=0;
  cublasStatus_t cubl_status;
  cublasHandle_t cubl_handle=0;
  cudaDeviceReset();
  cudaThreadExit();


  int devCount = 0;	
  cudaGetDeviceCount(&devCount);
  cudaSetDevice(0);

  /* create a sparse and dense vector */

  cudaStreamCreate(&stream0);
  cudaStreamCreate(&stream1);
  cudaStreamCreate(&stream2);
  cudaStreamCreate(&stream3);
  cudaStreamCreate(&stream4);
  cudaStreamCreate(&stream5);
  cudaStreamCreate(&stream6);
  cudaStreamCreate(&stream7);

  cudaError_t cudaStat[33];
  for (counter=1; counter < 33; counter++) cudaStat[counter] = cudaSuccess;

  cudaStat[1]  = cudaMalloc((void**) &colptr_device_real, size_mat_real*sizeof(int));
  cudaStat[2]  = cudaMalloc((void**) &rowptr_device_real, (n_ham+1)*sizeof(int));
  cudaStat[3]  = cudaMalloc((void**) &valptr_device_real, size_mat_real*sizeof(double));
  cudaStat[4]  = cudaMalloc((void**) &p_current_device, n_ham*sizeof(cuDoubleComplex));
  cudaStat[5]  = cudaMalloc((void**) &p_previous_device, n_ham*sizeof(cuDoubleComplex));
  cudaStat[6]  = cudaMalloc((void**) &p_next_device, n_ham*sizeof(cuDoubleComplex));
  cudaStat[7]  = cudaMalloc((void**) &p_aux_device, n_ham*sizeof(cuDoubleComplex));
  cudaStat[23] = cudaMalloc((void**) &p_aux1_device, n_ham*sizeof(cuDoubleComplex));
  cudaStat[8]  = cudaMalloc((void**) &p_aux2_device, n_ham*sizeof(cuDoubleComplex));
  cudaStat[10] = cudaMalloc((void**) &eigen_vec_device, ldz*2*sizeof(double));
  cudaStat[11] = cudaMalloc((void**) &colptr_device_img, size_mat_img*sizeof(int));
  cudaStat[12] = cudaMalloc((void**) &rowptr_device_img, (n_ham+1)*sizeof(int));
  cudaStat[13] = cudaMalloc((void**) &valptr_device_img, size_mat_img*sizeof(double));
  cudaStat[14] = cudaMalloc((void**) &eigen_seed_device, n_ham*sizeof(cuDoubleComplex));
  cudaStat[15] = cudaMalloc((void**) &norm_eigen_seed_device, sizeof(double));
  cudaStat[18] = cudaMalloc((void**) &temp1_device, sizeof(cuDoubleComplex));
  cudaStat[20] = cudaMalloc((void**) &temp_1_device, sizeof(cuDoubleComplex));
  cudaStat[21] = cudaMalloc((void**) &temp_2_device, sizeof(cuDoubleComplex));
  cudaStat[22] = cudaMalloc((void**) &eigen_v_device, nr_eigv*n_ham*sizeof(cuDoubleComplex));
  for (counter=1; counter < 33; counter++)
  { 
    if (cudaStat[counter] != cudaSuccess){CLEANUP("Device malloc failed\n");}
  }  



  test = 0; // FALSE
  nstep = min_step;

  for(counter=0; counter < n_ham; counter++){
    eigen_seed[counter].x = rand()%100;
    eigen_seed[counter].y = rand()%100;
  }

  cudaStat[1] = cudaMemcpyAsync(colptr_device_real, colptr_real, (size_t)(size_mat_real*sizeof(int)), cudaMemcpyHostToDevice, stream0);
  cudaStat[2] = cudaMemcpyAsync(rowptr_device_real, rowptr_real, (size_t)((n_ham+1)*sizeof(int)), cudaMemcpyHostToDevice, stream1);
  cudaStat[3] = cudaMemcpyAsync(valptr_device_real, valptr_real, (size_t)(size_mat_real*sizeof(double)), cudaMemcpyHostToDevice, stream2);
  cudaStat[4] = cudaMemcpyAsync(eigen_seed_device, eigen_seed, (size_t)(n_ham*sizeof(cuDoubleComplex)), cudaMemcpyHostToDevice, stream3);
  cudaStat[5] = cudaMemcpyAsync(eigen_v_device, eigen_v, (size_t)(nr_eigv*n_ham*sizeof(cuDoubleComplex)), cudaMemcpyHostToDevice, stream4);
  cudaStat[6] = cudaMemcpyAsync(colptr_device_img, colptr_img, (size_t)(size_mat_img*sizeof(int)), cudaMemcpyHostToDevice, stream5);
  cudaStat[7] = cudaMemcpyAsync(rowptr_device_img, rowptr_img, (size_t)((n_ham+1)*sizeof(int)), cudaMemcpyHostToDevice, stream6);
  cudaStat[8] = cudaMemcpyAsync(valptr_device_img, valptr_img, (size_t)(size_mat_img*sizeof(double)), cudaMemcpyHostToDevice, stream7);
  
  for (counter=1; counter < 9; counter++)
    if (cudaStat[counter] != cudaSuccess){CLEANUP("Device malloc failed\n");}

  cudaSetDeviceFlags(cudaDeviceMapHost);
  cudaHostAlloc((void **) &alpha, sizeof(cuDoubleComplex), cudaHostAllocWriteCombined | cudaHostAllocMapped | cudaHostAllocPortable);
  cudaHostGetDevicePointer((void **)&alpha_device, (void *)alpha, 0);
  //cudaSetDeviceFlags(cudaDeviceMapHost);
  cudaHostAlloc((void **) &beta, sizeof(cuDoubleComplex), cudaHostAllocWriteCombined | cudaHostAllocMapped | cudaHostAllocPortable); 
  cudaHostGetDevicePointer((void **)&beta_device, (void *)beta, 0);
  //*******************************************************************************************************************

  cudaStat[1]=cudaDeviceSynchronize();
  if (cudaStat[1] != cudaSuccess){ CLEANUP("Device Sync Error failed A\n");}

  //*******************************************************************************************************************


  // initalization of CUBLAS library
  cubl_status = cublasCreate(&cubl_handle);
  if(cubl_status != CUBLAS_STATUS_SUCCESS){CLEANUP("CUBLAS initilization failed \n");}

  cubl_status = cublasSetPointerMode(cubl_handle, CUBLAS_POINTER_MODE_DEVICE);
  if(cubl_status != CUBLAS_STATUS_SUCCESS){CLEANUP("CUBLAS setting device pointer mode failed\n");}

  /* initialize cusparse library */
  cusp_status= cusparseCreate(&cusp_handle);
  if (cusp_status != CUSPARSE_STATUS_SUCCESS){CLEANUP("CUSPARSE Library initialization failed\n");}

  /* create and setup matrix descriptor */
  cusp_status= cusparseCreateMatDescr(&cusp_descra);
  if (cusp_status != CUSPARSE_STATUS_SUCCESS){CLEANUP("Matrix descriptor initialization failed 1");}
  
  cusp_status=cusparseSetMatType(cusp_descra, CUSPARSE_MATRIX_TYPE_GENERAL);
  if (cusp_status != CUSPARSE_STATUS_SUCCESS){CLEANUP("Matrix descriptor initialization failed 2");}

  cusp_status=cusparseSetMatIndexBase(cusp_descra, CUSPARSE_INDEX_BASE_ONE);
  if (cusp_status != CUSPARSE_STATUS_SUCCESS){CLEANUP("Matrix descriptor initialization failed 3");}


  //cublasSetStream(cubl_handle, stream0);
  cubl_status = cublasDznrm2(cubl_handle, n_ham, eigen_seed_device, 1, norm_eigen_seed_device);
  if (cubl_status != CUBLAS_STATUS_SUCCESS){CLEANUP("CUBLAS initilization failed \n");}

  cudaStat[1] = cudaMemcpy(norm_eigen_seed, norm_eigen_seed_device,(size_t)(sizeof(double)), cudaMemcpyDeviceToHost);
  if (cudaStat[1] != cudaSuccess){ CLEANUP("copy to host failed 1\n"); }

  scl_mul_vct_kernel<<<numBlocks, threadsPerBlock>>>( eigen_seed_device, eigen_seed_device, n_ham, 1.0/(*norm_eigen_seed));

  if (nr_eigv > 1) project_out(cubl_handle, eigen_seed_device, eigen_v_device, nr_eigv-1, n_ham, 
                                                                                         numBlocks, threadsPerBlock);
       
       //"--------+-----+-----+------------------+------------------+-------------------+-------------------+"
  printf("         niters nsteps       test_E             energy              error            ort_chk\n"); 

  while(!test)
  {
    alpha[0].x = 0.0;
    beta[0].x = 0.0;
  
    //cubl_status = cublasDznrm2(cubl_handle, n_ham, eigen_seed_device, 1, norm_eigen_seed_device);
    //cudaStat[1] = cudaMemcpy(norm_eigen_seed, norm_eigen_seed_device,(size_t)(sizeof(double)), cudaMemcpyDeviceToHost);
    //if (cudaStat[1] != cudaSuccess){ CLEANUP("copy to host failed 1\n"); }

    vct_init_zero_kernel<<<numBlocks, threadsPerBlock>>>( p_previous_device, n_ham);
    vct_init_zero_kernel<<<numBlocks, threadsPerBlock>>>( p_next_device, n_ham);
    vct_init_zero_kernel<<<numBlocks, threadsPerBlock>>>( p_aux_device, n_ham);
    cpy_vct1_to_vct2_kernel<<<numBlocks, threadsPerBlock>>>( eigen_seed_device, p_current_device, n_ham);

    //cudaStat[1]=cudaDeviceSynchronize();
    //if (cudaStat[1] != cudaSuccess){CLEANUP("Device Sync Error failed A\n");}

    test_fast = 0; 

    for(j = 0; j < nstep; j++)
    {
      
      /* 
      printf("Timing spMV \n");
      time_t time_start = time(NULL);

      for(counter =0; counter <10000; counter++)
      {
        spmv_csr_vector_kernel<<<numBlocks, threadsPerBlock>>>(n_ham, rowptr_device_real, colptr_device_real, valptr_device_real, 
                                                             p_device_real, p_device_real_a, repeat, coop);
 
        //cusparseSetStream(cusp_handle, stream1);
        //cusp_status=cusparseDcsrmv(cusp_handle,CUSPARSE_OPERATION_NON_TRANSPOSE, n_ham, n_ham, size_mat_real, 
        //                           &scalar_1, cusp_descra, valptr_device_real, rowptr_device_real, 
        //                          colptr_device_real, p_device_real, &scalar_2, p_device_real_a);

      }     
      time_t time_end = time(NULL);
      printf("time elapsed  %d in sec\n", (time_end-time_start));
      exit(0);
      */      

      /* 
      spmv(rowptr_device_real, colptr_device_real, valptr_device_real, size_mat_real, 
           rowptr_device_img, colptr_device_img, valptr_device_img,    size_mat_img,
           p_current_device, p_aux_device, n_ham, cusp_handle, stream1, stream2, stream3, stream4);

      spmv(rowptr_device_real, colptr_device_real, valptr_device_real,  size_mat_real,
           rowptr_device_img, colptr_device_img, valptr_device_img,     size_mat_img,
           p_aux_device, p_aux2_device, n_ham, cusp_handle,stream1, stream2, stream3, stream4);
      */
 
      //********************************************************************************************************************

      // update | J + 1 > to  | J + 1 > + A * A | j >
      //-------------------------------------------------------------------
      spmv_csr_hybrid_kernel<<<numBlocksMul, BLOCK_SIZE_MUL, 0, stream1>>>(n_ham, rowptr_device_real, colptr_device_real, 
                                                                        valptr_device_real, 
                                                                         p_current_device, p_aux_device, repeat, coop);
 
      spmv_csr_hybrid_kernel<<<numBlocks, BLOCK_SIZE,0, stream2>>>(n_ham, rowptr_device_img, colptr_device_img, 
                                                                       valptr_device_img, 
                                                                       p_current_device,  p_aux1_device);

      cudaDeviceSynchronize();

      vct_pls_scl_mul_vct_kernel<<<numBlocks, threadsPerBlock>>>(p_aux_device, p_aux1_device, jcmpx, n_ham);
      
      spmv_csr_hybrid_kernel<<<numBlocksMul, BLOCK_SIZE_MUL, 0, stream1>>>(n_ham, rowptr_device_real, colptr_device_real, 
                                                                        valptr_device_real, 
                                                                         p_aux_device, p_aux2_device, repeat, coop);
 
      spmv_csr_hybrid_kernel<<<numBlocks, BLOCK_SIZE, 0, stream2>>>(n_ham, rowptr_device_img, colptr_device_img, 
                                                                     valptr_device_img, 
                                                                      p_aux_device,  p_aux1_device);

      cudaDeviceSynchronize();

      vct_pls_scl_mul_vct_kernel<<<numBlocks, threadsPerBlock>>>(p_aux2_device, p_aux1_device, jcmpx, n_ham);
      //********************************************************************************************************************

      cubl_status = cublasZdotc(cubl_handle, n_ham, p_aux_device, 1, p_aux_device, 1, alpha_device);

      vct_pls_vct_kernel<<<numBlocks, threadsPerBlock>>>( p_next_device, p_aux2_device, n_ham);

      cudaStat[1] = cudaMemcpy(alpha, alpha_device,(size_t)(sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);
      if (cudaStat[1] != cudaSuccess){ CLEANUP("copy to host failed 4\n");}

      vct_pls_scl_mul_vct_kernel<<<numBlocks, threadsPerBlock>>>( p_next_device, p_current_device, -alpha[0].x, n_ham);

      if (nr_eigv > 1) project_out(cubl_handle, p_next_device, eigen_v_device, nr_eigv-1, n_ham, 
                                                                                         numBlocks, threadsPerBlock);

      cubl_status = cublasZdotc(cubl_handle, n_ham, p_next_device, 1, p_next_device, 1, beta_device);

      cudaStat[1] = cudaMemcpy(beta, beta_device,(size_t)(sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);
      if (cudaStat[1] != cudaSuccess){ CLEANUP("copy to host failed 5\n");}

      beta[0].x = sqrt(beta[0].x);

      if(j == nstep-1 && nstep > 2 && nstep < max_step-10)
      {
        cubl_status = cublasZdotc(cubl_handle, n_ham, p_next_device, 1, eigen_seed_device, 1, temp1_device);

        cudaStat[1] = cudaMemcpy(temp1, temp1_device,(size_t)(sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);
        if (cudaStat[1] != cudaSuccess){ CLEANUP("copy to host failed 6\n"); }

        ort_chk = fabs(temp1[0].x/beta[0].x);

        if (ort_chk < ort_tol) nstep = nstep + 10;
      }

      ad2[j] = alpha[0].x;
      ad1[j] = beta[0].x;

      cpy_vct1_to_vct2_kernel<<<numBlocks, threadsPerBlock>>>( p_current_device, p_previous_device, n_ham);
      scl_mul_vct_kernel<<<numBlocks, threadsPerBlock>>>( p_current_device, p_next_device, n_ham, 1.0/beta[0].x);
      scl_mul_vct_kernel<<<numBlocks, threadsPerBlock>>>( p_next_device, p_previous_device, n_ham, -beta[0].x);

      cudaStat[1]=cudaDeviceSynchronize();
      if (cudaStat[1] != cudaSuccess){ CLEANUP("Device Sync Error failed D\n"); }

    } // end of lanczos loop

    //********************************************************************************************************************

    COUNTER = COUNTER + nstep;

    //diagonalize the T matrix
    //NOTE: mkl 10.2 does not define LAPACKE so using standard version
    //error_code = LAPACKE_dstevx( LAPACK_COL_MAJOR, 'V', 'I', nstep, ad2, ad1, VL, VU, IL, IU, ABSTOL, 
    //                          &M_out, eigen_val, eigen_vec, LZ, ifail );

    dstevx( &JOBZ, &RANGE, &nstep, ad2, ad1, &VL, &VU, &IL, &IU, &ABSTOL, 
            &M_out, eigen_val, eigen_vec, &LZ, work, iwork, ifail, &error_code);

    cudaStat[1] = cudaMemcpy(eigen_vec_device, eigen_vec, (size_t)(2*ldz*sizeof(double)), cudaMemcpyHostToDevice);
    if (cudaStat[1] != cudaSuccess){CLEANUP("copy to device failed \n");}
    //********************************************************************************************************************

    energy = sqrt(fabs(eigen_val[0]));

    double mean_en = 0.0;
    for(counter=0; counter < 9; counter++)
    { 
      energy_test[counter] = energy_test[counter+1];
      mean_en += energy_test[counter];
    }
    energy_test[9] = energy;
    mean_en += energy_test[9]; mean_en /= 10.0;   

    temp=0.0;
    for(counter=0; counter < 10; counter++)
      temp += pow((energy_test[counter] - mean_en),2);

    temp=temp/10;
    double dev_en = sqrt(temp);

    if(nstep == min_step){
      if( (dev_en/mean_en) < fast_tol) test_fast = 1; // true
    }

    //********************************************************************************************************************
    // SECOND LOOP: project Ritz eigenvector. 
    //********************************************************************************************************************
    
    // RESTART ITERATION FROM  eigen_seed 
    // p_current_device = eigen_seed_device
    cpy_vct1_to_vct2_kernel<<<numBlocks, threadsPerBlock,0, stream0>>>( eigen_seed_device, p_current_device, n_ham);

    // p_previous = 0; 
    // p_next = 0; 
    // p_aux = 0; 
    // eigen_seed = 0;
    vct_init_zero_kernel<<<numBlocks, threadsPerBlock>>>( p_previous_device, n_ham);
    vct_init_zero_kernel<<<numBlocks, threadsPerBlock>>>( p_next_device, n_ham);
    vct_init_zero_kernel<<<numBlocks, threadsPerBlock>>>( p_aux_device, n_ham);
    vct_init_zero_kernel<<<numBlocks, threadsPerBlock>>>( eigen_seed_device, n_ham);
    //********************************************************************************************************************

    cudaStat[1]=cudaDeviceSynchronize();
    if (cudaStat[1] != cudaSuccess){ CLEANUP("Device Sync Error failed A\n"); }

    alpha[0].x = 0.0;
    beta[0].x = 0.0;

    for( j= 0; j < nstep; j++)
    {

      //spmv(rowptr_device_real, colptr_device_real, valptr_device_real,  size_mat_real,
      //     rowptr_device_img, colptr_device_img, valptr_device_img,     size_mat_img,
      //     p_current_device, p_aux_device, n_ham,cusp_handle, stream1, stream2, stream3, stream4);

      //spmv(rowptr_device_real, colptr_device_real, valptr_device_real,  size_mat_real,
      //     rowptr_device_img, colptr_device_img, valptr_device_img,     size_mat_img,
      //     p_aux_device, p_aux2_device, n_ham,cusp_handle, stream1, stream2, stream3, stream4);
 
      //********************************************************************************************************************
      spmv_csr_hybrid_kernel<<<numBlocksMul, BLOCK_SIZE_MUL, 0, stream1>>>(n_ham, rowptr_device_real, colptr_device_real, 
                                                                        valptr_device_real, 
                                                                         p_current_device, p_aux_device, repeat, coop);
 
      spmv_csr_hybrid_kernel<<<numBlocks, BLOCK_SIZE, 0, stream2>>>(n_ham, rowptr_device_img, colptr_device_img, 
                                                                         valptr_device_img, 
                                                                         p_current_device,  p_aux1_device);

      cudaDeviceSynchronize();

      vct_pls_scl_mul_vct_kernel<<<numBlocks, threadsPerBlock>>>(p_aux_device, p_aux1_device, jcmpx, n_ham);
      
      spmv_csr_hybrid_kernel<<<numBlocksMul, BLOCK_SIZE_MUL, 0, stream1>>>(n_ham, rowptr_device_real, colptr_device_real, 
                                                                        valptr_device_real, 
                                                                         p_aux_device, p_aux2_device, repeat, coop);
 
      spmv_csr_hybrid_kernel<<<numBlocks, BLOCK_SIZE, 0, stream2>>>(n_ham, rowptr_device_img, colptr_device_img, 
                                                                          valptr_device_img, 
                                                                          p_aux_device,  p_aux1_device);

      cudaDeviceSynchronize();

      vct_pls_scl_mul_vct_kernel<<<numBlocks, threadsPerBlock>>>(p_aux2_device, p_aux1_device, jcmpx, n_ham);
      //********************************************************************************************************************

      cubl_status = cublasZdotc(cubl_handle, n_ham, p_aux_device, 1, p_aux_device, 1, alpha_device);

  
      // |j+1> = <j| A A |j>
      vct_pls_vct_kernel<<<numBlocks, threadsPerBlock>>>( p_next_device, p_aux2_device, n_ham);
   
      cubl_status = cublasZdotc(cubl_handle, n_ham, p_aux_device, 1, p_aux_device, 1, alpha_device);
      
      cudaStat[1] = cudaMemcpy(alpha, alpha_device,(size_t)(sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);
      if (cudaStat[1] != cudaSuccess){ CLEANUP("copy to host failed 7\n"); }
  
      vct_pls_scl_mul_vct_kernel<<<numBlocks, threadsPerBlock>>>( p_next_device, p_current_device, -alpha[0].x, n_ham);
  
      // project out existing eigenvectors
      if (nr_eigv > 1) project_out(cubl_handle, p_next_device, eigen_v_device, nr_eigv-1, n_ham, 
                                                                                         numBlocks, threadsPerBlock);
  
      cubl_status = cublasZdotc(cubl_handle, n_ham, p_next_device, 1, p_next_device, 1, beta_device);
  
      cudaStat[1] = cudaMemcpy(beta, beta_device,(size_t)(sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);
      if (cudaStat[1] != cudaSuccess){ CLEANUP("copy to host failed 8\n");}
  
      beta[0].x = sqrt(beta[0].x);
  
      //vct_1_add_vct_2__mul_vct_3_kernel<<<numBlocks, threadsPerBlock>>>( eigen_seed_device, eigen_vec_device, 
      //                                                                   p_current_device, n_ham, j);
      vct_pls_scl_mul_vct_kernel<<<numBlocks, threadsPerBlock>>>( eigen_seed_device,  p_current_device, 
                                                                  eigen_vec[j], n_ham);

      cpy_vct1_to_vct2_kernel<<<numBlocks, threadsPerBlock>>>( p_current_device, p_previous_device, n_ham);
      scl_mul_vct_kernel<<<numBlocks, threadsPerBlock>>>( p_current_device, p_next_device, n_ham, 1.0/beta[0].x);
      scl_mul_vct_kernel<<<numBlocks, threadsPerBlock>>>( p_next_device, p_previous_device, n_ham,-beta[0].x);
  
      //*********************************************************************************************************************
  
      cudaStat[1]=cudaDeviceSynchronize();
  
      if (cudaStat[1] != cudaSuccess){CLEANUP("Device Sync Error failed A\n");}
  
      //*********************************************************************************************************************

    } // end second loop

    // project out existing eigenvectors
    if (nr_eigv > 1) project_out(cubl_handle, eigen_seed_device, eigen_v_device, nr_eigv-1, n_ham, 
                                                                                         numBlocks, threadsPerBlock);

    // temp_1 = <u|u> 
    cubl_status = cublasDznrm2(cubl_handle, n_ham, eigen_seed_device, 1, norm_eigen_seed_device);

    cudaStat[1] = cudaMemcpy(norm_eigen_seed, norm_eigen_seed_device,(size_t)(sizeof(double)), cudaMemcpyDeviceToHost);
    if (cudaStat[1] != cudaSuccess){ CLEANUP("copy to host failed 9\n"); }

    // normalize |u>  
    scl_mul_vct_kernel<<<numBlocks, threadsPerBlock>>>( eigen_seed_device, eigen_seed_device, n_ham, 1.0/(*norm_eigen_seed));

    cudaStat[1] = cudaMemcpy(eigen_seed, eigen_seed_device, (size_t)(n_ham*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);
    if (cudaStat[1] != cudaSuccess){ CLEANUP("copy to device failed 2\n"); }
 
    //**********************************************************************************************************************
 
    // CALCULATION of <u| H |u>    
    //**********************************************************************************************************************
    // |v> = H |u> 
    //spmv(rowptr_device_real, colptr_device_real, valptr_device_real,  size_mat_real,
    //     rowptr_device_img, colptr_device_img, valptr_device_img,     size_mat_img,
    //     eigen_seed_device, p_aux_device, n_ham, cusp_handle,stream1, stream2, stream3, stream4);
    //********************************************************************************************************************
    spmv_csr_hybrid_kernel<<<numBlocksMul, BLOCK_SIZE_MUL, 0, stream1>>>(n_ham, rowptr_device_real, colptr_device_real, 
                                                                        valptr_device_real, 
                                                                         eigen_seed_device, p_aux_device, repeat, coop);
 
    spmv_csr_hybrid_kernel<<<numBlocks, BLOCK_SIZE,0, stream2>>>(n_ham, rowptr_device_img, colptr_device_img, 
                                                                       valptr_device_img, 
                                                                       eigen_seed_device, p_aux1_device);
    cudaDeviceSynchronize();

    vct_pls_scl_mul_vct_kernel<<<numBlocks, threadsPerBlock>>>(p_aux_device, p_aux1_device, jcmpx, n_ham);
    //********************************************************************************************************************

    // temp_2 = <u|v>
    cubl_status = cublasZdotc(cubl_handle, n_ham, eigen_seed_device, 1, p_aux_device, 1, temp_2_device);
    cudaStat[1] = cudaMemcpyAsync(temp_2, temp_2_device, (size_t)(sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost, stream1);
    if (cudaStat[1] != cudaSuccess){ CLEANUP("copy to host failed 10\n");}
 
    test_E = temp_2[0].x;
 
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
      test = 1; 
      printf("Eig value = %f, Error = %18.8g\n", test_E, deltaE);
      *Energy = test_E;
      *DeltaE = deltaE;
      *COUNTER_inout = COUNTER;
    }
    else 
    {
      test = 0; 
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
 
  } //while (!test)

  //cudaStat[1] = cudaMemcpy(eigen_seed, eigen_seed_device, (size_t)(n_ham*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);
 	//if (cudaStat[1] != cudaSuccess){ CLEANUP("copy to device failed 2\n");}

  printf("(CUDA) ENDS\n");

}


void spmv(int* rowptr_device_real, int* colptr_device_real, double* valptr_device_real, int size_mat_real,
          int* rowptr_device_img, int* colptr_device_img, double* valptr_device_img, int size_mat_img,         
          cuDoubleComplex* p_current_device, cuDoubleComplex* p_aux_device, int n_ham,
          cusparseHandle_t& cusp_handle, cudaStream_t& stream1, cudaStream_t& stream2, 
          cudaStream_t& stream3, cudaStream_t& stream4)
{
  
   int coop=16;   //16 
   int repeat=5;  //5  
   int threadsPerBlock = BLOCK_SIZE_MUL;
   int numBlocks=(n_ham*coop-1)/(repeat*BLOCK_SIZE_MUL)+1;  //gridSize
 
   /*
   cusparseMatDescr_t cusp_descra=0;
   cusparseCreateMatDescr(&cusp_descra);
   cusparseSetMatType(cusp_descra, CUSPARSE_MATRIX_TYPE_GENERAL);
   cusparseSetMatIndexBase(cusp_descra, CUSPARSE_INDEX_BASE_ONE);

   const double scalar_1 = 1.0;
   const double scalar_2 = 0.0;
   */
   /* 
   double * p_device_real_a, * p_device_real_b, * p_device_img_a, * p_device_img_b;
   double * p_device_real, * p_device_img;
   
   cudaMalloc((void**) &p_device_real_a, n_ham*sizeof(double));
   cudaMalloc((void**) &p_device_real_b, n_ham*sizeof(double));
   cudaMalloc((void**) &p_device_img_a, n_ham*sizeof(double));
   cudaMalloc((void**) &p_device_img_b, n_ham*sizeof(double));
   cudaMalloc((void**) &p_device_real, n_ham*sizeof(double));
   cudaMalloc((void**) &p_device_img, n_ham*sizeof(double));
   */   

   const cuDoubleComplex jcmpx=make_cuDoubleComplex(0.0,1.0);   

   cuDoubleComplex* p_aux1;
   cudaMalloc((void**) &p_aux1, n_ham*sizeof(cuDoubleComplex));

   //***************************************************************************************************************************
   //call kernel split vector into real vector and img vector
   //split_vct_kernel<<<numBlocks, threadsPerBlock>>>( p_current_device, p_device_real, p_device_img, n_ham);
   //spmv_csr_vector_kernel<<<numBlocks, threadsPerBlock, 0, stream1>>>(n_ham, rowptr_device_real, colptr_device_real, valptr_device_real, 
   //                                                       p_device_real, p_device_real_a, repeat, coop);
 
   //spmv_csr_vector_kernel<<<numBlocks, threadsPerBlock,0, stream2>>>(n_ham, rowptr_device_real, colptr_device_real, valptr_device_real, 
   //                                                       p_device_img, p_device_img_a, repeat, coop);

   //spmv_csr_vector_kernel<<<numBlocks, threadsPerBlock,0, stream3>>>(n_ham, rowptr_device_img, colptr_device_img, valptr_device_img, 
    //                                                      p_device_img, p_device_real_b, repeat, coop);

   //spmv_csr_vector_kernel<<<numBlocks, threadsPerBlock, 0, stream4>>>(n_ham, rowptr_device_img, colptr_device_img, valptr_device_img, 
   //                                                       p_device_real, p_device_img_b, repeat, coop);
   //cudaDeviceSynchronize();
   
   //***************************************************************************************************************************


   spmv_csr_hybrid_kernel<<<numBlocks, threadsPerBlock, 0, stream1>>>(n_ham, rowptr_device_real, colptr_device_real, valptr_device_real, 
                                                          p_current_device, p_aux_device, repeat, coop);
 
   spmv_csr_hybrid_kernel<<<numBlocks, threadsPerBlock,0, stream2>>>(n_ham, rowptr_device_img, colptr_device_img, valptr_device_img, 
                                                          p_current_device,  p_aux1);

   cudaDeviceSynchronize();

   //cusparseSetStream(cusp_handle, stream1);
   //cusparseDcsrmv(cusp_handle,CUSPARSE_OPERATION_NON_TRANSPOSE, n_ham, n_ham, size_mat_real, 
   //                           &scalar_1, cusp_descra, valptr_device_real, rowptr_device_real, 
   //                         colptr_device_real, p_device_real, &scalar_2, p_device_real_a);

   //***************************************************************************************************************************
   //cusparseSetStream(cusp_handle, stream2);
   //cusparseDcsrmv(cusp_handle,CUSPARSE_OPERATION_NON_TRANSPOSE, n_ham, n_ham, size_mat_real, 
   //                          &scalar_1, cusp_descra, valptr_device_real, rowptr_device_real, 
   //                        colptr_device_real, p_device_img, &scalar_2, p_device_img_a);
   //***************************************************************************************************************************
   //cusparseSetStream(cusp_handle, stream3);
   //cusparseDcsrmv(cusp_handle,CUSPARSE_OPERATION_NON_TRANSPOSE, n_ham, n_ham, size_mat_img, 
   //                        &scalar_1, cusp_descra, valptr_device_img, rowptr_device_img, 
   //                          colptr_device_img, p_device_img, &scalar_2, p_device_real_b);
   //***************************************************************************************************************************
   //cusparseSetStream(cusp_handle, stream4);
   //cusparseDcsrmv(cusp_handle,CUSPARSE_OPERATION_NON_TRANSPOSE, n_ham, n_ham, size_mat_img, 
   //                          &scalar_1, cusp_descra, valptr_device_img, rowptr_device_img, 
   //                         colptr_device_img, p_device_real, &scalar_2, p_device_img_b);
   //******************************************************************************************************************
   
   vct_pls_scl_mul_vct_kernel<<<numBlocks, threadsPerBlock>>>(p_aux_device, p_aux1, jcmpx, n_ham);
   cudaFree(p_aux1);

   //call kernel form complex vector
   //form_complex_vct_kernel<<<numBlocks, threadsPerBlock>>>( p_device_real_a, p_device_real_b, 
   //                                                         p_device_img_a, p_device_img_b, p_aux_device, n_ham);
   /*
   cudaFree(p_device_real_a);
   cudaFree(p_device_real_b);
   cudaFree(p_device_img_a);
   cudaFree(p_device_img_b);
   cudaFree(p_device_real);
   cudaFree(p_device_img);
   */
}


//********************************************************************************************************************
//********************************************************************************************************************
//********************************************************************************************************************
//********************************************************************************************************************

void fast_lanczos_ev_cusparse_(int* Size_Mat, cuDoubleComplex* valptr, int* rowptr, int* colptr, char* sparse_fmt, 
                               int* Min_step, int* Long_step, int* Max_step, cuDoubleComplex* eigen_seed, int* N_ham, 
                               double* Fast_tol, double* Long_tol, double* Ort_tol_in, int* COUNTER_inout, int* res_flag, 
                               cuDoubleComplex* eigen_v, int* NR_eigv, double* Energy, double* DeltaE) 
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
  
  const cuDoubleComplex  scalar_1 = make_cuDoubleComplex(1.0000,0.0000);
  const cuDoubleComplex  scalar_2 = make_cuDoubleComplex(0.0000,0.0000);

  cuDoubleComplex* alpha = (cuDoubleComplex *) malloc(sizeof(cuDoubleComplex)); 
  cuDoubleComplex* beta = (cuDoubleComplex *) malloc(sizeof(cuDoubleComplex)); 
  cuDoubleComplex* alpha_device, *beta_device;

  // Tolerance test variables
  double tol, errors;
  double test_E;
  double ort_tol = ort_tol_in;
  double ort_chk = 0;
  int test_fast, test;
  // Counters
  int counter, j, nstep;
  
  cuDoubleComplex* temp1 = (cuDoubleComplex *) malloc(sizeof(cuDoubleComplex));
  cuDoubleComplex* temp1_device;
  cuDoubleComplex* temp_1, * temp_2, * temp_1_device, * temp_2_device;
  curandGenerator_t gen;

  // Matrix on device
  cuDoubleComplex* valptr_device;
  int* colptr_device;
  int* rowptr_device;
  cuDoubleComplex* p_current;
  cuDoubleComplex* p_previous;
  cuDoubleComplex* p_next;
  cuDoubleComplex* p_aux;
  cuDoubleComplex* p_aux2;
  cuDoubleComplex* p_current_device;
  cuDoubleComplex* p_previous_device;
  cuDoubleComplex* p_next_device;
  cuDoubleComplex* p_aux_device;
  cuDoubleComplex* p_aux2_device;

  cuDoubleComplex* eigen_seed_device;
  cuDoubleComplex* eigen_v_device;
  double* norm_eigen_seed = (double*) malloc(sizeof(double));
  double* norm_eigen_seed_device;

  p_current = (cuDoubleComplex*) malloc(n_ham*sizeof(cuDoubleComplex));
  memset(p_current, 0, n_ham*sizeof(cuDoubleComplex));
  p_previous = (cuDoubleComplex*) malloc(n_ham*sizeof(cuDoubleComplex));
  memset(p_previous, 0, n_ham*sizeof(cuDoubleComplex));
  p_next = (cuDoubleComplex*) malloc(n_ham*sizeof(cuDoubleComplex));
  memset(p_next, 0, n_ham*sizeof(cuDoubleComplex));
  p_aux = (cuDoubleComplex*) malloc(n_ham*sizeof(cuDoubleComplex));
  memset(p_aux, 0, n_ham*sizeof(cuDoubleComplex));
  p_aux2 = (cuDoubleComplex*) malloc(n_ham*sizeof(cuDoubleComplex));
  memset(p_aux2, 0, n_ham*sizeof(cuDoubleComplex));

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

  set_gpu_id();

  int devCount = 0;	
  //cudaGetDeviceCount(&devCount);
  //printf("Device count = %d\n", devCount);
  //int list[4] = {2, 3, 0, 1};
  //cudaSetValidDevices(list, devCount);
  //scanf("%d", &devCount);
  //printf("Device chosen = %d\n", devCount);

  cudaDeviceReset();
  //cudaThreadExit();
  //cudaSetDevice(devCount);
  cudaGetDevice(&devCount);
  printf("Device used = %d\n", devCount);
  /* create a sparse and dense vector */

  temp_1 = (cuDoubleComplex *) malloc(sizeof(cuDoubleComplex));
  temp_2 = (cuDoubleComplex *) malloc(sizeof(cuDoubleComplex));

  cudaError_t cudaStat[23];
  for (counter=0; counter < 23; counter++) cudaStat[counter] = cudaSuccess;

  cudaStat[1]  = cudaMalloc((void**) &colptr_device, size_mat*sizeof(int));
  cudaStat[2]  = cudaMalloc((void**) &rowptr_device, (n_ham+1)*sizeof(int));
  cudaStat[3]  = cudaMalloc((void**) &valptr_device, size_mat*sizeof(cuDoubleComplex));
  cudaStat[4]  = cudaMalloc((void**) &p_current_device, n_ham*sizeof(cuDoubleComplex));
  cudaStat[5]  = cudaMalloc((void**) &p_previous_device, n_ham*sizeof(cuDoubleComplex));
  cudaStat[6]  = cudaMalloc((void**) &p_next_device, n_ham*sizeof(cuDoubleComplex));
  cudaStat[7]  = cudaMalloc((void**) &p_aux_device, n_ham*sizeof(cuDoubleComplex));
  cudaStat[8]  = cudaMalloc((void**) &p_aux2_device, n_ham*sizeof(cuDoubleComplex));
  cudaStat[9]  = cudaMalloc((void**) &eigen_seed_device, n_ham*sizeof(cuDoubleComplex));
  cudaStat[10] = cudaMalloc((void**) &eigen_v_device, n_ham*nr_eigv*sizeof(cuDoubleComplex));
  cudaStat[14] = cudaMalloc((void**) &eigen_vec_device, ldz*2*sizeof(double));
  cudaStat[15] = cudaMalloc((void**) &norm_eigen_seed_device, sizeof(double));
  cudaStat[16] = cudaMalloc((void**) &alpha_device, sizeof(cuDoubleComplex));
  cudaStat[17] = cudaMalloc((void**) &beta_device, sizeof(cuDoubleComplex));
  cudaStat[18] = cudaMalloc((void**) &temp1_device, sizeof(cuDoubleComplex));
  cudaStat[20] = cudaMalloc((void**) &temp_1_device, sizeof(cuDoubleComplex));
  cudaStat[21] = cudaMalloc((void**) &temp_2_device, sizeof(cuDoubleComplex));

  for (counter=1; counter < 23; counter++)
    if (cudaStat[counter] != cudaSuccess){CLEANUP("Device malloc failed\n")};

  test = 0; // FALSE
  nstep = min_step;
  for(counter=0; counter < n_ham; counter++)
  {
    eigen_seed[counter].x = rand()%100/100.0;
    eigen_seed[counter].y = rand()%100/100.0;
  }                        
    
  cudaStat[1] = cudaMemcpy(colptr_device, colptr, (size_t)(size_mat*sizeof(int)), cudaMemcpyHostToDevice);
  cudaStat[2] = cudaMemcpy(rowptr_device, rowptr, (size_t)((n_ham+1)*sizeof(int)), cudaMemcpyHostToDevice);
  cudaStat[3] = cudaMemcpy(valptr_device, valptr, (size_t)(size_mat*sizeof(cuDoubleComplex)), cudaMemcpyHostToDevice);
  cudaStat[4] = cudaMemcpy(eigen_seed_device, eigen_seed, (size_t)(n_ham*sizeof(cuDoubleComplex)), cudaMemcpyHostToDevice);
  cudaStat[5] = cudaMemcpy(eigen_v_device, eigen_v, (size_t)(n_ham*nr_eigv*sizeof(cuDoubleComplex)), cudaMemcpyHostToDevice);

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

  // norm = sqrt(< seed | seed >)
  cubl_status = cublasDznrm2(cubl_handle, n_ham, eigen_seed_device, 1, norm_eigen_seed_device);

  if(cubl_status != CUBLAS_STATUS_SUCCESS){ CLEANUP("CUBLAS initilization failed \n");}

  // create the sparse matrix in CSR format
  cusparseSpMatDescr_t matH;
  cusp_status = cusparseCreateCsr(&matH, n_ham, n_ham, size_mat,
                                  rowptr_device, colptr_device, valptr_device,
                                  CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
				  CUSPARSE_INDEX_BASE_ONE, CUDA_C_64F);
  if (cusp_status != CUSPARSE_STATUS_SUCCESS) { CLEANUP("Failed to create CSR matrix on device\n");}


  // create the dense vectors
  cusparseDnVecDescr_t vec_current, vec_previous, vec_next, vec_aux, vec_aux2, vec_eigen_seed;
  cusp_status = cusparseCreateDnVec(&vec_current, n_ham, p_current_device, CUDA_C_64F);
  if (cusp_status != CUSPARSE_STATUS_SUCCESS) { CLEANUP("Failed to create dense vector on device\n");}

  cusp_status = cusparseCreateDnVec(&vec_previous, n_ham, p_previous_device, CUDA_C_64F);
  if (cusp_status != CUSPARSE_STATUS_SUCCESS) { CLEANUP("Failed to create dense vector on device\n");}

  cusp_status = cusparseCreateDnVec(&vec_next, n_ham, p_next_device, CUDA_C_64F);
  if (cusp_status != CUSPARSE_STATUS_SUCCESS) { CLEANUP("Failed to create dense vector on device\n");}

  cusp_status = cusparseCreateDnVec(&vec_aux, n_ham, p_aux_device, CUDA_C_64F);
  if (cusp_status != CUSPARSE_STATUS_SUCCESS) { CLEANUP("Failed to create dense vector on device\n");}

  cusp_status = cusparseCreateDnVec(&vec_aux, n_ham, p_aux2_device, CUDA_C_64F);
  if (cusp_status != CUSPARSE_STATUS_SUCCESS) { CLEANUP("Failed to create dense vector on device\n");}

  cusp_status = cusparseCreateDnVec(&vec_eigen_seed, n_ham, eigen_seed_device, CUDA_C_64F);
  if (cusp_status != CUSPARSE_STATUS_SUCCESS) { CLEANUP("Failed to create dense vector on device\n");}


  cudaStat[1] = cudaMemcpy(norm_eigen_seed, norm_eigen_seed_device,(size_t)(sizeof(double)), cudaMemcpyDeviceToHost);
  if (cudaStat[1] != cudaSuccess){CLEANUP("copy to host failed 1\n"); }

  // |seed> = |seed>/norm
  scl_mul_vct_kernel<<<numBlocks, threadsPerBlock>>>(eigen_seed_device, eigen_seed_device, n_ham, 1.0/(*norm_eigen_seed));

  if (nr_eigv > 1) project_out(cubl_handle, eigen_seed_device, eigen_v_device, nr_eigv-1, n_ham, 
                                                                                         numBlocks, threadsPerBlock);

  // LOOP ////////////////////////////////////////////////////////////////////////////////////////////////////////////
       //"--------+-----+-----+------------------+------------------+-------------------+-------------------+"
  printf("         niters nsteps       test_E             energy              error            ort_chk\n"); 


  while(!test)
  {
 
    alpha[0].x = 0.0;
    beta[0].x = 0.0;

    vct_init_zero_kernel<<<numBlocks, threadsPerBlock>>>( p_previous_device, n_ham);
    vct_init_zero_kernel<<<numBlocks, threadsPerBlock>>>( p_next_device, n_ham);
    vct_init_zero_kernel<<<numBlocks, threadsPerBlock>>>( p_aux_device, n_ham);
    cpy_vct1_to_vct2_kernel<<<numBlocks, threadsPerBlock>>>( eigen_seed_device, p_current_device, n_ham);

    test_fast = 0; 

    // Lanczos Loop	
    for(j = 0; j < nstep; j++)
    {

      // potentially we need buffers for matrix operations
      void* dBuffer     = NULL;
      size_t bufferSize = 0;

      // update | J + 1 > to  | J + 1 > + A * A | j >
      //-------------------------------------------------------------------
      cusp_status = cusparseSpMV_bufferSize(cusp_handle,CUSPARSE_OPERATION_NON_TRANSPOSE, 
                                    &scalar_1, matH, vec_current,
				    &scalar_2, vec_aux, CUDA_C_64F,
				    CUSPARSE_MV_ALG_DEFAULT, &bufferSize);
      if (cusp_status != CUSPARSE_STATUS_SUCCESS){ CLEANUP("Failed to figure out needed buffer size\n"); }

      cudaStat[1] = cudaMalloc(&dBuffer, bufferSize);
      if (cudaStat[1] != cudaSuccess){CLEANUP("Failed to alloc buffer\n"); }

      cusp_status = cusparseSpMV(cusp_handle,CUSPARSE_OPERATION_NON_TRANSPOSE, 
                                    &scalar_1, matH, vec_current,
				    &scalar_2, vec_aux, CUDA_C_64F,
				    CUSPARSE_MV_ALG_DEFAULT, dBuffer);

      if (cusp_status != CUSPARSE_STATUS_SUCCESS){ CLEANUP("Matrixvector multiplication failed 1 \n"); }

      cusp_status = cusparseSpMV(cusp_handle,CUSPARSE_OPERATION_NON_TRANSPOSE, 
                                    &scalar_1, matH, vec_aux,
				    &scalar_2, vec_aux2, CUDA_C_64F,
				    CUSPARSE_MV_ALG_DEFAULT, dBuffer);

      if (cusp_status != CUSPARSE_STATUS_SUCCESS){ CLEANUP("Matrixvector multiplication failed 2"); }


      vct_pls_vct_kernel<<<numBlocks, threadsPerBlock>>>( p_next_device, p_aux2_device, n_ham);
 
      // alpha = < j | A * A | j >
      //-------------------------------------------------------------------
      cubl_status = cublasZdotc(cubl_handle, n_ham, p_aux_device, 1, p_aux_device, 1, alpha_device);
      if(cubl_status != CUBLAS_STATUS_SUCCESS){ CLEANUP("CUBLAS dot product failed hello1  %d\n"); }

      cudaStat[1] = cudaMemcpy(alpha, alpha_device,(size_t)(sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);
      if (cudaStat[1] != cudaSuccess){ CLEANUP("copy to host failed 4\n"); }

      //  update | J + 1 > to | J + 1 > - alpha | j >
      //-------------------------------------------------------------------    
      vct_pls_scl_mul_vct_kernel<<<numBlocks, threadsPerBlock>>>( p_next_device, p_current_device, -alpha[0].x, n_ham);
 
      //  project out existing eigenvectors
      //-------------------------------------------------------------------    
      if (nr_eigv > 1) project_out(cubl_handle, p_next_device, eigen_v_device, nr_eigv-1, n_ham, 
                                                                                         numBlocks, threadsPerBlock);
      // beta = SQRT( < J + 1 | J + 1 > )
      //------------------------------------------------------------------- 
      cubl_status = cublasZdotc(cubl_handle, n_ham, p_next_device, 1, p_next_device, 1, beta_device);
      if(cubl_status != CUBLAS_STATUS_SUCCESS){ CLEANUP("CUBLAS dot product failed hello2\n");}

      cudaStat[1] = cudaMemcpy(beta, beta_device,(size_t)(sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);
      if (cudaStat[1] != cudaSuccess){ CLEANUP("copy to host failed 5\n"); }

      beta[0].x = sqrt(beta[0].x);
      //------------------------------------------------------------------- 
      // orthogonality check
      if (j == nstep-1 && nstep > 2 && nstep < max_step-10)
      {
        cubl_status = cublasZdotc(cubl_handle, n_ham, p_next_device, 1, eigen_seed_device, 1, temp1_device);
        if(cubl_status != CUBLAS_STATUS_SUCCESS){CLEANUP("CUBLAS dot product failed 3\n");}

        cudaStat[1] = cudaMemcpy(temp1, temp1_device,(size_t)(sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);
        if (cudaStat[1] != cudaSuccess){ CLEANUP("copy to host failed 6\n");}

        ort_chk = fabs(temp1[0].x/beta[0].x);

        if(ort_chk < ort_tol)  nstep = nstep + 10;
      }
      //------------------------------------------------------------------- 
      ad2[j] = alpha[0].x;
      ad1[j] = beta[0].x;
      //------------------------------------------------------------------- 
      // update |j - 1> to | j >
      cpy_vct1_to_vct2_kernel<<<numBlocks, threadsPerBlock>>>( p_current_device, p_previous_device, n_ham);
      // update |j> to |j + 1>/beta
      scl_mul_vct_kernel<<<numBlocks, threadsPerBlock>>>( p_current_device, p_next_device, n_ham, 1.0/beta[0].x);
      // initialize |j + 1> = -beta |j> (= -beta |j-1>)
      scl_mul_vct_kernel<<<numBlocks, threadsPerBlock>>>( p_next_device, p_previous_device, n_ham, -beta[0].x);

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

    //cudaStat[1] = cudaMemcpy(eigen_vec_device, eigen_vec, (size_t)(2*ldz*sizeof(double)), cudaMemcpyHostToDevice);
    //if (cudaStat[1] != cudaSuccess){ CLEANUP("copy to device failed \n"); }

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

    alpha[0].x = 0.0;
    beta[0].x = 0.0;
    
    for( j= 0; j < nstep; j++)
    { 
      void* dBuffer     = NULL;
      size_t bufferSize = 0;

      cusp_status = cusparseSpMV_bufferSize(cusp_handle,CUSPARSE_OPERATION_NON_TRANSPOSE, 
                                    &scalar_1, matH, vec_current,
				    &scalar_2, vec_aux, CUDA_C_64F,
				    CUSPARSE_MV_ALG_DEFAULT, &bufferSize);
      if (cusp_status != CUSPARSE_STATUS_SUCCESS){ CLEANUP("Failed to figure out needed buffer size\n"); }

      cudaStat[1] = cudaMalloc(&dBuffer, bufferSize);
      if (cudaStat[1] != cudaSuccess){CLEANUP("Failed to alloc buffer\n"); }

      // |aux> = A |j> 
      cusp_status = cusparseSpMV(cusp_handle,CUSPARSE_OPERATION_NON_TRANSPOSE, 
                                    &scalar_1, matH, vec_current,
				    &scalar_2, vec_aux, CUDA_C_64F,
				    CUSPARSE_MV_ALG_DEFAULT, dBuffer);

      if (cusp_status != CUSPARSE_STATUS_SUCCESS){ CLEANUP("Matrixvector multiplication failed 1 \n"); }

      // |aux2> = A |aux> 
      cusp_status = cusparseSpMV(cusp_handle,CUSPARSE_OPERATION_NON_TRANSPOSE, 
                                    &scalar_1, matH, vec_aux,
				    &scalar_2, vec_aux2, CUDA_C_64F,
				    CUSPARSE_MV_ALG_DEFAULT, dBuffer);

      if (cusp_status != CUSPARSE_STATUS_SUCCESS){ CLEANUP("Matrixvector multiplication failed 2"); }


      // |j+1> = |j+1> + |aux2> = |j+1> + A^2 |j> 
      vct_pls_vct_kernel<<<numBlocks, threadsPerBlock>>>( p_next_device, p_aux2_device, n_ham);

      //  alpha = <j| A A |j> = <j | j+1>  
      cubl_status = cublasZdotc(cubl_handle, n_ham, p_aux_device, 1, p_aux_device, 1, alpha_device);
      if(cubl_status != CUBLAS_STATUS_SUCCESS){ CLEANUP("CUBLAS dot product failed hello4\n"); }

      cudaStat[1] = cudaMemcpy(alpha, alpha_device,(size_t)(sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);
      if (cudaStat[1] != cudaSuccess){ CLEANUP("copy to host failed 7\n");}

      // |j+1> = |j+1> - alpha |j>
      vct_pls_scl_mul_vct_kernel<<<numBlocks, threadsPerBlock>>>( p_next_device, p_current_device, -alpha[0].x, n_ham);

      if (nr_eigv > 1) project_out(cubl_handle, p_next_device, eigen_v_device, nr_eigv-1, n_ham, 
                                                                                         numBlocks, threadsPerBlock);

      // beta = <j+1|j+1> 
      cubl_status = cublasZdotc(cubl_handle, n_ham, p_next_device, 1, p_next_device, 1, beta_device);
      if (cubl_status != CUBLAS_STATUS_SUCCESS){ CLEANUP("CUBLAS dot product failed hello5\n");}

      cudaStat[1] = cudaMemcpy(beta, beta_device,(size_t)(sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);
      if (cudaStat[1] != cudaSuccess){ CLEANUP("copy to host failed 8\n"); }

      beta[0].x = sqrt(beta[0].x);

      // Ritz Projection: |u> = |u> + Vj*|j> 
      //vct_1_add_vct_2__mul_vct_3_kernel<<<numBlocks, threadsPerBlock>>>( eigen_seed_device, eigen_vec_device, 
      //                                                                    p_current_device, n_ham, j);

      vct_pls_scl_mul_vct_kernel<<<numBlocks, threadsPerBlock>>>( eigen_seed_device,  p_current_device, 
                                                                  eigen_vec[j], n_ham);
      // |j-1>  <- |j> 
      cpy_vct1_to_vct2_kernel<<<numBlocks, threadsPerBlock>>>( p_current_device, p_previous_device, n_ham);

      // |j>  <- |j+1>/beta 
      scl_mul_vct_kernel<<<numBlocks, threadsPerBlock>>>( p_current_device, p_next_device, n_ham, 1.0/beta[0].x);

      // Initialize |j+1> = - beta |j>  ( = -beta |j-1> )
      scl_mul_vct_kernel<<<numBlocks, threadsPerBlock>>>( p_next_device, p_previous_device, n_ham, -beta[0].x);

    } // end second loop


    if (nr_eigv > 1) project_out(cubl_handle, eigen_seed_device, eigen_v_device, nr_eigv-1, n_ham, 
                                                                                         numBlocks, threadsPerBlock);
    // temp_1 = <u|u> 
    cubl_status = cublasDznrm2(cubl_handle, n_ham, eigen_seed_device, 1, norm_eigen_seed_device);

    cudaStat[1] = cudaMemcpy(norm_eigen_seed, norm_eigen_seed_device,(size_t)(sizeof(double)), cudaMemcpyDeviceToHost);
    if (cudaStat[1] != cudaSuccess){ CLEANUP("copy to host failed 9\n"); }

    // normalize |u>  
    scl_mul_vct_kernel<<<numBlocks, threadsPerBlock>>>(eigen_seed_device, eigen_seed_device, n_ham, 1.0/(*norm_eigen_seed));

    cudaStat[1] = cudaMemcpy(eigen_seed, eigen_seed_device, (size_t)(n_ham*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);
    if (cudaStat[1] != cudaSuccess){ CLEANUP("copy to device failed 2\n");}

    //************************************************************************************************************************
 
    // CALCULATION of <u| H |u>    
    //************************************************************************************************************************
    // |v> = H |u> 
    void* dBuffer     = NULL;
    size_t bufferSize = 0;

    cusp_status = cusparseSpMV_bufferSize(cusp_handle,CUSPARSE_OPERATION_NON_TRANSPOSE, 
                                    &scalar_1, matH, vec_eigen_seed,
				    &scalar_2, vec_aux, CUDA_C_64F,
				    CUSPARSE_MV_ALG_DEFAULT, &bufferSize);
    if (cusp_status != CUSPARSE_STATUS_SUCCESS){ CLEANUP("Failed to figure out needed buffer size\n"); }

    cudaStat[1] = cudaMalloc(&dBuffer, bufferSize);
    if (cudaStat[1] != cudaSuccess){CLEANUP("Failed to alloc buffer\n"); }

    cusp_status = cusparseSpMV(cusp_handle,CUSPARSE_OPERATION_NON_TRANSPOSE, 
                                    &scalar_1, matH, vec_eigen_seed,
				    &scalar_2, vec_aux, CUDA_C_64F,
				    CUSPARSE_MV_ALG_DEFAULT, dBuffer);

    // temp_2 = <u|v>
    cubl_status = cublasZdotc(cubl_handle, n_ham, eigen_seed_device, 1, p_aux_device, 1, temp_2_device);
    cudaStat[1] = cudaMemcpy(temp_2, temp_2_device,(size_t)(sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);
    if (cudaStat[1] != cudaSuccess){ CLEANUP("copy to host failed 10\n");}

    //printf("temp_1=%f, temp_2=%f\n",temp_1[0].x,temp_2[0].x);
    test_E = temp_2[0].x;

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
void project_out(cublasHandle_t cubl_handle, cuDoubleComplex* vec, cuDoubleComplex* eigen_v_device, int nr_eigv, 
                 int n_ham, int numBlocks, int threadsPerBlock)
{
  
  cublasStatus_t cubl_status;
  
  //cuDoubleComplex* val = (cuDoubleComplex *) malloc(sizeof(cuDoubleComplex));
 
  cuDoubleComplex* alpha;  
  cudaMalloc((void**) &alpha, sizeof(cuDoubleComplex));

  for (unsigned int k=0; k < nr_eigv; k++)
  {
     //cudaMemcpy(val, &eigen_v_device[k*n_ham], (size_t)(sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);
     //printf("eigv: %g %g \n",val[0].x,val[0].y); 

     cubl_status = cublasZdotc(cubl_handle, n_ham,  &eigen_v_device[k*n_ham], 1, vec, 1, alpha);
     project_out_kernel<<<numBlocks, threadsPerBlock>>>(vec, &eigen_v_device[k*n_ham], alpha, n_ham);
  }

  cudaFree(alpha);
}


/*
__global__ void vct_1_add_vct_2__mul_vct_3_kernel(cuDoubleComplex* vct_1, double* vct_2, cuDoubleComplex* vct_3, int n_ham, int j)
{
  int ID;
  ID = blockDim.x * blockIdx.x + threadIdx.x;
  if(ID < n_ham)
  {
    vct_1[ID].x = vct_1[ID].x + vct_2[j] * vct_3[ID].x;
    vct_1[ID].y = vct_1[ID].y + vct_2[j] * vct_3[ID].y;
  }
}
*/

__global__ void scl_mul_vct_kernel(cuDoubleComplex* v1, cuDoubleComplex* v2, int vct_size, double slr)
{
  int ID;
  //cuDoubleComplex * dptr;
  ID = blockDim.x * blockIdx.x + threadIdx.x;
  //dptr = (cuDoubleComplex *)((char *)V+j*pitch);
  //for(m=0; m < Vct_Size; m++){
  if(ID < vct_size)
  {
    v1[ID].x = slr * v2[ID].x;
    v1[ID].y = slr * v2[ID].y;
  }
}



__global__ void scl_mul_vct_kernel(cuDoubleComplex* v1, cuDoubleComplex* v2, int vct_size, cuDoubleComplex slr)
{
  int ID;
  //cuDoubleComplex * dptr;
  ID = blockDim.x * blockIdx.x + threadIdx.x;
  //dptr = (cuDoubleComplex *)((char *)V+j*pitch);
  //for(m=0; m < Vct_Size; m++){
  if(ID < vct_size)
  {
    v1[ID].x = slr.x * v2[ID].x - slr.y * v2[ID].y;
    v1[ID].y = slr.x * v2[ID].y + slr.y * v2[ID].x;
  }
}



__global__ void vct_init_zero_kernel(cuDoubleComplex* vct, int n_ham)
{
  int ID;
  ID = blockDim.x * blockIdx.x + threadIdx.x;
  if(ID < n_ham)
  {
    vct[ID].x = 0.00000000000000;
    vct[ID].y = 0.00000000000000;
  }
}


__global__ void cpy_vct1_to_vct2_kernel(cuDoubleComplex* vct_scr, cuDoubleComplex* vct_des, int n_ham)
{
  int ID;
  ID = blockDim.x * blockIdx.x + threadIdx.x;
  if(ID < n_ham)
  {
    vct_des[ID].x = vct_scr[ID].x;
    vct_des[ID].y = vct_scr[ID].y;
  }
}


__global__ void vct_pls_vct_kernel(cuDoubleComplex* vct_1, cuDoubleComplex* vct_2, int n_ham)
{
  int ID;
  ID = blockDim.x * blockIdx.x + threadIdx.x;
  if(ID < n_ham)
  {
    vct_1[ID].x = vct_1[ID].x + vct_2[ID].x;
    vct_1[ID].y = vct_1[ID].y + vct_2[ID].y;
  }
}



__global__ void vct_pls_scl_mul_vct_kernel(cuDoubleComplex* vct1, cuDoubleComplex* vct2, double scalar, int n_ham)
{
  int ID;
  ID = blockDim.x * blockIdx.x + threadIdx.x;
  if(ID < n_ham)
  {
    vct1[ID].x =  vct1[ID].x + (scalar * vct2[ID].x);
    vct1[ID].y =  vct1[ID].y + (scalar * vct2[ID].y);
  }
}


__global__ void vct_pls_scl_mul_vct_kernel(cuDoubleComplex* vct1, cuDoubleComplex* vct2, cuDoubleComplex scalar, int n_ham)
{
  int ID;
  ID = blockDim.x * blockIdx.x + threadIdx.x;
  if(ID < n_ham)
  {
    vct1[ID].x =  vct1[ID].x + scalar.x * vct2[ID].x - scalar.y * vct2[ID].y;
    vct1[ID].y =  vct1[ID].y + scalar.x * vct2[ID].y + scalar.y * vct2[ID].x;
  }
}



__global__ void project_out_kernel(cuDoubleComplex* vec, cuDoubleComplex* eigen_vector, cuDoubleComplex* alpha, int n_ham)
{

  int  ID;
  ID = blockDim.x * blockIdx.x + threadIdx.x;

  if (ID < n_ham)
  {
      vec[ID].x = vec[ID].x - alpha[0].x * eigen_vector[ID].x + alpha[0].y * eigen_vector[ID].y;
      vec[ID].y = vec[ID].y - alpha[0].x * eigen_vector[ID].y - alpha[0].y * eigen_vector[ID].x;
  }
}





__global__ void split_vct_kernel( cuDoubleComplex * scr_vct, double * dest_real, double * dest_img, int n_ham)
{
  int ID;
  ID = blockDim.x * blockIdx.x + threadIdx.x;
  if(ID < n_ham)
  {
    dest_real[ID] = scr_vct[ID].x;
    dest_img[ID]  = scr_vct[ID].y;
  }
}


  __global__ void form_complex_vct_kernel( double * real_a, double * real_b, double * img_a, double * img_b, 
                                           cuDoubleComplex * result, int n_ham)
  {
    int ID;
    ID = blockDim.x * blockIdx.x + threadIdx.x;
    if(ID < n_ham)
    {
      result[ID].x = real_a[ID] - real_b[ID];
      result[ID].y = img_a[ID] + img_b[ID];
    }
  }




  //__global__ double dot_product_kernel(cuDoubleComplex* vec1, cuDoubleComplex* vec2)
  // {
  //  int ID;
  //  ID = blockDim.x * blockIdx.x + threadIdx.x;
  //      
  //  __shared__ double cache[BLOCK_SIZE_MUL];
  //
  // 
  //}



__global__ void spmv_csr_vector_kernel(int dimRow, const int* rowPtrs, const int* colIdxs, const double* values, 
                                        const double* x, double* y, int repeat, int coop)
{
    // blockDim = gridSize
    
    int i = (repeat*blockIdx.x * blockDim.x + threadIdx.x)/coop;
    int coopIdx = threadIdx.x%coop;
    int tid = threadIdx.x;
    __shared__  double sdata[ BLOCK_SIZE_MUL ];

    for (int r = 0; r<repeat; r++) 
    {
      double localSum = 0.0;
      if (i<dimRow) 
      {
         // do multiplication
         int rowPtr = rowPtrs[i]-1;
         for (int j = coopIdx; j<(rowPtrs[i+1]-1-rowPtr); j+=coop) 
             localSum += values[rowPtr+j] * x[colIdxs[rowPtr+j]-1];
         
         // do reduction in shared mem
         sdata[tid] = localSum;
         for (unsigned int s=coop/2; s>0; s>>=1) 
         {
           if (coopIdx < s) sdata[tid] += sdata[tid + s];
         }
         if (coopIdx == 0) y[i] = sdata[tid];
         i += blockDim.x/coop;
      }
    }
}

__global__ void spmv_csr_hybrid_kernel(int dimRow, const int* rowPtrs, const int* colIdxs, const double* values, 
                                       const cuDoubleComplex* x, cuDoubleComplex* y, int repeat, int coop)
{
    // blockDim = gridSize
    
    int i = (repeat*blockIdx.x * blockDim.x + threadIdx.x)/coop;
    int coopIdx = threadIdx.x%coop;
    int tid = threadIdx.x;
    __shared__  cuDoubleComplex sdata[ BLOCK_SIZE_MUL ];

    for (int r = 0; r<repeat; r++) 
    {
      cuDoubleComplex localSum;
      localSum.x = 0.0; localSum.y = 0.0;
      if (i<dimRow) 
      {
         // do multiplication
         int rowPtr = rowPtrs[i]-1;
         int stop = rowPtrs[i+1]-1-rowPtr;
         for (int j = coopIdx; j < stop; j+=coop)
         {
             localSum.x += values[rowPtr+j] * x[colIdxs[rowPtr+j]-1].x;
             localSum.y += values[rowPtr+j] * x[colIdxs[rowPtr+j]-1].y;
         }

         // do reduction in shared mem
         sdata[tid] = localSum;
         for (unsigned int s=coop/2; s>0; s>>=1) 
         {
           if (coopIdx < s){ sdata[tid].x += sdata[tid + s].x; sdata[tid].y += sdata[tid + s].y;}
         }
         if (coopIdx == 0) y[i] = sdata[tid];
         i += blockDim.x/coop;
      }
    }

}


__global__ void spmv_csr_simple_kernel(int num_rows, const int* rowPtrs, const int* colIdxs, const double* values, 
                                       const double* x, double* y)
{
   int i = blockIdx.x * blockDim.x + threadIdx.x;
   if (i<num_rows)
   {
     double rowSum = 0.0;
     int row_start = rowPtrs[i]-1;
     int row_end = rowPtrs[i+1]-1;
     for (int j=row_start; j<row_end; j++)
     {
        rowSum += values[j] * x[colIdxs[j]-1];
     }
     y[i] = rowSum;
   }

} 

__global__ void spmv_csr_hybrid_kernel(int num_rows, const int* rowPtrs, const int* colIdxs, const double* values, 
                                       const cuDoubleComplex* x, cuDoubleComplex* y)
{
   int i = blockIdx.x * blockDim.x + threadIdx.x;
   if (i<num_rows)
   {
     cuDoubleComplex rowSum;
     rowSum.x = 0.0; rowSum.y = 0.0;
     int row_start = rowPtrs[i]-1;
     int row_end = rowPtrs[i+1]-1;
     for (int j=row_start; j<row_end; j++)
     {
        rowSum.x += values[j] * x[colIdxs[j]-1].x;
        rowSum.y += values[j] * x[colIdxs[j]-1].y;
     }
     y[i] = rowSum;
   }
} 


  // __global__ void spmv_csr_vector_kernel(int num_rows, int * ptr, int * indices, double * data, double * x, double * y){
  //         __shared__ double vals[(VECTORS_PER_BLOCK+1) * THREADS_PER_VECTOR];
  // 	int jj, row_start, row_end;        
  // 	int thread_id = blockDim.x * blockIdx.x + threadIdx.x ;                     // global thread index
  //         int warp_id   = thread_id / 32;                                             // global warp index
  //         int lane      = thread_id & (32 - 1);                                       // thread index within the warp
  //         // one warp per row
  //         int row = warp_id ;
  //         //y[ row ] = 0;
  //         if ( row < num_rows ){
  //                 row_start = ptr[ row ]-1 ;
  //                 row_end   = ptr[ row+1] -1;
  //                 // compute running sum per thread
  //                 vals[ threadIdx.x ] = 0.0f;
  //                 for ( jj = row_start+lane  ; jj < row_end ; jj += 32)
  //                         vals[ threadIdx.x ] += data[ jj ] * x[indices[jj]-1];
  //                 //    pa r a l l e l r e d u c t i o n in shared memory
  //                 if    ( lane < 16) vals[ threadIdx.x ] += vals[ threadIdx.x + 16];
  //                 if    ( lane < 8) vals[ threadIdx.x ] += vals[ threadIdx.x  + 8];
  //                 if    ( lane < 4) vals[ threadIdx.x ] += vals[ threadIdx.x  + 4];
  //                 if    ( lane < 2) vals[ threadIdx.x ] += vals[ threadIdx.x  + 2];
  //                 if    ( lane < 1) vals[ threadIdx.x ] += vals[ threadIdx.x  + 1];
  //                 // first thread writes the result
  //                 if ( lane == 0) y[ row ] = vals[ threadIdx.x ];
  //         }
  // }




  // __global__ void
  // spmv_csr_vector_kernel( int num_rows,
  //                         int * Ap, 
  //                         int * Aj, 
  //                         double * Ax, 
  //                         double * x, 
  //                         double * y){
  //     __shared__  volatile double sdata[BLOCK_SIZE + 16];                          // padded to avoid reduction ifs
  //     __shared__  volatile int ptrs[BLOCK_SIZE/WARP_SIZE][2];
  //     int thread_id   = BLOCK_SIZE * blockIdx.x + threadIdx.x;  // global thread index
  //     int thread_lane = threadIdx.x & (WARP_SIZE-1);            // thread index within the warp
  //     int warp_id     = thread_id   / WARP_SIZE;                // global warp index
  //     int warp_lane   = threadIdx.x / WARP_SIZE;                // warp index within the CTA
  //     int num_warps   = (BLOCK_SIZE / WARP_SIZE) * gridDim.x;   // total number of active warps
  // 
  //     for(int row = warp_id; row < num_rows; row += num_warps){
  //         // use two threads to fetch Ap[row] and Ap[row+1]
  //         // this is considerably faster than the straightforward version
  //         if(thread_lane < 2)
  //             ptrs[warp_lane][thread_lane] = Ap[row + thread_lane];
  //          int row_start = ptrs[warp_lane][0]-1;                   //same as: row_start = Ap[row];
  //          int row_end   = ptrs[warp_lane][1]-1;                   //same as: row_end   = Ap[row+1];
  // 
  //         // compute local sum
  //         double sum = 0;
  //         for(int jj = row_start + thread_lane; jj < row_end; jj += WARP_SIZE)
  //             sum += Ax[jj] * x[ Aj[jj]];    //fetch_x<UseCache>(Aj[jj], x);
  // 
  //         // reduce local sums to row sum (ASSUME: warpsize 32)
  //         sdata[threadIdx.x] = sum;
  //         sdata[threadIdx.x] = sum = sum + sdata[threadIdx.x + 16]; //EMUSYNC; 
  //         sdata[threadIdx.x] = sum = sum + sdata[threadIdx.x +  8]; //EMUSYNC;
  //         sdata[threadIdx.x] = sum = sum + sdata[threadIdx.x +  4]; //EMUSYNC;
  //         sdata[threadIdx.x] = sum = sum + sdata[threadIdx.x +  2]; //EMUSYNC;
  //         sdata[threadIdx.x] = sum = sum + sdata[threadIdx.x +  1]; //EMUSYNC;
  //        
  // //// Alternative method (slightly slower)
  // //        // compute local sum
  // //        sdata[threadIdx.x] = 0;
  // //        for(int jj = row_start + thread_lane; jj < row_end; jj += WARP_SIZE)
  // //            sdata[threadIdx.x] += Ax[jj] * fetch_x<UseCache>(Aj[jj], x);
  // //
  // //        // reduce local sums to row sum (ASSUME: warpsize 32)
  // //        sdata[threadIdx.x] += sdata[threadIdx.x + 16]; EMUSYNC;
  // //        sdata[threadIdx.x] += sdata[threadIdx.x +  8]; EMUSYNC;
  // //        sdata[threadIdx.x] += sdata[threadIdx.x +  4]; EMUSYNC;
  // //        sdata[threadIdx.x] += sdata[threadIdx.x +  2]; EMUSYNC;
  // //        sdata[threadIdx.x] += sdata[threadIdx.x +  1]; EMUSYNC;
  // 
  //         // first thread writes warp result
  //         if (thread_lane == 0)
  //             y[row] += sdata[threadIdx.x];
  //     }
  // }







  //Unrolled Nvidia CSR kernel from the latest version of CUSP 
  //  __global__ void spmv_csr_vector_kernel(int num_rows,
  // 			int * Ap, 
  // 			int * Aj, 
  // 			double * Ax, 
  // 			double * x, 
  // 			double * y)
  //  {
  //      __shared__  double sdata[(VECTORS_PER_BLOCK+1) * THREADS_PER_VECTOR];      // padded to avoid reduction conditionals
  //      __shared__  int ptrs[VECTORS_PER_BLOCK][2];
  // 
  //      int THREADS_PER_BLOCK = VECTORS_PER_BLOCK * THREADS_PER_VECTOR;
  // 
  //      int thread_id   = THREADS_PER_BLOCK * blockIdx.x + threadIdx.x;    // global thread index
  //      int thread_lane = threadIdx.x & (THREADS_PER_VECTOR - 1);          // thread index within the vector
  //      int vector_id   = thread_id   /  THREADS_PER_VECTOR;               // global vector index
  //      int vector_lane = threadIdx.x /  THREADS_PER_VECTOR;               // vector index within the block
  //      int num_vectors = VECTORS_PER_BLOCK * gridDim.x;                   // total number of active vectors
  // 
  //      for(int row = vector_id; row < num_rows; row += num_vectors)
  //      {
  // 	 use two threads to fetch Ap[row] and Ap[row+1]
  // 	 this is considerably faster than the straightforward version
  // 	 if(thread_lane < 2)
  // 	     ptrs[vector_lane][thread_lane] = Ap[row + thread_lane];
  // 
  // 	 int row_start = ptrs[vector_lane][0]-1;                   //same as: row_start = Ap[row];
  // 	 int row_end   = ptrs[vector_lane][1]-1;                   //same as: row_end   = Ap[row+1];
  // 
  // 	 initialize local sum
  // 	 double sum = 0.00;
  // 
  // 	 accumulate local sums
  // 	 for(int jj = row_start + thread_lane; jj < row_end; jj += THREADS_PER_VECTOR)
  // 	     sum += Ax[jj] * x[ Aj[jj]-1];
  // 
  // 	 store local sum in shared memory
  // 	 sdata[threadIdx.x] = sum;
  // 
  // 	 reduce local sums to row sum
  // 	 if (THREADS_PER_VECTOR > 16) sdata[threadIdx.x] = sum = sum + sdata[threadIdx.x + 16];
  // 	 if (THREADS_PER_VECTOR >  8) sdata[threadIdx.x] = sum = sum + sdata[threadIdx.x +  8];
  // 	 if (THREADS_PER_VECTOR >  4) sdata[threadIdx.x] = sum = sum + sdata[threadIdx.x +  4];
  // 	 if (THREADS_PER_VECTOR >  2) sdata[threadIdx.x] = sum = sum + sdata[threadIdx.x +  2];
  // 	 if (THREADS_PER_VECTOR >  1) sdata[threadIdx.x] = sum = sum + sdata[threadIdx.x +  1];
  // 
  // 	 first thread writes the result
  // 	 if (thread_lane == 0)
  // 	     y[row] += sdata[threadIdx.x];
  //      }
  //  }






  // __global__ void spmv_csr_vector_kernel(int num_rows, int * Ap, int * Aj, double * Ax, double * x, double * y)
  // {
  //     __shared__ volatile double sdata[(VECTORS_PER_BLOCK+1) * THREADS_PER_VECTOR];  // padded to avoid reduction conditionals
  //     __shared__ volatile int ptrs[VECTORS_PER_BLOCK][2];
  //     
  //     const int THREADS_PER_BLOCK = 128;
  // 	//int THREADS_PER_VECTOR = 32;
  // int row_start, row_end;
  //     const int thread_id   = THREADS_PER_BLOCK * blockIdx.x + threadIdx.x;    // global thread index
  //     const int thread_lane = threadIdx.x & (32 - 1);          // thread index within the vector
  //     const int vector_id   = thread_id   /32;               // global vector index
  //     const int vector_lane = threadIdx.x /32;               // vector index within the block
  //     const int num_vectors = 32 * gridDim.x;                   // total number of active vectors
  // 
  //     for(int row = vector_id; row < num_rows; row += num_vectors)
  //     {
  //         // use two threads to fetch Ap[row] and Ap[row+1]
  //         // this is considerably faster than the straightforward version
  //         if(thread_lane < 2)
  //             ptrs[vector_lane][thread_lane] = Ap[row + thread_lane];
  // 
  //         row_start = (ptrs[vector_lane][0])-1;                   //same as: row_start = Ap[row];
  //         row_end   = (ptrs[vector_lane][1])-1;                   //same as: row_end   = Ap[row+1];
  // 
  //         // initialize local sum
  //         double sum = 0;
  //      
  //         if (THREADS_PER_VECTOR == 32 && row_end - row_start > 32)
  //         {
  //             // ensure aligned memory access to Aj and Ax
  // 
  //             int jj = row_start-(row_start&(THREADS_PER_VECTOR-1))+thread_lane;
  // 
  //             // accumulate local sums
  //             if(jj >= row_start && jj < row_end)
  //                 sum += Ax[jj] * x[Aj[jj]-1];
  // 
  //             // accumulate local sums
  //             for(jj += THREADS_PER_VECTOR; jj < row_end; jj += THREADS_PER_VECTOR)
  //                 sum += Ax[jj] * x[Aj[jj]-1];
  //         }
  //         else
  //         {
  //             // accumulate local sums
  //             for(int jj = row_start + thread_lane; jj < row_end; jj += THREADS_PER_VECTOR)
  //                 sum += Ax[jj] * x[Aj[jj]-1];
  //         }
  // 
  //         // store local sum in shared memory
  //         sdata[threadIdx.x] = sum;
  //         
  //         // reduce local sums to row sum
  //         if (THREADS_PER_VECTOR > 16) sdata[threadIdx.x] = sum = sum + sdata[threadIdx.x + 16];
  //         if (THREADS_PER_VECTOR >  8) sdata[threadIdx.x] = sum = sum + sdata[threadIdx.x +  8];
  //         if (THREADS_PER_VECTOR >  4) sdata[threadIdx.x] = sum = sum + sdata[threadIdx.x +  4];
  //         if (THREADS_PER_VECTOR >  2) sdata[threadIdx.x] = sum = sum + sdata[threadIdx.x +  2];
  //         if (THREADS_PER_VECTOR >  1) sdata[threadIdx.x] = sum = sum + sdata[threadIdx.x +  1];
  //        
  //         // first thread writes the result
  //         if (thread_lane == 0)
  //             y[row] += sdata[threadIdx.x];
  //     }
  // }


  // __global__ void spmv_csr_vector_kernel(  int      num_rows,
  //                          int    * start, 
  //                          int    * colid, 
  //                          double * data, 
  //                          double * x, 
  //                          double * y)
  //  {/*{{{*/
  //          int BLOCK_SIZE = 128;
  //  	 int WARP_SIZE = 16;
  //          //typedef typename _Field::Element double ;
  //          __shared__ volatile double sdata[BLOCK_SIZE + 16];                          // padded to avoid reduction ifs
  //          __shared__ volatile int ptrs[BLOCK_SIZE/WARP_SIZE][2];
  //  
  //           int thread_id   = BLOCK_SIZE * blockIdx.x + threadIdx.x;  // global thread index
  //           int thread_lane = threadIdx.x & (WARP_SIZE-1);            // thread index within the warp
  //           int warp_id     = thread_id   / WARP_SIZE;                // global warp index
  //           int warp_lane   = threadIdx.x / WARP_SIZE;                // warp index within the CTA
  //           int num_warps   = (BLOCK_SIZE / WARP_SIZE) * gridDim.x;   // total number of active warps
  //  
  //          for(int row = warp_id; row < num_rows; row += num_warps){
  //                  // use two threads to fetch start[row] and start[row+1]
  //                  // this is considerably faster than the straightforward version
  //                  if(thread_lane < 2)
  //                          ptrs[warp_lane][thread_lane] = start[row + thread_lane];
  //                   int row_start = ptrs[warp_lane][0]-1;                   //same as: row_start = start[row];
  //                   int row_end   = ptrs[warp_lane][1]-1;                   //same as: row_end   = start[row+1];
  //  
  //                  // compute local sum
  //                  double sum = 0;
  //                  int blk_nb  = DIVIDE_INTO(row_end-row_start-thread_lane,WARP_SIZE) ;        
  //                  blk_nb          /= block ;
  //                  int jj       = row_start ;           
  //                  int l        = 0 ;
  //                  int loc_end = row_start ;
  //                  for ( ; l < blk_nb  ; l++) {
  //                          loc_end += block*WARP_SIZE ;
  //                          for ( ; jj < loc_end ; jj += WARP_SIZE )
  //                                  sum += data[jj] * fetch_x<UseCache>(colid[jj], x);
  //                          F_ReduceIn(&sum,p,_Representation());
  //                  } 
  //                  for ( ; jj < row_end ; jj += WARP_SIZE) 
  //                          sum += data[jj] * fetch_x<UseCache>(colid[jj], x);
  //                  F_ReduceIn(&sum,p,_Representation());
  //  
  //                  // reduce local sums to row sum (ASSUME: warpsize 32)
  //                  sdata[threadIdx.x] = sum;
  //                  sdata[threadIdx.x] = sum = sum + sdata[threadIdx.x + 16];
  //                  sdata[threadIdx.x] = sum = sum + sdata[threadIdx.x +  8];
  //                  sdata[threadIdx.x] = sum = sum + sdata[threadIdx.x +  4];
  //                  sdata[threadIdx.x] = sum = sum + sdata[threadIdx.x +  2];
  //                  sdata[threadIdx.x] = sum = sum + sdata[threadIdx.x +  1];
  //  
  //                  // first thread writes warp result
  //                  if (thread_lane == 0){
  //                          y[row] += sdata[threadIdx.x];
  //                          F_ReduceIn(&y[row],p,_Representation());
  //                  }
  //          }
  //  }/*}}}*/
  // 
  // 
  // __inline__ __device__ void F_ReduceIn (double * x,double p, Right)
  // {/*{{{*/
  //          *x = fmodf(*x,p);
  //          if (*x<0) *x+=p ;
  //          return;
  //  }/*}}}*/

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


