// !=============================================================================
// !
// !                              JD SOLVER
// !
// !=============================================================================
// !
// ! Walter Rodrigues
// ! Dipartimento di Ingegneria Elettronica
// ! Universita` di Roma "Tor Vergata"
// ! 25-02-2014
// !
// !=============================================================================
// !
// ! This subroutine finds k harmonic eigenvalues of the operator:
// !









#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<mpi.h>
#include<omp.h>
#include<math.h>
//#include <complex.h>



#include <cuda.h>
#include <cusparse_v2.h>
#include "cublas_v2.h"
#include <cuda_runtime.h>
#include <curand.h>
#include <cuComplex.h>
//#include <cuda_profiler_api.h>


#include <mkl.h>
#include <mkl_lapacke.h>


#include <unistd.h>
#include <sys/unistd.h>
#include <ctype.h>

#include <time.h>
#include <sys/time.h>


#define ENV_LOCAL_RANK		"OMPI_COMM_WORLD_LOCAL_RANK"


#define WARP_SIZE         32 
#define BLOCK_SIZE        1024 
#define BLOCK_SIZE_MUL    256 


#define CLEANUP(s) \
do { \
printf ("%s\n", s); \
if (V0) free(V0);\
fflush (stdout); \
} while (0);



#define CLEANUP_LS(s) \
do { \
printf ("%s\n", s); \
if (h_min) free(h_min);\
fflush (stdout); \
} while (0);



void gmres(float *valptr_device_real, int *rowptr_device_real, int *colptr_device_real, float *valptr_device_img, int *rowptr_device_img, int *colptr_device_img, int n_ham, int size_mat_real, cuDoubleComplex * r_device, double ls_tol, int restart, int maxit, cuDoubleComplex *Q_bar_device, cuDoubleComplex *t_device, int numBlock, int threadPerBlock, int coop, int repeat, int numBlocksMul, int k, cusparseHandle_t cusp_handle_ls, cusparseMatDescr_t cusp_descra_ls, cublasHandle_t cubl_handle_ls, int *ls_counter);


__global__ void vct_div_slr_kernel_1_GPU(cuDoubleComplex * scr, cuDoubleComplex * des, int n_ham, double Slr);
__global__ void vct1_sub_mul_vct_kernel_1_GPU(cuDoubleComplex * vct1, cuDoubleComplex * vct2, int n_ham, cuDoubleComplex * dot);
__global__ void vct_sub_scl_mul_vct_kernel_1_GPU(cuDoubleComplex *vct1, cuDoubleComplex * vct2, double scalar, int n_ham);
__global__ void shift_A_kernel(float * val, int * row, int * col, int n_ham, double shift);
__global__ void cpy_vct_1_to_vct_2_kernel_1_GPU(cuDoubleComplex * vct_scr, cuDoubleComplex * vct_des, int n_ham);
__global__ void vct_1_div_asg_to_vct_2_kernel_1_GPU(cuDoubleComplex * vct_1, cuDoubleComplex * vct_2, int n_ham, double *beta);
__global__ void vct_div_slr_minus_scl_vct_kernel_1_GPU(cuDoubleComplex * scr1, cuDoubleComplex * scr2, cuDoubleComplex * des, int n_ham, double Slr1, double Slr2);
__global__ void vct1_add_vct2_asg_vct3_kernel_1_GPU(cuDoubleComplex *vct1, cuDoubleComplex * vct2, cuDoubleComplex * vct3, int n_ham);
__global__ void vct1_neg_asg_vct2_kernel_1_GPU(cuDoubleComplex *vct1, cuDoubleComplex * vct2, int n_ham);
__global__ void vct1_sub_vct2_asg_vct3_kernel_1_GPU(cuDoubleComplex *vct1, cuDoubleComplex * vct2, cuDoubleComplex * vct3, int n_ham);
__global__ void mv_kernel_1_GPU(int n_ham, int col, cuDoubleComplex * xVal_kr, cuDoubleComplex * y_kr, cuDoubleComplex * Finalans_kr);
__global__ void vct_pls_scl_mul_vct_kernel_1_GPU(cuDoubleComplex *vct1, cuDoubleComplex * vct2, cuDoubleComplex scalar, int n_ham);
__global__ void spmv_csr_hybrid_kernel_1_GPU(int num_rows, const int* rowPtrs, const int* colIdxs, const float* values, const cuDoubleComplex* x, cuDoubleComplex* y, const int offset);
__global__ void spmv_csr_hybrid_kernel_1_GPU(int dimRow, const int* rowPtrs, const int* colIdxs, const float* values, const cuDoubleComplex* x, cuDoubleComplex* y, int repeat, int coop, int offset);


extern "C" {
void jd_single_gpu_no_pc_split_mxprec_(int * N_ham, int * Size_Mat_real, int * Size_Mat_img, float * valptr_real, int * rowptr_real, int * colptr_real, float * valptr_img, int * rowptr_img, int * colptr_img, char * sparse_fmt, int * Band_type, double * JD_tol, double * Shift,  int * JD_Min_step, int * JD_Max_step, int * Num_ev, double * lambda_out, cuDoubleComplex * eigen_vec_out, double * LS_tol, int * LS_restart, int *LS_maxit);

void setdeviceinit_();

} 



void jd_single_gpu_no_pc_split_mxprec_(int * N_ham, int * Size_Mat_real, int * Size_Mat_img, float * valptr_real, int * rowptr_real, int * colptr_real, float * valptr_img, int * rowptr_img, int * colptr_img, char * sparse_fmt, int * Band_type, double * JD_tol, double * Shift,  int *JD_Min_step, int * JD_Max_step, int * Num_ev, double * lambda_out, cuDoubleComplex * eigen_vec_out, double * LS_tol, int * LS_restart, int *LS_maxit)
{

//int * N_ham, *Size_Mat, *JD_Min_step, *JD_Max_step, *Num_ev, *Shift, *LS_tol, *LS_restart, *LS_maxit, *Band_type, *rowptr, *colptr;
//cuDoubleComplex *valptr, eigen_vec_out;
//double * lambda_out, *JD_tol;


int n_ham = * N_ham;
int size_mat_real = * Size_Mat_real;
int size_mat_img = * Size_Mat_img;
int jd_min_step = * JD_Min_step;
int jd_max_step = * JD_Max_step;
int num_ev = * Num_ev;
double shift = * Shift;
double ls_tol = * LS_tol;
int ls_restart = * LS_restart;
int ls_maxit = * LS_maxit;
int band_type = * Band_type;

// printf("n_ham %d\n", n_ham);
// printf("size_mat %d\n", size_mat);
// printf("jd_min_step %d\n", jd_min_step);
// printf("jd_max_step %d\n", jd_max_step);
// printf("num_ev %d\n", num_ev);
// printf("ls_restart %d\n", ls_restart);
// printf("ls_maxit %d\n", ls_maxit);
// printf("band_type %d\n", band_type);
// printf("shift %f\n", shift);
// printf("ls_tol %f\n", ls_tol);
printf("JD_tol %f\n", *JD_tol);

cuDoubleComplex * V0;
cuDoubleComplex * X_bar;
float * valptr_device_real;
int * rowptr_device_real, * colptr_device_real;
float * valptr_device_img;
int * rowptr_device_img, * colptr_device_img;
//cuDoubleComplex * V0_device;
cuDoubleComplex * t_device; 
cuDoubleComplex * X_bar_device;
cuDoubleComplex * vct_device; 
cuDoubleComplex * w_device; 
cuDoubleComplex * v_device;
double * norm_V0_device;
double * norm_vct_device;
double * norm_t_device;
double * norm_t_in_device;
int numBlock, threadPerBlock;
double tol_shift;
cuDoubleComplex  scalar_1, scalar_2;
cuDoubleComplex * scalar1_device;
cuDoubleComplex * scalar2_device;
cuDoubleComplex * scalar1;
cuDoubleComplex * scalar2;
cuDoubleComplex * y_device;
cuDoubleComplex * temp_ev;
cuDoubleComplex * eigen_vec_device;
cuDoubleComplex * u_device;
cuDoubleComplex * eigen_vec;
double * eigen_val;
cuDoubleComplex * M;
cuDoubleComplex * M_device;
double * norm_t_in;
double * norm_vct;
double  * eigen_val_device;
double * s_u_device;
double * s_u;
double lambda[num_ev];
double * norm_r;
double * norm_t;
double * norm_r_device;
cuDoubleComplex * Q_bar;
cuDoubleComplex * Q_bar_device;
//cuDoubleComplex * w_temp_device;
cuDoubleComplex * temp_device;
cuDoubleComplex * w_bar_device;
cuDoubleComplex * r_device;
cuDoubleComplex * u;
double * s_v;
double * s_v_device;
cuDoubleComplex * ubar_device;
double * norm_V0;
const cuDoubleComplex jcmpx=make_cuDoubleComplex(0.0,1.0);
cuDoubleComplex * mxv_temp_device;

////////////////this part needs to be deleted////////////////////////////////
//cuDoubleComplex * u_bar;
//u_bar = (cuDoubleComplex *) malloc(sizeof(cuDoubleComplex)*n_ham);
cuDoubleComplex * v, * w;
//v = (cuDoubleComplex *) malloc(n_ham*jd_max_step*sizeof(cuDoubleComplex));
//w = (cuDoubleComplex *) malloc(n_ham*jd_max_step*sizeof(cuDoubleComplex));
cudaMallocHost((void**)&v, n_ham*jd_max_step*sizeof(cuDoubleComplex));
cudaMallocHost((void**)&w, n_ham*jd_max_step*sizeof(cuDoubleComplex));
//cuDoubleComplex *w_bar;
//w_bar = (cuDoubleComplex *) malloc(sizeof(cuDoubleComplex)*n_ham);
//cuDoubleComplex *r_cpu;
//r_cpu = (cuDoubleComplex *) malloc(sizeof(cuDoubleComplex)*n_ham);
//cuDoubleComplex *t_cpu;
//t_cpu = (cuDoubleComplex *) malloc(sizeof(cuDoubleComplex)*n_ham);
/////////////////////////////////////////////////////

int *ls_counter;
ls_counter = (int *) malloc(sizeof(int));
*ls_counter = 0;

lapack_int IL=0, IU=0, *M_out, *isuppz;
double  ABSTOL = 0.001;
double VL=0, VU=0;


norm_r = (double *) malloc(sizeof(double));
norm_t = (double *) malloc(sizeof(double));
u = (cuDoubleComplex *) malloc(sizeof(cuDoubleComplex)*n_ham);
norm_vct = (double *) malloc(sizeof(double));
norm_t_in = (double *) malloc(sizeof(double));
//norm_r = (double *) malloc(sizeof(double));
M_out = (lapack_int *) malloc(sizeof(lapack_int));
isuppz = (lapack_int *) malloc(2*jd_max_step*sizeof(lapack_int));
s_u = (double *) malloc(sizeof(double));
s_v = (double *) malloc(sizeof(double));
scalar1 = (cuDoubleComplex *) malloc(sizeof(cuDoubleComplex));
scalar2 = (cuDoubleComplex *) malloc(sizeof(cuDoubleComplex));
norm_V0 = (double *) malloc(sizeof(double));

cudaError_t cudaStat1,cudaStat2,cudaStat3,cudaStat4,cudaStat5,cudaStat6,cudaStat7,cudaStat8,cudaStat9,cudaStat10,cudaStat11,cudaStat12,cudaStat13,cudaStat14,cudaStat15,cudaStat16, cudaStat17,cudaStat18,cudaStat19,cudaStat20,cudaStat21,cudaStat22,cudaStat23,cudaStat24,cudaStat25,cudaStat26,cudaStat27;

lapack_int error_code=0;

cusparseStatus_t cusp_status;
cusparseHandle_t cusp_handle=0;
cusparseMatDescr_t cusp_descra=0;
cublasStatus_t cubl_status;
cublasHandle_t cubl_handle=0;

setdeviceinit_();

//cudaDeviceReset();
//cudaThreadExit();	
//int devCount = 0;	
//cudaGetDeviceCount(&devCount);
//printf("device count %d %d %d\n", devCount, id, id%devCount);
//cudaStat1 = cudaSetDevice(1);
//if(cudaStat1 != cudaSuccess)
//printf("ERROR DEVICE SET FAILED\n");
//cudaDeviceReset();

//cudaProfilerStart();

//threadPerBlock = 1024;
//numBlock=(n_ham/threadPerBlock)+1;

threadPerBlock = BLOCK_SIZE;                      // blockSize
numBlock=(n_ham/threadPerBlock)+1;                // gridSize

int coop= 16  ;    
int repeat= 2 ;   // repeat = nrows * coop/ BLOCK_SIZE_MUL 
int numBlocksMul =(n_ham*coop-1)/(repeat*BLOCK_SIZE_MUL)+1;  //gridSize
 
  
printf(" (CUDA) JD single GPU, threadsPerBlock %d \n", threadPerBlock);
printf(" (CUDA) gridSize  %d \n",numBlock);
printf(" (CUDA) JD single, threadsPerBlock %d \n", BLOCK_SIZE_MUL);
printf(" (CUDA) gridSize for Mxv %d \n",numBlocksMul);



scalar_1.x = 1.00;
scalar_1.y = 0.00;
scalar_2.x = 0.00;
scalar_2.y = 0.00;

scalar1[0].x = 1.00;
scalar1[0].y = 0.00;
scalar2[0].x = 0.00;
scalar2[0].y = 0.00;


int flag =1;
V0 = (cuDoubleComplex *) malloc(sizeof(cuDoubleComplex)*n_ham);

for(int counter=0; counter < n_ham; counter++)
{
	V0[counter].x = rand();
	V0[counter].y = rand();
}

X_bar = (cuDoubleComplex *) malloc(sizeof(cuDoubleComplex)*n_ham*num_ev);



printf("inside JD\n");

double memory_device = (sizeof(float)*(size_mat_real+size_mat_img)+sizeof(int)*(size_mat_real+n_ham+size_mat_img+n_ham)+sizeof(double)*(jd_max_step)+sizeof(cuDoubleComplex)*(n_ham*12+n_ham*num_ev+n_ham*jd_max_step*3)+(n_ham*((ls_restart+1)+ls_restart))+(jd_max_step*jd_max_step*2)+((ls_restart+1)*(ls_restart+1))+((ls_restart+1)*ls_restart))/1000000;

printf("memory needed on GPU is ~ %f MB\n", memory_device); 


cudaStat1 = cudaMalloc((void**)&colptr_device_real,size_mat_real*sizeof(int));
cudaStat2 = cudaMalloc((void**)&rowptr_device_real,(n_ham+1)*sizeof(int));
cudaStat3 = cudaMalloc((void**)&valptr_device_real,size_mat_real*sizeof(float));
//cudaStat4 = cudaMalloc((void**)&V0_device,n_ham*sizeof(cuDoubleComplex));
cudaStat6 = cudaMalloc((void**)&t_device,n_ham*sizeof(cuDoubleComplex));
cudaStat9 = cudaMalloc((void**)&y_device,sizeof(cuDoubleComplex));
//cudaStat10 = cudaMalloc((void**)&w_device,n_ham*jd_max_step*sizeof(cuDoubleComplex));
cudaStat11 = cudaMalloc((void**)&colptr_device_img,size_mat_img*sizeof(int));
cudaStat14 = cudaMalloc((void**)&rowptr_device_img,(n_ham+1)*sizeof(int));
cudaStat16 = cudaMalloc((void**)&valptr_device_img,size_mat_img*sizeof(float));
//cudaStat17 = cudaMalloc((void**)&w_bar_device,n_ham*sizeof(cuDoubleComplex));
//cudaStat18 = cudaMalloc((void**)&w_temp_device,n_ham*jd_max_step*sizeof(cuDoubleComplex));
//cudaStat22 = cudaMalloc((void**)&ubar_device,n_ham*sizeof(cuDoubleComplex));
cudaStat25 = cudaMalloc((void**)&scalar1_device,sizeof(cuDoubleComplex));
cudaStat26 = cudaMalloc((void**)&scalar2_device,sizeof(cuDoubleComplex));
cudaStat27 = cudaMalloc((void**)&s_v_device,sizeof(double));
cudaStat24 = cudaMalloc((void**)&norm_r_device,sizeof(double));
cudaStat12 = cudaMalloc((void**)&norm_vct_device,sizeof(double));
cudaStat13 = cudaMalloc((void**)&norm_t_in_device,sizeof(double));
cudaStat15 = cudaMalloc((void**)&s_u_device,sizeof(double));
cudaStat21 = cudaMalloc((void**)&norm_t_device,sizeof(double));
cudaStat5 = cudaMalloc((void**)&norm_V0_device,sizeof(double));
cudaStat23 = cudaMalloc((void**)&r_device,n_ham*sizeof(cuDoubleComplex));
if( 
(cudaStat1 != cudaSuccess) ||
(cudaStat2 != cudaSuccess) ||
(cudaStat3 != cudaSuccess) ||
//(cudaStat4 != cudaSuccess) ||
(cudaStat5 != cudaSuccess) ||
(cudaStat6 != cudaSuccess) ||
//(cudaStat7 != cudaSuccess) ||
//(cudaStat8 != cudaSuccess) ||
(cudaStat9 != cudaSuccess) ||
//(cudaStat10 != cudaSuccess) ||
(cudaStat11 != cudaSuccess) ||
(cudaStat14 != cudaSuccess) ||
(cudaStat12 != cudaSuccess) ||
(cudaStat13 != cudaSuccess) ||
(cudaStat15 != cudaSuccess) ||
(cudaStat16 != cudaSuccess) ||
//(cudaStat17 != cudaSuccess) ||
//(cudaStat18 != cudaSuccess) ||
//(cudaStat19 != cudaSuccess) ||
(cudaStat21 != cudaSuccess) ||
//(cudaStat22 != cudaSuccess) ||
(cudaStat23 != cudaSuccess) ||
(cudaStat24 != cudaSuccess) ||
(cudaStat25 != cudaSuccess) ||
(cudaStat26 != cudaSuccess) ||
(cudaStat27 != cudaSuccess))
{
CLEANUP("Device malloc failed loc 1\n");
}

//printf("cuda memory alloc done\n");

cudaStat1 = cudaMemcpy(colptr_device_real, colptr_real, (size_t)(size_mat_real*sizeof(int)), cudaMemcpyHostToDevice);
cudaStat2 = cudaMemcpy(rowptr_device_real, rowptr_real, (size_t)((n_ham+1)*sizeof(int)), cudaMemcpyHostToDevice);
cudaStat3 = cudaMemcpy(valptr_device_real, valptr_real, (size_t)(size_mat_real*sizeof(float)), cudaMemcpyHostToDevice);
cudaStat4 = cudaMemcpy(t_device, V0, (size_t)(n_ham*sizeof(cuDoubleComplex)), cudaMemcpyHostToDevice);
cudaStat5 = cudaMemcpy(scalar1_device, scalar1, (size_t)(sizeof(cuDoubleComplex)), cudaMemcpyHostToDevice);
cudaStat6 = cudaMemcpy(scalar2_device, scalar2, (size_t)(sizeof(cuDoubleComplex)), cudaMemcpyHostToDevice);
cudaStat7 = cudaMemcpy(colptr_device_img, colptr_img, (size_t)(size_mat_img*sizeof(int)), cudaMemcpyHostToDevice);
cudaStat8 = cudaMemcpy(rowptr_device_img, rowptr_img, (size_t)((n_ham+1)*sizeof(int)), cudaMemcpyHostToDevice);
cudaStat9 = cudaMemcpy(valptr_device_img, valptr_img, (size_t)(size_mat_img*sizeof(float)), cudaMemcpyHostToDevice);

if(
(cudaStat1 != cudaSuccess) ||
(cudaStat2 != cudaSuccess) ||
(cudaStat3 != cudaSuccess) ||
(cudaStat4 != cudaSuccess) ||
(cudaStat5 != cudaSuccess) ||
(cudaStat6 != cudaSuccess) ||
(cudaStat7 != cudaSuccess) ||
(cudaStat8 != cudaSuccess) ||
(cudaStat9 != cudaSuccess))
{
     CLEANUP("copy to device failed loc 2\n");
}


//printf("cuda memory copy done\n");

// initalization of CUBLAS library
cubl_status = cublasCreate(&cubl_handle);
if(cubl_status != CUBLAS_STATUS_SUCCESS)
{
    CLEANUP("CUBLAS initilization failed \n");
}

cubl_status = cublasSetPointerMode(cubl_handle, CUBLAS_POINTER_MODE_DEVICE);
if(cubl_status != CUBLAS_STATUS_SUCCESS)
{
    CLEANUP("CUBLAS setting device pointer mode failed\n");
}

/* initialize cusparse library */
cusp_status= cusparseCreate(&cusp_handle);
if (cusp_status != CUSPARSE_STATUS_SUCCESS) 
{
    CLEANUP("CUSPARSE Library initialization failed\n");
}

 /* create and setup matrix descriptor */
cusp_status= cusparseCreateMatDescr(&cusp_descra);
if (cusp_status != CUSPARSE_STATUS_SUCCESS) 
{
    CLEANUP("Matrix descriptor initialization failed 1");
}

cusp_status=cusparseSetMatType(cusp_descra, CUSPARSE_MATRIX_TYPE_GENERAL);
if (cusp_status != CUSPARSE_STATUS_SUCCESS) 
{
    CLEANUP("Matrix descriptor initialization failed 2");
}
        
cusp_status=cusparseSetMatIndexBase(cusp_descra, CUSPARSE_INDEX_BASE_ONE);
if (cusp_status != CUSPARSE_STATUS_SUCCESS) 
{
    CLEANUP("Matrix descriptor initialization failed 3");
}

cublasOperation_t transaction;
transaction = CUBLAS_OP_N;


//printf("cublas and cusparse initilized\n");



// V0 = V0/norm(V0); t_device = V0
cubl_status = cublasDznrm2(cubl_handle, n_ham, t_device, 1, norm_V0_device);
cudaStat1 = cudaMemcpy(norm_V0, norm_V0_device,(size_t)(sizeof(double)), cudaMemcpyDeviceToHost);

//printf("Norm of V0 done\n");

vct_div_slr_kernel_1_GPU<<<numBlock, threadPerBlock>>>( t_device, t_device, n_ham, *norm_V0);

//printf("V0 = V0/norm(V0) done\n");

// t=v0; k = 0; m = 0; X_bar = [];
int k = 0;
int m = 0;

//cpy_vct_1_to_vct_2_kernel_1_GPU<<<numBlock, threadPerBlock>>>( V0_device, t_device, n_ham);

//printf("t_device = V0 done\n");

*s_v = 0;
int count = 0;
double teta = 0;


while(k < num_ev)
{

//printf("inside while\n");



     // A = (A - shift * I)

    tol_shift = teta-shift;
    shift_A_kernel<<<numBlock, threadPerBlock>>>(valptr_device_real, rowptr_device_real, colptr_device_real, n_ham, tol_shift);

// printf("total shift = %f\n", tol_shift);
// 
// cuDoubleComplex *val_A;
// int *row_A, *col_A;
// val_A = (cuDoubleComplex *) malloc(size_mat*sizeof(cuDoubleComplex));
// col_A = (int *) malloc(size_mat*sizeof(int));
// row_A = (int *) malloc((n_ham+1)*sizeof(int));
// cudaStat1 = cudaMemcpy(val_A, valptr_device,(size_t)(size_mat*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);
// 
// 
// printf("val_A\n");
// printf("val_A[0].x=%f val[0].y=%f\n", val_A[0].x, val_A[0].y); 
// printf("val_A[1].x=%f val[1].y=%f\n", val_A[1].x, val_A[1].y);
// printf("val_A[2].x=%f val[2].y=%f\n", val_A[2].x, val_A[2].y);
// printf("val_A[3].x=%f val[3].y=%f\n", val_A[3].x, val_A[3].y);
// printf("val_A[4].x=%f val[4].y=%f\n", val_A[4].x, val_A[4].y);
// printf("val_A[5].x=%f val[5].y=%f\n", val_A[5].x, val_A[5].y);
// printf("val_A[6].x=%f val[6].y=%f\n", val_A[6].x, val_A[6].y);
// printf("val_A[7].x=%f val[7].y=%f\n", val_A[7].x, val_A[7].y);
// printf("val_A[8].x=%f val[8].y=%f\n", val_A[8].x, val_A[8].y);
// printf("val_A[9].x=%f val[9].y=%f\n", val_A[9].x, val_A[9].y);
// printf("val_A[10].x=%f val[10].y=%f\n", val_A[10].x, val_A[10].y);
// printf("val_A[11].x=%f val[11].y=%f\n", val_A[11].x, val_A[11].y);

    

//      for(int i=1; i <= n_ham; i++)
//      {
//       for(int j=rowptr_device(i); j <= rowptr_device(i+1)-1)
//       {
//        if(colptr_device(j) == i)
//        {
//          valptr_device(j) = valptr_device(j) - shift + teta; 
//          break;
//        }
//       }
//      }
        
      cudaStat1 = cudaMalloc((void**)&vct_device,n_ham*sizeof(cuDoubleComplex));
      cudaStat2 = cudaMalloc((void**)&mxv_temp_device,n_ham*sizeof(cuDoubleComplex));
      if(cudaStat1 != cudaSuccess || cudaStat2 != cudaSuccess)
      {
         CLEANUP("Device malloc failed loc 5\n");
      }

      // vct = A * t;
      //cusp_status= cusparseZcsrmv(cusp_handle,CUSPARSE_OPERATION_NON_TRANSPOSE, n_ham, n_ham, size_mat_real, &scalar_1, cusp_descra, valptr_device, rowptr_device, colptr_device, t_device, &scalar_2, vct_device);

     spmv_csr_hybrid_kernel_1_GPU<<<numBlocksMul, BLOCK_SIZE_MUL>>>(n_ham, rowptr_device_real, colptr_device_real, valptr_device_real, t_device, vct_device, repeat, coop, 0);

     spmv_csr_hybrid_kernel_1_GPU<<<numBlock, BLOCK_SIZE>>>(n_ham, rowptr_device_img, colptr_device_img, valptr_device_img, t_device, mxv_temp_device, 0);

     vct_pls_scl_mul_vct_kernel_1_GPU<<<numBlock, threadPerBlock>>>(vct_device, mxv_temp_device, jcmpx, n_ham);

     cudaFree(mxv_temp_device);

//printf("mxv done\n");


     // t_in = norm(vct);
     cubl_status = cublasDznrm2(cubl_handle, n_ham, vct_device, 1, norm_t_in_device);

//printf("norm(vct) done\n");

     cudaStat1 = cudaMemcpy(norm_t_in, norm_t_in_device,(size_t)(sizeof(double)), cudaMemcpyDeviceToHost);

//printf("memcpy done\n");
//printf("norm_t_in %f\n", *norm_t_in);


    cudaStat1 = cudaMalloc((void**)&w_device,n_ham*jd_max_step*sizeof(cuDoubleComplex));
    cudaStat2 = cudaMalloc((void**)&v_device,n_ham*jd_max_step*sizeof(cuDoubleComplex));
    if(cudaStat1 != cudaSuccess || cudaStat1 != cudaSuccess)
      {
         CLEANUP("Device malloc failed loc 5\n");
      }

    cudaStat1 = cudaMemcpy(v_device, v,(size_t)(sizeof(cuDoubleComplex)*n_ham*m), cudaMemcpyHostToDevice);
    cudaStat1 = cudaMemcpy(w_device, w,(size_t)(sizeof(cuDoubleComplex)*n_ham*m), cudaMemcpyHostToDevice);
     
    for(int i = 0; i < m; i++)
    {
     // y = w(:,i)' * vct;
     cubl_status = cublasZdotc(cubl_handle, n_ham, &w_device[i*n_ham], 1, vct_device, 1, y_device); 

//printf("y = w(:,i)' * vct\n");

     // vct = vct - y * w(:,i);
     vct1_sub_mul_vct_kernel_1_GPU<<<numBlock, threadPerBlock>>>(vct_device, &w_device[i*n_ham], n_ham, y_device);


//printf("vct = vct - y * w(:,i)\n");

     // t = t - y * v(:,i);
     vct1_sub_mul_vct_kernel_1_GPU<<<numBlock, threadPerBlock>>>(t_device, &v_device[i*n_ham], n_ham, y_device);

//printf("t = t - y * v(:,i)\n");
     }

     // norm(vct);
     cubl_status = cublasDznrm2(cubl_handle, n_ham, vct_device, 1, norm_vct_device);

     cudaStat1 = cudaMemcpy(norm_vct, norm_vct_device,(size_t)(sizeof(double)), cudaMemcpyDeviceToHost); 

//printf("norm_vct %f\n", *norm_vct);

     if( ((*norm_vct)/(*norm_t_in)) < 0.25)
     {
         for(int i = 0; i < m; i++)
         {
         // y = w(:,i)' * vct;
         cubl_status = cublasZdotc(cubl_handle, n_ham, &w_device[i*n_ham], 1, vct_device, 1, y_device); 


        // vct = vct - y * w(:,i);
        vct1_sub_mul_vct_kernel_1_GPU<<<numBlock, threadPerBlock>>>(vct_device, &w_device[i*n_ham], n_ham, y_device);


       // t = t - y * v(:,i);
       vct1_sub_mul_vct_kernel_1_GPU<<<numBlock, threadPerBlock>>>(t_device, &v_device[i*n_ham], n_ham, y_device);
        }
     }

// cudaStat1 = cudaMemcpy(t_cpu, t_device,(size_t)(sizeof(cuDoubleComplex)*n_ham), cudaMemcpyDeviceToHost);
// printf("t_cpu[0].x=%f t_cpu[0].y=%f\n", t_cpu[0].x, t_cpu[0].y); 
// printf("t_cpu[1].x=%f t_cpu[1].y=%f\n", t_cpu[1].x, t_cpu[1].y);
// printf("t_cpu[2].x=%f t_cpu[2].y=%f\n", t_cpu[2].x, t_cpu[2].y);
// printf("t_cpu[3].x=%f t_cpu[3].y=%f\n", t_cpu[3].x, t_cpu[3].y);

//printf("MGS done\n");

     m = m+1;
   
     // norm(vct);cudaError_t cudaStat1,cudaStat2,cudaStat3,cudaStat4,cudaStat5,cudaStat6, cudaStat7,cudaStat8,cudaStat9,cudaStat10,cudaStat11,cudaStat12,cudaStat13,cudaStat14,cudaStat15
     cubl_status = cublasDznrm2(cubl_handle, n_ham, vct_device, 1, norm_vct_device);

     cudaStat1 = cudaMemcpy(norm_vct, norm_vct_device,(size_t)(sizeof(double)), cudaMemcpyDeviceToHost);

//printf("norm_vct %f\n", *norm_vct);


    // w(:,m) = vct/norm(vct);vct_div_slr_kernel_1_GPU(cuDoubleComplex * scr, cuDoubleComplex * des, int n_ham, double * Slr)
    // v(:,m) = t/norm(vct);
    vct_div_slr_kernel_1_GPU<<<numBlock, threadPerBlock>>>(vct_device, &w_device[(m-1)*n_ham], n_ham, *norm_vct);
    vct_div_slr_kernel_1_GPU<<<numBlock, threadPerBlock>>>(t_device, &v_device[(m-1)*n_ham], n_ham, *norm_vct);

    cudaFree(vct_device);

// cudaStat1 = cudaMemcpy(v, v_device,(size_t)(sizeof(cuDoubleComplex)*n_ham*m), cudaMemcpyDeviceToHost);    
// for(int i =0; i < m*n_ham; i++){ 
// printf("v[%d].x=%f v[%d].y=%f\n", i, v[i].x, i, v[i].y); 
// }

//cudaStat2 = cudaMemcpy(w, w_device,(size_t)(n_ham*m*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);
// if(cudaStat1 != cudaSuccess){
// CLEANUP("cudaMemcpy failed\n");}
 
//  for(int i =0; i < m*n_ham; i++){ 
//  printf("w[%d].x=%f w[%d].y=%f\n", i, w[i].x, i, w[i].y); 
//  }
    
//     if(m != 1) 
//       {
//          cudaFree(M_device);
//          cudaFree(eigen_vec_device);
//          cudaFree(eigen_val_device);
//          free(M);
//          free(eigen_val);
//          free(eigen_vec);
//       }
    
    // alloc M(m,m);
    cudaStat1 = cudaMalloc((void**)&M_device,m*m*sizeof(cuDoubleComplex));

    M = (cuDoubleComplex *) malloc(sizeof(cuDoubleComplex)*m*m);

    for(int j=0; j < m; j++)
    { 

    for(int i = 0; i < m; i++)
    {
     // M(i,m) = w(:,i)' * v(:,m);
     cubl_status = cublasZdotc(cubl_handle, n_ham, &w_device[i*n_ham], 1, &v_device[j*n_ham], 1, &M_device[j*m+i]);
     // M(m,i) = w(:,m)' * v(:,i);
     cubl_status = cublasZdotc(cubl_handle, n_ham, &w_device[j*n_ham], 1, &v_device[i*n_ham], 1, &M_device[m*i+j]);
    }

    // M(m,m) = w(:,m)' * v(:,m);
    cubl_status = cublasZdotc(cubl_handle, n_ham, &w_device[j*n_ham], 1, &v_device[j*n_ham], 1, &M_device[j*m+j]);

    }

    cudaStat1 = cudaMemcpy(M, M_device,(size_t)(sizeof(cuDoubleComplex)*m*m), cudaMemcpyDeviceToHost);
    if(cudaStat1 != cudaSuccess)
    {
     CLEANUP("cudaMemcpy failed loc 3\n");
    }

// for(int ii = 0; ii < m*m; ii++){
// printf("M[%d] = %f %f\n", ii, M[ii].x, M[ii].y); }
    
    eigen_vec = (cuDoubleComplex *) malloc(sizeof(cuDoubleComplex)*m*m);
    eigen_val = (double *) malloc(sizeof(double)*m);
    cudaStat1 = cudaMalloc((void**)&eigen_vec_device,sizeof(cuDoubleComplex)*m*m);
    cudaStat2 = cudaMalloc((void**)&eigen_val_device,sizeof(double)*m);
    if(cudaStat1 != cudaSuccess || cudaStat2 != cudaSuccess){
      CLEANUP("cudaMemcpy failed loc 4\n");}
//printf("memory for M allocated done\n");

    
    //error_code = LAPACKE_zheev('V', 'U', m, M, m, VL, VU, IL, IU, ABSTOL, M_out, eigen_val, eigen_vec, m, Ifail );
 
      //error_code = LAPACKE_zheev(LAPACK_COL_MAJOR, 'V', 'U', m, M, m, m, eigen_val);

     // schur(M)
     // lapack_int LAPACKE_zheevr( int matrix_layout, char jobz, char range, char uplo, lapack_int n, lapack_complex_double* a, lapack_int lda, double vl, double vu, lapack_int il, lapack_int iu, double abstol, lapack_int* m, double* w, lapack_complex_double* z, lapack_int ldz, lapack_int* isuppz );

     //if(m != 1)
     //{
       error_code = LAPACKE_zheevr( LAPACK_COL_MAJOR, 'V', 'A', 'L', m, (MKL_Complex16 *)M, m, VL, VU, IL, IU, ABSTOL, M_out, eigen_val, (MKL_Complex16 *)eigen_vec, m, isuppz );
     //}
//      else
//      {
//        *M_out = 1;
//        eigen_val[0] = M[0].x;
//        eigen_vec[0].x = 1.000000; eigen_vec[0].y = 0.000000;
//      }

// printf("schur done\n");
// printf("error code %d\n", error_code); 

     if(m != 1)
     {
     int index[m], temp_index;
     double temp;
     temp_ev = (cuDoubleComplex *) malloc(sizeof(cuDoubleComplex)*m);

     for(int i =0; i < m; i++)
     {
       index[i] = i;
     }
     
// printf("inside sort\n");
// for(int i =0; i < m; i++){
// printf("index[%d] = %d eigen_val[%d] = %f\n",i, index[i], i, eigen_val[i]);
// }
// for(int i =0; i < m*m; i++){
// printf("eigen_vec[%d].x = %f eigen_vec[%d].y = %f\n", i, eigen_vec[i].x, i, eigen_vec[i].y );
// }

     // sort(eigen_val - shift)  band_type = 1 for CB and band_type = 2 for VB
     if(band_type == 1)
     {
     for(int i =0; i < m-1; i++)
     {
      for(int j = 0; j < m-1-i; j++)
      {
         if((eigen_val[j]-shift) < (eigen_val[j+1]-shift))
         {
            temp = eigen_val[j];
            eigen_val[j] = eigen_val[j+1];
            eigen_val[j+1] = temp;
            temp_index = index[j];
            index[j] = index[j+1];
            index[j+1] = temp_index;
         }
      }
     }
     }
     else if(band_type == 2)
     {
      for(int i =0; i < m-1; i++)
     {
      for(int j = 0; j < m-1-i; j++)
      {
         if((eigen_val[j]-shift) > (eigen_val[j+1]-shift))
         {
            temp = eigen_val[j];
            eigen_val[j] = eigen_val[j+1];
            eigen_val[j+1] = temp;
            temp_index = index[j];
            index[j] = index[j+1];
            index[j+1] = temp_index;
         }
      }
     }
     }

// printf("inside sort\n");
// for(int i =0; i < m; i++){
// printf("index[%d] = %d\n", i, index[i]);
// }

      // sort eigen_vector
      for(int i =0; i < m; i++)
      {
         //temp_index = index[i];
         //memcpy( &temp_ev[0], &eigen_vec[i*m], (sizeof(cuDoubleComplex)*m) );
         //memcpy( &eigen_vec[i*m], &eigen_vec[temp_index*m], (sizeof(cuDoubleComplex)*m) );
         //memcpy( &eigen_vec[temp_index*m], &temp_ev[0], (sizeof(cuDoubleComplex)*m) );

         cudaStat1 = cudaMemcpy(&eigen_vec_device[i*m], &eigen_vec[index[i]*m],(size_t)(sizeof(cuDoubleComplex)*m), cudaMemcpyHostToDevice);
         if(cudaStat1 != cudaSuccess)
         {
         CLEANUP("cudaMemcpy failed loc 5\n");
         }

      }
      }

//cudaStat1=cudaDeviceSynchronize();
//if(cudaStat1 != cudaSuccess){
//CLEANUP("device sync failed 1\n");
//}

      if(m == 1) 
      { 
        cudaStat1 = cudaMemcpy(eigen_vec_device, eigen_vec,(size_t)(sizeof(cuDoubleComplex)*m*m), cudaMemcpyHostToDevice);
        if(cudaStat1 != cudaSuccess){
        CLEANUP("cudaMemcpy failed loc 6\n");}
      }
      else {
      cudaStat1 = cudaMemcpy(eigen_vec, eigen_vec_device,(size_t)(sizeof(cuDoubleComplex)*m*m), cudaMemcpyDeviceToHost);
      if(cudaStat1 != cudaSuccess){
        CLEANUP("cudaMemcpy failed loc 7\n");}
      }
      cudaStat2 = cudaMemcpy(eigen_val_device, eigen_val,(size_t)(sizeof(double)*m), cudaMemcpyHostToDevice);
      if(cudaStat2 != cudaSuccess)
      {
         CLEANUP("cudaMemcpy failed loc 8\n");
      }


//printf("outside sort\n");
// for(int i =0; i < m; i++){
// printf("eigen_val[%d] = %f\n", i, eigen_val[i]);
// }
// for(int i =0; i < (m*m); i++){
// printf("eigen_vec[%d].x = %f eigen_vec[%d].y = %f\n", i, eigen_vec[i].x, i, eigen_vec[i].y );
// }
// 
// 
// printf(" i m here 1\n");
// 
// 
//  cudaStat1 = cudaMemcpy(v, v_device,(size_t)(n_ham*m*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);
//  cudaStat2 = cudaMemcpy(w, w_device,(size_t)(n_ham*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);
//  if(cudaStat1 != cudaSuccess){
//  CLEANUP("cudaMemcpy failed\n");}
//  
//  for(int i =0; i < m*n_ham; i++){ 
//  printf("v[%d].x=%f v[%d].y=%f\n", i, v[i].x, i, v[i].y); 
//  }
//  
//  cudaStat1 = cudaMemcpy(eigen_vec, eigen_vec_device,(size_t)(m*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);
//  if(cudaStat1 != cudaSuccess){
//  CLEANUP("cudaMemcpy failed\n");}
//  for(int i =0; i < m; i++){ 
//  printf("eigen_vec[%d].x=%f eigen_vec[%d].y=%f\n", i, eigen_vec[i].x, i, eigen_vec[i].y); 
//  }

// cudaStat1=cudaDeviceSynchronize();
//        if(cudaStat1 != cudaSuccess){
//       CLEANUP("device sync failed 1\n");
//       }
// 
//  cudaStat1=cudaDeviceSynchronize();
//        if(cudaStat1 != cudaSuccess){
//       CLEANUP("device sync failed 2\n");
//       }

//  printf("w[0].x=%f w[0].y=%f\n", w[0].x, w[0].y); 
//  printf("w[1].x=%f w[1].y=%f\n", w[1].x, w[1].y);
//  printf("w[2].x=%f w[2].y=%f\n", w[2].x, w[2].y);
//  printf("w[3].x=%f w[3].y=%f\n", w[3].x, w[3].y);


      cudaStat1 = cudaMalloc((void**)&ubar_device,n_ham*sizeof(cuDoubleComplex));
      if(cudaStat1 != cudaSuccess)
      {
      CLEANUP("Device malloc failed\n");
      }

      // u_bar = v * eig_mat(:,1);
      cubl_status = cublasZgemv(cubl_handle, CUBLAS_OP_N, n_ham, m, scalar1_device, v_device, n_ham, eigen_vec_device, 1, scalar2_device, ubar_device, 1);
      if(cubl_status != CUBLAS_STATUS_SUCCESS){
      //printf("%d\n",cubl_status);
      CLEANUP("CUBLAS Zgemv failed loc 9\n");
      }




        //mv_kernel_1_GPU<<<numBlock, threadPerBlock>>>(n_ham, m, v_device, eigen_vec_device, ubar_device);

//       cudaStat1=cudaDeviceSynchronize();
//       if(cudaStat1 != cudaSuccess){
//       printf("%s\n", cudaGetErrorString(cudaStat1));
//        CLEANUP("device sync failed 3\n");
//      }
        
//          cudaStat1 = cudaMemcpy(ubar_device, result_device, (size_t)(n_ham*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToDevice);
//          if(cudaStat1 != cudaSuccess)
//          {
//           //printf("%s\n", cudaGetErrorString(cudaStat1));
//           CLEANUP("cudaMemcpy failed test\n");
//          }

//          cudaStat1 = cudaMemcpy(u_bar, ubar_device,(size_t)(n_ham*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);
//          if(cudaStat1 != cudaSuccess)
//          {
//           //printf("%s\n", cudaGetErrorString(cudaStat1));
//           CLEANUP("cudaMemcpy failed\n");
//          }
//  
//          printf("u_bar[0].x=%f u_bar[0].y=%f\n", u_bar[0].x, u_bar[0].y); 
//          printf("u_bar[1].x=%f u_bar[1].y=%f\n", u_bar[1].x, u_bar[1].y);
//          printf("u_bar[2].x=%f u_bar[2].y=%f\n", u_bar[2].x, u_bar[2].y);
//          printf("u_bar[3].x=%f u_bar[3].y=%f\n", u_bar[3].x, u_bar[3].y);
//         //free(u_bar);
// 
// 
// printf(" i m here 2\n");

      // s_u = norm(u_bar);
      cubl_status = cublasDznrm2(cubl_handle, n_ham, ubar_device, 1, s_u_device);

      cudaStat1 = cudaMemcpy(s_u, s_u_device,(size_t)(sizeof(double)), cudaMemcpyDeviceToHost);
      if(cudaStat1 != cudaSuccess)
      {
         CLEANUP("cudaMemcpy failed loc 10\n");
      }

//printf("s_u %f\n", *s_u);

      cudaStat1 = cudaMalloc((void**)&u_device,n_ham*sizeof(cuDoubleComplex));
      if(cudaStat1 != cudaSuccess)
      {
         CLEANUP("Device malloc failed loc 10\n");
      }

      // u = u_bar/s_u; 
      vct_div_slr_kernel_1_GPU<<<numBlock, threadPerBlock>>>(ubar_device, u_device, n_ham, *s_u);

      cudaFree(ubar_device);

//         cuDoubleComplex *u_cpu;
//         u_cpu = (cuDoubleComplex *) malloc(sizeof(cuDoubleComplex)*n_ham);
//         cudaStat1 = cudaMemcpy(u_cpu, u_device,(size_t)(sizeof(cuDoubleComplex)*n_ham), cudaMemcpyDeviceToHost);
//         printf("u[0].x=%f u[0].y=%f\n", u_cpu[0].x, u_cpu[0].y); 
//         printf("u[1].x=%f u[1].y=%f\n", u_cpu[1].x, u_cpu[1].y);
//         printf("u[2].x=%f u[2].y=%f\n", u_cpu[2].x, u_cpu[2].y);
//         printf("u[3].x=%f u[3].y=%f\n", u_cpu[3].x, u_cpu[3].y);

      // s_v = eig_val(1)/ (s_u * s_u);
      *s_v = eigen_val[0]/( (*s_u) * (*s_u));

//printf("s_v %f\n", *s_v);

      cudaStat1 = cudaMemcpy(s_v_device, s_v,(size_t)(sizeof(double)), cudaMemcpyHostToDevice);
      if(cudaStat1 != cudaSuccess)
      {
         CLEANUP("cudaMemcpy failed loc 11\n");
      }

      cudaStat1 = cudaMalloc((void**)&w_bar_device,n_ham*sizeof(cuDoubleComplex));
      if(cudaStat1 != cudaSuccess)
      {
         CLEANUP("Device malloc failed loc 11\n");
      }

      // w_bar = w * eig_mat(:,1);
      cubl_status = cublasZgemv(cubl_handle, transaction, n_ham, m, scalar1_device, w_device, n_ham, eigen_vec_device, 1, scalar2_device, w_bar_device, 1); 
      if(cubl_status != CUBLAS_STATUS_SUCCESS)
      {
         CLEANUP("cublasZgemv failed loc 12\n");
      }

        
//          cudaStat1 = cudaMemcpy(w_bar, w_bar_device,(size_t)(sizeof(cuDoubleComplex)*n_ham), cudaMemcpyDeviceToHost);
//          if(cudaStat1 != cudaSuccess)
//           {
//           CLEANUP("cudaMemcpy failed\n");
//           }
//          printf("w_bar[0].x=%f w_bar[0].y=%f\n", w_bar[0].x, w_bar[0].y); 
//          printf("w_bar[1].x=%f w_bar[1].y=%f\n", w_bar[1].x, w_bar[1].y);
//          printf("w_bar[2].x=%f w_bar[2].y=%f\n", w_bar[2].x, w_bar[2].y);
//          printf("w_bar[3].x=%f w_bar[3].y=%f\n", w_bar[3].x, w_bar[3].y);

      // r = (w_bar / s_u) - (s_v * u);
      vct_div_slr_minus_scl_vct_kernel_1_GPU<<<numBlock, threadPerBlock>>>(w_bar_device, u_device, r_device, n_ham, *s_u, *s_v);

      cudaFree(w_bar_device);

//cudaStat1=cudaDeviceSynchronize();
//if(cudaStat1 != cudaSuccess){
//CLEANUP("device sync failed 7\n");
//}


//          printf("i m here 1st r_cpu\n");
//          cudaStat1 = cudaMemcpy(r_cpu, r_device,(size_t)(n_ham*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);
//          printf("r_cpu[0].x=%f r_cpu[0].y=%f\n", r_cpu[0].x, r_cpu[0].y); 
//          printf("r_cpu[1].x=%f r_cpu[1].y=%f\n", r_cpu[1].x, r_cpu[1].y);
//          printf("r_cpu[2].x=%f r_cpu[2].y=%f\n", r_cpu[2].x, r_cpu[2].y);
//          printf("r_cpu[3].x=%f r_cpu[3].y=%f\n", r_cpu[3].x, r_cpu[3].y);


      // norm(r_bar);
      cubl_status = cublasDznrm2(cubl_handle, n_ham, r_device, 1, norm_r_device);

      if(cubl_status != CUBLAS_STATUS_SUCCESS){
      //printf("%d\n",cubl_status);
      CLEANUP("CUBLAS Dznrm2 failed loc 13\n");
      }

//cudaStat1=cudaDeviceSynchronize();
//if(cudaStat1 != cudaSuccess){
//CLEANUP("device sync failed 1\n");
//}

// printf("i m here 1st r_cpu\n");
// cudaStat1 = cudaMemcpy(r_cpu, r_device,(size_t)(n_ham*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);
// printf("r_cpu[0].x=%f r_cpu[0].y=%f\n", r_cpu[0].x, r_cpu[0].y); 
// printf("r_cpu[1].x=%f r_cpu[1].y=%f\n", r_cpu[1].x, r_cpu[1].y);
// printf("r_cpu[2].x=%f r_cpu[2].y=%f\n", r_cpu[2].x, r_cpu[2].y);
// printf("r_cpu[3].x=%f r_cpu[3].y=%f\n", r_cpu[3].x, r_cpu[3].y);



      cudaStat1 = cudaMemcpy(norm_r, norm_r_device,(size_t)(sizeof(double)), cudaMemcpyDeviceToHost);
      if(cudaStat1 != cudaSuccess)
          {
          CLEANUP("cudaMemcpy failed loc 14\n");
          }

         


//printf("norm_r %f\n", *norm_r);

// if(m == 4){
// cudaDeviceReset();
// exit(0);
// }

      //while(norm(r) < tol)
      while(*norm_r < *JD_tol)
      {

//printf("inside harmonic part\n");        

        k = k + 1;

        // X_bar = [X_bar, u];
        cudaStat1 = cudaMemcpy(u, u_device,(size_t)(sizeof(cuDoubleComplex)*n_ham), cudaMemcpyDeviceToHost);
        memcpy ( &X_bar[(k-1)*n_ham], u, sizeof(cuDoubleComplex)*n_ham );

        // lambda(k) = s_v + shift
        lambda[k-1] = *s_v + shift;
//         printf("lambda[%d] = %f \n", k, lambda[k-1]);
//         printf("JD count = %d\n", count);

        if(k == num_ev)
        {
          double temp;
            int index[k], temp_index;
            for(int i =0; i < k; i++)
            {
              index[i] = i;
            }
            
            if(band_type == 1)
            {
              for(int i =0; i < k-1; i++)
              {
                for(int j = 0; j < k-1-i; j++)
                {
                  if(lambda[j] > lambda[j+1])
                  {
                    temp = lambda[j];
                    lambda[j] = lambda[j+1];
                    lambda[j+1] = temp;
                    temp_index = index[j];
                    index[j] = index[j+1];
                    index[j+1] = temp_index;
                  }
                }
              }
            }
            
            if(band_type == 2)
            {
              for(int i =0; i < k-1; i++)
              {
                for(int j = 0; j < k-1-i; j++)
                {
                  if(lambda[j] < lambda[j+1])
                  {
                    temp = lambda[j];
                    lambda[j] = lambda[j+1];
                    lambda[j+1] = temp;
                    temp_index = index[j];
                    index[j] = index[j+1];
                    index[j+1] = temp_index;
                  }
                }
              }
            }
            
            // sort eigen_vector
            //double complex temp_vec;
            
            //temp_vec = (double complex *) malloc(sizeof(double complex)*n_ham);
         
            for(int i =0; i < k; i++){ memcpy(&eigen_vec_out[i*n_ham], &X_bar[index[i]*n_ham], sizeof(cuDoubleComplex)*n_ham); }

            memcpy(lambda_out, lambda, k*sizeof(double));

            for(int i =0; i < k; i++){
              printf("lambda[%d] = %f \n", i, lambda[i]);
            }
//           for(int i =0; i < k*n_ham; i++){ 
//           printf("X_bar[%d].x=%f X_bar[%d].y=%f\n", i, X_bar[i].x, i, X_bar[i].y); }
          //cudaProfilerStop();
          //cudaDeviceReset();
          printf("Total JD count = %d\n", count);
          printf("Total MxV = %d\n", (*ls_counter+count+(*ls_counter/ls_restart))); 
          return; // call cleanup
        }

        m = m -1;

        // no need of keeping M bcoz its calculated on every run
        //M = [];
        //cudaFree(M_device);
        //free(M);

        //cudaStat1 = cudaMalloc((void**)&M_device,m*m*sizeof(cuDoubleComplex));
        //cudaStat2 = cudaMemset(M_device, 0, m*m*sizeof(cuDoubleComplex));	

        // M is kept only on CPU here
        //M = (cuDoubleComplex *) calloc(m*m, sizeof(cuDoubleComplex));
        cudaStat1 = cudaMalloc((void**)&temp_device,n_ham*jd_max_step*sizeof(cuDoubleComplex));
        if(cudaStat1 != cudaSuccess)
        {
           CLEANUP("Device malloc failed loc 14.1\n");
        }

        for(int i = 0; i < m; i++)
        {
          //v_temp(:,i) = v * eig_mat(:,i+1);
          cubl_status = cublasZgemv(cubl_handle, CUBLAS_OP_N, n_ham, m+1, scalar1_device, v_device, n_ham, &eigen_vec_device[(i+1)*(m+1)], 1, scalar2_device, &temp_device[i*n_ham], 1);
if(cubl_status != CUBLAS_STATUS_SUCCESS)
{
 CLEANUP("cublasZgemv failed loc 15\n");
}

         }

        // v = v_temp(:,(1:m));
        cudaStat1 = cudaMemcpy(v_device, temp_device,(size_t)(sizeof(cuDoubleComplex)*n_ham*m), cudaMemcpyDeviceToDevice);

        for(int i = 0; i < m; i++)
        {

          // w_temp(:,i) = w * eig_mat(:,i+1);
          cubl_status = cublasZgemv(cubl_handle, CUBLAS_OP_N, n_ham, m+1, scalar1_device, w_device, n_ham, &eigen_vec_device[(i+1)*(m+1)], 1, scalar2_device, &temp_device[i*n_ham], 1);
if(cubl_status != CUBLAS_STATUS_SUCCESS)
{
 CLEANUP("cublasZgemv failed loc 16\n");
}

         }

         // w = w_temp(:,(1:m)); 
         cudaStat2 = cudaMemcpy(w_device, temp_device,(size_t)(sizeof(cuDoubleComplex)*n_ham*m), cudaMemcpyDeviceToDevice);

         cudaFree(temp_device);

          // M(i,i) = eig_val(i+1);
          //M[(i*m)+i].x = eigen_val[i+1];

         for(int i = 0; i < m; i++)
         {
          // eig_val(i) = eig_val(i+1);
          eigen_val[i] = eigen_val[i+1];

         }

         
        cudaFree(eigen_vec_device);
        free(eigen_vec);

        cudaStat1 = cudaMalloc((void**)&eigen_vec_device,m*m*sizeof(cuDoubleComplex));
        cudaStat2 = cudaMemset(eigen_vec_device, 0, m*m*sizeof(cuDoubleComplex));	

        eigen_vec = (cuDoubleComplex *) calloc(m*m, sizeof(cuDoubleComplex));
         
        for(int i = 0; i < m; i++)
        {
          // eig_mat(:,i) = e(:,i);
          eigen_vec[(i*m)+i].x = 1.0000; eigen_vec[(i*m)+i].y = 0.0000;
        }
        
        cudaStat1 = cudaMemcpy(eigen_vec_device, eigen_vec,(size_t)(sizeof(cuDoubleComplex)*m*m), cudaMemcpyHostToDevice);

        // s_u = norm(v(:,1));
        cubl_status = cublasDznrm2(cubl_handle, n_ham, v_device, 1, s_u_device);
        cudaStat1 = cudaMemcpy(s_u, s_u_device,(size_t)(sizeof(double)), cudaMemcpyDeviceToHost);

        // s_v = eig_val(1) / (s_u * s_u); 
        *s_v = eigen_val[0]/((*s_u) * (*s_u));
        cudaStat1 = cudaMemcpy(s_v_device, s_v,(size_t)(sizeof(double)), cudaMemcpyHostToDevice);
        
        // u = v(:,1)/s_u;
        vct_div_slr_kernel_1_GPU<<<numBlock, threadPerBlock>>>(v_device, u_device, n_ham, *s_u);

        // r = (w(:,1) / s_u) - (s_v * u); 
        vct_div_slr_minus_scl_vct_kernel_1_GPU<<<numBlock, threadPerBlock>>>(w_device, u_device, r_device, n_ham, *s_u, *s_v);

        // norm(r_bar);
        cubl_status = cublasDznrm2(cubl_handle, n_ham, r_device, 1, norm_r_device);

        cudaStat1 = cudaMemcpy(norm_r, norm_r_device,(size_t)(sizeof(double)), cudaMemcpyDeviceToHost);

      }// end while(norm(r) < tol)

       
      // restart
      if (m >= jd_max_step) 
      {

//printf("inside restart part\n");
        // no need of keeping M bcoz its calculated on every run
        // M = [];
        //cudaFree(M_device);
        //free(M);

        //cudaStat1 = cudaMalloc((void**)&M_device,jd_min_step*jd_min_step*sizeof(cuDoubleComplex));
        //cudaStat2 = cudaMemset(M_device, 0, jd_min_step*jd_min_step*sizeof(cuDoubleComplex));

        // M is kept only on CPU here
        //M = (cuDoubleComplex *) calloc(jd_min_step*jd_min_step, sizeof(cuDoubleComplex));

        cudaStat1 = cudaMalloc((void**)&temp_device,n_ham*jd_max_step*sizeof(cuDoubleComplex));
        if(cudaStat1 != cudaSuccess)
        {
           CLEANUP("Device malloc failed loc 14.1\n");
        }

        // for i = 1 : m_min
        for(int i = 0; i < jd_min_step; i++)
        {
          //v_temp(:,i) = v * eig_mat(:,i);
          cubl_status = cublasZgemv(cubl_handle, CUBLAS_OP_N, n_ham, m, scalar1_device, v_device, n_ham, &eigen_vec_device[i*m], 1, scalar2_device, &temp_device[i*n_ham], 1);
if(cubl_status != CUBLAS_STATUS_SUCCESS)
{
 CLEANUP("cublasZgemv failed loc 17\n");
}

         }

       // v = v_temp(:,(1:m_min));
       cudaStat1 = cudaMemcpy(v_device, temp_device,(size_t)(n_ham*jd_min_step*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToDevice);

        for(int i = 0; i < jd_min_step; i++)
        {

          // w_temp(:,i) = w * eig_mat(:,i);
          cubl_status = cublasZgemv(cubl_handle, CUBLAS_OP_N, n_ham, m, scalar1_device, w_device, n_ham, &eigen_vec_device[i*m], 1, scalar2_device, &temp_device[i*n_ham], 1);
if(cubl_status != CUBLAS_STATUS_SUCCESS)
{
 CLEANUP("cublasZgemv failed loc 18\n");
}

          // M(i,i) = eig_val(i);
          //M[(i*jd_min_step)+i].x = eigen_val[i];

         }

         // w = w_temp(:,(1:m_min)); 
         cudaStat2 = cudaMemcpy(w_device, temp_device,(size_t)(n_ham*jd_min_step*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToDevice);

//cudaStat1=cudaDeviceSynchronize();
//if(cudaStat1 != cudaSuccess){
//CLEANUP("device sync failed 2\n");
//}

       cudaFree(temp_device);

       m = jd_min_step;

       
       
// cudaStat1 = cudaMemcpy(v, v_device,(size_t)(n_ham*m*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);
// cudaStat2 = cudaMemcpy(w, w_device,(size_t)(n_ham*m*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);
// if(cudaStat1 != cudaSuccess){
// CLEANUP("cudaMemcpy failed\n");}
//  
// for(int i =0; i < m*n_ham; i++){ 
// printf("v[%d].x=%f v[%d].y=%f\n", i, v[i].x, i, v[i].y); }
// for(int i =0; i < m*n_ham; i++){ 
// printf("w[%d].x=%f w[%d].y=%f\n", i, w[i].x, i, w[i].y); }


      } // end if(restart)

// printf("i m here 2nd r_cpu\n");
// cudaStat1 = cudaMemcpy(r_cpu, r_device,(size_t)(n_ham*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);
// printf("r_cpu[0].x=%f r_cpu[0].y=%f\n", r_cpu[0].x, r_cpu[0].y); 
// printf("r_cpu[1].x=%f r_cpu[1].y=%f\n", r_cpu[1].x, r_cpu[1].y);
// printf("r_cpu[2].x=%f r_cpu[2].y=%f\n", r_cpu[2].x, r_cpu[2].y);
// printf("r_cpu[3].x=%f r_cpu[3].y=%f\n", r_cpu[3].x, r_cpu[3].y);

      cudaStat1 = cudaMemcpy(v, v_device,(size_t)(sizeof(cuDoubleComplex)*n_ham*m), cudaMemcpyDeviceToHost);
      cudaStat1 = cudaMemcpy(w, w_device,(size_t)(sizeof(cuDoubleComplex)*n_ham*m), cudaMemcpyDeviceToHost);


      cudaFree(v_device);
      cudaFree(w_device);
      
      cudaFree(M_device);
      cudaFree(eigen_vec_device);
      cudaFree(eigen_val_device);
      free(M);
      free(eigen_val);
      free(eigen_vec);


      cubl_status = cublasDznrm2(cubl_handle, n_ham, r_device, 1, norm_r_device);
      cudaStat1   = cudaMemcpy(norm_r, norm_r_device,(size_t)(sizeof(double)), cudaMemcpyDeviceToHost);
      if(cudaStat1 != cudaSuccess)
      {
         CLEANUP("cudaMemcpy failed loc 19\n");
      }

//printf("norm_r %f\n", *norm_r);

//printf("i am going towards gmres\n");
      
      if(*norm_r > 0.1 && flag ==1)
      {
        teta = shift;
      }
      else
      {
        teta = *s_v+shift;
        flag =0;
      }
      
// for(int i =0; i < n_ham; i++){
// printf("u[%d].x = %f u[%d].y = %f\n", i, u[i].x, i, u[i].y );
// }

      //cudaDeviceSynchronize(); 
      cudaStat1 = cudaMemcpy(u, u_device,(size_t)(sizeof(cuDoubleComplex)*n_ham), cudaMemcpyDeviceToHost);
      if(cudaStat1 != cudaSuccess)
      {
         CLEANUP("cudaMemcpy failed loc 20\n");
      }

      cudaFree(u_device);

// for(int i =0; i < n_ham; i++){
// printf("u[%d].x = %f u[%d].y = %f\n", i, u[i].x, i, u[i].y );
// }
// 
// printf("i am going towards gmres i m here\n");

      if(k !=0)
      {
//printf("i am going towards gmres i m here in if\n");
         Q_bar = (cuDoubleComplex *) malloc(sizeof(cuDoubleComplex)*n_ham*(k+1));
         cudaStat1 = cudaMalloc((void**)&Q_bar_device,sizeof(cuDoubleComplex)*n_ham*(k+1)); 
         memcpy (Q_bar, X_bar, sizeof(cuDoubleComplex)*n_ham*k);
         memcpy (&Q_bar[k*n_ham], u, sizeof(cuDoubleComplex)*n_ham);
         cudaStat2 = cudaMemcpy(Q_bar_device, Q_bar,(size_t)(sizeof(cuDoubleComplex)*n_ham*(k+1)), cudaMemcpyHostToDevice);
      }
      else
      {
//printf("i am going towards gmres i m here in else\n");
        Q_bar = (cuDoubleComplex *) malloc(n_ham*sizeof(cuDoubleComplex));
        cudaStat1 = cudaMalloc((void**)&Q_bar_device,n_ham*sizeof(cuDoubleComplex)); 
        memcpy (Q_bar, u, n_ham*sizeof(cuDoubleComplex));
        cudaStat2 = cudaMemcpy(Q_bar_device, Q_bar,(size_t)(n_ham*sizeof(cuDoubleComplex)), cudaMemcpyHostToDevice);
      }


//printf("form Q_bar\n");

     // A = (A - shift * I)
     tol_shift = shift-teta;
     shift_A_kernel<<<numBlock, threadPerBlock>>>(valptr_device_real, rowptr_device_real, colptr_device_real, n_ham, tol_shift);

//printf("calling gmres\n");

     gmres(valptr_device_real, rowptr_device_real, colptr_device_real, valptr_device_img, rowptr_device_img, colptr_device_img, n_ham, size_mat_real, r_device, ls_tol, ls_restart, ls_maxit, Q_bar_device, t_device, numBlock, threadPerBlock, coop, repeat, numBlocksMul, k, cusp_handle, cusp_descra, cubl_handle, ls_counter);


// cudaStat1 = cudaMemcpy(t_cpu, t_device,(size_t)(sizeof(cuDoubleComplex)*n_ham), cudaMemcpyDeviceToHost);
// printf("t_cpu[0].x=%f t_cpu[0].y=%f\n", t_cpu[0].x, t_cpu[0].y); 
// printf("t_cpu[1].x=%f t_cpu[1].y=%f\n", t_cpu[1].x, t_cpu[1].y);
// printf("t_cpu[2].x=%f t_cpu[2].y=%f\n", t_cpu[2].x, t_cpu[2].y);
// printf("t_cpu[3].x=%f t_cpu[3].y=%f\n", t_cpu[3].x, t_cpu[3].y);

////////////delete this/////////////////////////////////////////////////////////////////////////////////////
//printf("gmres done\n");
//cubl_status = cublasDznrm2(cubl_handle, n_ham, t_device, 1, norm_t_in_device);
//cudaStat1   = cudaMemcpy(norm_t_in, norm_t_in_device,(size_t)(sizeof(double)), cudaMemcpyDeviceToHost);
//printf("norm t after gmres = %f\n", *norm_t_in);
////////////////////////////////////////////////////////////////////////////////////////////////////////////



     free(Q_bar);
     cudaFree(Q_bar_device);

     if(k!=0) // MGS after linear solver
     {
       cudaStat1 = cudaMalloc((void**)&X_bar_device,n_ham*k*sizeof(cuDoubleComplex));
       if(cudaStat1 != cudaSuccess)
       {
         CLEANUP("cudaMemcpy device malloc falied loc 21\n");
       }
       cubl_status = cublasDznrm2(cubl_handle, n_ham, t_device, 1, norm_t_in_device);
       cudaStat1   = cudaMemcpy(norm_t_in, norm_t_in_device,(size_t)(sizeof(double)), cudaMemcpyDeviceToHost);
       cudaStat1   = cudaMemcpy(X_bar_device, X_bar,(size_t)(sizeof(cuDoubleComplex)*n_ham*k), cudaMemcpyHostToDevice);

       for(int i = 0; i < k; i++)
       {
         // y = X_bar(:,i)' * t;
         cubl_status = cublasZdotc(cubl_handle, n_ham, &X_bar_device[i*n_ham], 1, t_device, 1, y_device); 

         // t = t - y * X_bar(:,i);
         vct1_sub_mul_vct_kernel_1_GPU<<<numBlock, threadPerBlock>>>(t_device, &X_bar_device[i*n_ham], n_ham, y_device);

       }

        // norm(t);
        cubl_status = cublasDznrm2(cubl_handle, n_ham, t_device, 1, norm_t_device);

        cudaStat1 = cudaMemcpy(norm_t, norm_t_device,(size_t)(sizeof(double)), cudaMemcpyDeviceToHost); 

        if( ((*norm_t)/(*norm_t_in)) < 0.25)
        {
          for(int i = 0; i < k; i++)
          {
            // y = X_bar(:,i)' * t;
            cubl_status = cublasZdotc(cubl_handle, n_ham, &X_bar_device[i*n_ham], 1, t_device, 1, y_device); 

            // t = t - y * X_bar(:,i);
            vct1_sub_mul_vct_kernel_1_GPU<<<numBlock, threadPerBlock>>>(t_device, &X_bar_device[i*n_ham], n_ham, y_device);

           }
         }
         cudaFree(X_bar_device);
      }

      //free(eigen_vec);
      //free(Q_bar);
      //cudaFree(Q_bar_device);
count = count+1;
//printf("count = %d\n", count);
//if(count == 3)
//exit(0);


} // end while(k < num_ev)

printf("Total JD count = %d\n", count);
printf("Total GMRES iterations = %d\n", *ls_counter); 

//cudaProfilerStop();

//cudaDeviceReset();

}// end of JD



















void gmres(float *valptr_device_real, int *rowptr_device_real, int *colptr_device_real, float *valptr_device_img, int *rowptr_device_img, int *colptr_device_img, int n_ham, int size_mat_real, cuDoubleComplex * r_device, double ls_tol, int restart, int maxit, cuDoubleComplex *Q_bar_device, cuDoubleComplex *t_device, int numBlock, int threadPerBlock, int coop, int repeat, int numBlocksMul, int k, cusparseHandle_t cusp_handle_ls, cusparseMatDescr_t cusp_descra_ls, cublasHandle_t cubl_handle_ls, int *ls_counter)
{

//printf("inside gmres\n");
  //cuDoubleComplex *r_device_ls;
  cuDoubleComplex *X0_device, *Xmin_device, *dot_device, *P, *P_device;
  cuDoubleComplex *rq_device, *h_min, *g, cudoublecomplex_temp1, cudoublecomplex_temp2, cudoublecomplex_temp, sqrt_cudoublecomplex_temp;
  cuDoubleComplex *temp_device, *temp2_device, *test, *temp1_device;
  double *beta_device, *norm_r_device_ls, *norm_r_ls, *beta;
  cuDoubleComplex *v_ls_device, *w_ls_device, *g_device, *g_min, *minimizer_device, *h_device, *h; 
  cuDoubleComplex *scalar1_device_ls, *scalar2_device_ls;
  cuDoubleComplex *scalar1_ls, *scalar2_ls;
  //X0 = (cuDoubleComplex *)
  lapack_int *jpvt;
  lapack_int error_code=0;
  cuDoubleComplex coso, sino, scalar_1, scalar_2;
  cudaError_t cudaStat1,cudaStat2,cudaStat3,cudaStat4,cudaStat5,cudaStat6,cudaStat7,cudaStat8,cudaStat9,cudaStat10,cudaStat11,cudaStat12,cudaStat13,cudaStat14,cudaStat15,cudaStat16,cudaStat17,cudaStat18;

  scalar_1.x = 1.00;
  scalar_1.y = 0.00;
  scalar_2.x = 0.00;
  scalar_2.y = 0.00;

  scalar1_ls = (cuDoubleComplex *) malloc(sizeof(cuDoubleComplex));
  scalar2_ls = (cuDoubleComplex *) malloc(sizeof(cuDoubleComplex));

  scalar1_ls[0].x = 1.00;
  scalar1_ls[0].y = 0.00;
  scalar2_ls[0].x = 0.00;
  scalar2_ls[0].y = 0.00;

  const cuDoubleComplex jcmpx=make_cuDoubleComplex(0.0,1.0);

  cusparseStatus_t cusp_status_ls;
  //cusparseHandle_t cusp_handle_ls=0;
  //cusparseMatDescr_t cusp_descra_ls=0;
  cublasStatus_t cubl_status_ls;
  //cublasHandle_t cubl_handle_ls=0;

//   // initalization of CUBLAS library
//   cubl_status_ls = cublasCreate(&cubl_handle_ls);
//   if(cubl_status_ls != CUBLAS_STATUS_SUCCESS)
//   {
//     CLEANUP_LS("CUBLAS initilization failed \n");
//   }
// 
//   cubl_status_ls = cublasSetPointerMode(cubl_handle_ls, CUBLAS_POINTER_MODE_DEVICE);
//   if(cubl_status_ls != CUBLAS_STATUS_SUCCESS)
//   {
//     CLEANUP_LS("CUBLAS setting device pointer mode failed\n");
//   }
// 
//   /* initialize cusparse library */
//   cusp_status_ls= cusparseCreate(&cusp_handle_ls);
//   if (cusp_status_ls != CUSPARSE_STATUS_SUCCESS) 
//   {
//     CLEANUP_LS("CUSPARSE Library initialization failed\n");
//   }
// 
//  /* create and setup matrix descriptor */
//  cusp_status_ls= cusparseCreateMatDescr(&cusp_descra_ls);
//  if (cusp_status_ls != CUSPARSE_STATUS_SUCCESS) 
//  {
//    CLEANUP_LS("Matrix descriptor initialization failed 1");
//  }
// 
//  cusp_status_ls=cusparseSetMatType(cusp_descra_ls, CUSPARSE_MATRIX_TYPE_GENERAL);
//  if (cusp_status_ls != CUSPARSE_STATUS_SUCCESS) 
//  {
//    CLEANUP_LS("Matrix descriptor initialization failed 2");
//  }
//         
//  cusp_status_ls=cusparseSetMatIndexBase(cusp_descra_ls, CUSPARSE_INDEX_BASE_ONE);
//  if (cusp_status_ls != CUSPARSE_STATUS_SUCCESS) 
//  {
//    CLEANUP_LS("Matrix descriptor initialization failed 3");
//  }

//printf("inside gmres pt 0\n");

  cudaStat1 = cudaMalloc((void**)&X0_device,sizeof(cuDoubleComplex)*n_ham);
  cudaStat2 = cudaMalloc((void**)&Xmin_device,sizeof(cuDoubleComplex)*n_ham); 
  cudaStat3 = cudaMemset(X0_device, 0, n_ham*sizeof(cuDoubleComplex));
  cudaStat4 = cudaMemset(Xmin_device, 0, n_ham*sizeof(cuDoubleComplex));
  cudaStat5 = cudaMalloc((void**)&dot_device,sizeof(cuDoubleComplex)); 
  cudaStat6 = cudaMalloc((void**)&temp_device,sizeof(cuDoubleComplex)*n_ham);  
  cudaStat7 = cudaMalloc((void**)&temp2_device,sizeof(cuDoubleComplex)*n_ham);
  cudaStat8 = cudaMalloc((void**)&rq_device,sizeof(cuDoubleComplex)*n_ham);
  cudaStat9 = cudaMalloc((void**)&beta_device,sizeof(double));
  cudaStat10 = cudaMalloc((void**)&v_ls_device,sizeof(cuDoubleComplex)*n_ham*(restart+1));
  cudaStat11 = cudaMalloc((void**)&w_ls_device,sizeof(cuDoubleComplex)*n_ham*restart);
  cudaStat12 = cudaMalloc((void**)&h_device,sizeof(cuDoubleComplex)*(restart+1)*restart);
  cudaStat13 = cudaMalloc((void**)&g_device,sizeof(cuDoubleComplex)*(restart+1));
  cudaStat14 = cudaMalloc((void**)&minimizer_device,sizeof(cuDoubleComplex)*restart);
  cudaStat15 = cudaMalloc((void**)&norm_r_device_ls,sizeof(double));
  cudaStat16 = cudaMalloc((void**)&scalar1_device_ls,sizeof(cuDoubleComplex)); 
  cudaStat17 = cudaMalloc((void**)&scalar2_device_ls,sizeof(cuDoubleComplex));  
  cudaStat18 = cudaMalloc((void**)&temp1_device,n_ham*sizeof(cuDoubleComplex));  


  if((cudaStat1 != cudaSuccess)||
  (cudaStat2 != cudaSuccess) ||
  (cudaStat3 != cudaSuccess) ||
  (cudaStat4 != cudaSuccess) ||
  (cudaStat5 != cudaSuccess) ||
  (cudaStat6 != cudaSuccess) ||
  (cudaStat7 != cudaSuccess) ||
  (cudaStat8 != cudaSuccess) ||
  (cudaStat9 != cudaSuccess) ||
  (cudaStat10 != cudaSuccess) ||
  (cudaStat11 != cudaSuccess) ||
  (cudaStat12 != cudaSuccess) ||
  (cudaStat13 != cudaSuccess) ||
  (cudaStat14 != cudaSuccess) ||
  (cudaStat15 != cudaSuccess) ||
  (cudaStat16 != cudaSuccess) ||
  (cudaStat17 != cudaSuccess) ||
  (cudaStat18 != cudaSuccess))
  {
    CLEANUP_LS("cudaMemcpy failed loc 20\n");
  }


  // r_ls = -r
  vct1_neg_asg_vct2_kernel_1_GPU<<<numBlock, threadPerBlock>>>(r_device, r_device, n_ham);

// cuDoubleComplex *r_cpu;
// r_cpu = (cuDoubleComplex *) malloc(n_ham*sizeof(cuDoubleComplex));
// cudaMemcpy(r_cpu, r_device_ls,(size_t)(n_ham*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);
// printf("r_cpu[0].x=%f r_cpu[0].y=%f\n", r_cpu[0].x, r_cpu[0].y); 
// printf("r_cpu[1].x=%f r_cpu[1].y=%f\n", r_cpu[1].x, r_cpu[1].y);
// printf("r_cpu[2].x=%f r_cpu[2].y=%f\n", r_cpu[2].x, r_cpu[2].y);
// printf("r_cpu[3].x=%f r_cpu[3].y=%f\n", r_cpu[3].x, r_cpu[3].y);


  double norm_min =10000000000;
  double rcond;
  lapack_int *eff_rank_h_min;

  norm_r_ls = (double *) malloc(sizeof(double));
  beta = (double *) malloc(sizeof(double));
  h = (cuDoubleComplex *) malloc(sizeof(cuDoubleComplex)*(restart+1)*restart);
  h_min = (cuDoubleComplex *) malloc(restart*restart*sizeof(cuDoubleComplex));
  g = (cuDoubleComplex *) malloc((restart+1)*sizeof(cuDoubleComplex));
  g_min = (cuDoubleComplex *) malloc(restart*sizeof(cuDoubleComplex));
  test = (cuDoubleComplex *) malloc(sizeof(cuDoubleComplex));
  jpvt = (lapack_int *) malloc(restart*sizeof(lapack_int));
  //norm_r = (double *) malloc(sizeof(double));
  eff_rank_h_min = (lapack_int *) malloc(sizeof(lapack_int));

  cudaStat1 = cudaMemcpy(scalar1_device_ls, scalar1_ls,(size_t)(sizeof(cuDoubleComplex)), cudaMemcpyHostToDevice);
  cudaStat2 = cudaMemcpy(scalar2_device_ls, scalar2_ls,(size_t)(sizeof(cuDoubleComplex)), cudaMemcpyHostToDevice);

  if(maxit > n_ham)
  {
      maxit = n_ham;
  }

  if (restart > n_ham)
  {
      restart = n_ham;
  }

  //tol =tol * norm(b);
  cubl_status_ls = cublasDznrm2(cubl_handle_ls, n_ham, r_device, 1, norm_r_device_ls);
  cudaStat1 = cudaMemcpy(norm_r_ls, norm_r_device_ls,(size_t)(sizeof(double)), cudaMemcpyDeviceToHost); 

  ls_tol = (*norm_r_ls) * ls_tol;

  int restart_count = 0;

//printf("inside gmres pt 1\n");

// double * norm_Q_bar, *norm_Q_bar_device_ls;
// cudaMalloc((void**)&norm_Q_bar_device_ls,sizeof(double));
// norm_Q_bar = (double *) malloc(sizeof(double));
  
  while(restart_count < maxit)
  {
    // temp = x0;
    cudaStat1 = cudaMemcpy(temp_device, X0_device,(size_t)(n_ham*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToDevice);

// cuDoubleComplex *X0_cpu;
// X0_cpu = (cuDoubleComplex *) malloc(n_ham*sizeof(cuDoubleComplex));
// cudaStat1 = cudaMemcpy(X0_cpu, X0_device,(size_t)(n_ham*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);
// printf("X0_cpu[0].x=%f X0_cpu[0].y=%f\n", X0_cpu[0].x, X0_cpu[0].y); 
// printf("X0_cpu[1].x=%f X0_cpu[1].y=%f\n", X0_cpu[1].x, X0_cpu[1].y);
// printf("X0_cpu[2].x=%f X0_cpu[2].y=%f\n", X0_cpu[2].x, X0_cpu[2].y);
// printf("X0_cpu[3].x=%f X0_cpu[3].y=%f\n", X0_cpu[3].x, X0_cpu[3].y);

    for(int i = 0; i < k+1; i++)
    {
// cubl_status_ls = cublasDznrm2(cubl_handle_ls, n_ham, &Q_bar_device[i*n_ham], 1, norm_Q_bar_device_ls);
// cudaStat1 = cudaMemcpy(norm_Q_bar, norm_Q_bar_device_ls,(size_t)(sizeof(double)), cudaMemcpyDeviceToHost);    
// printf("norm_Q_bar = %f\n", *norm_Q_bar);    
       // temp = temp - (dot(Q_bar(:,i),x0)) * Q_bar(:,i); norm(Q_bar) == 1
       //cubl_status_ls = cublasDznrm2(cubl_handle_ls, n_ham, Q_bar_device[i*n_ham], 1, norm_Q_bar_device);
       cubl_status_ls = cublasZdotc(cubl_handle_ls, n_ham, &Q_bar_device[i*n_ham], 1, X0_device, 1, dot_device);
       vct1_sub_mul_vct_kernel_1_GPU<<<numBlock, threadPerBlock>>>(temp_device, &Q_bar_device[i*n_ham], n_ham, dot_device);
    }

//cuDoubleComplex *temp_cpu;
//temp_cpu = (cuDoubleComplex *) malloc(n_ham*sizeof(cuDoubleComplex));
// cudaStat1 = cudaMemcpy(temp_cpu, temp_device,(size_t)(n_ham*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);
// printf("temp_cpu[0].x=%f temp_cpu[0].y=%f\n", temp_cpu[0].x, temp_cpu[0].y); 
// printf("temp_cpu[1].x=%f temp_cpu[1].y=%f\n", temp_cpu[1].x, temp_cpu[1].y);
// printf("temp_cpu[2].x=%f temp_cpu[2].y=%f\n", temp_cpu[2].x, temp_cpu[2].y);
// printf("temp_cpu[3].x=%f temp_cpu[3].y=%f\n", temp_cpu[3].x, temp_cpu[3].y);

    // temp2 = A * temp;
    //cusp_status_ls= cusparseZcsrmv(cusp_handle_ls,CUSPARSE_OPERATION_NON_TRANSPOSE, n_ham, n_ham, size_mat_real, &scalar_1, cusp_descra_ls, valptr_device, rowptr_device, colptr_device, temp_device, &scalar_2, temp2_device);


     spmv_csr_hybrid_kernel_1_GPU<<<numBlocksMul, BLOCK_SIZE_MUL>>>(n_ham, rowptr_device_real, colptr_device_real, valptr_device_real, temp_device, temp2_device, repeat, coop, 0);

     spmv_csr_hybrid_kernel_1_GPU<<<numBlock, BLOCK_SIZE>>>(n_ham, rowptr_device_img, colptr_device_img, valptr_device_img, temp_device, temp1_device, 0);

     vct_pls_scl_mul_vct_kernel_1_GPU<<<numBlock, threadPerBlock>>>(temp2_device, temp1_device, jcmpx, n_ham);


//cuDoubleComplex *temp2_cpu;
//temp2_cpu = (cuDoubleComplex *) malloc(n_ham*sizeof(cuDoubleComplex));
// cudaStat1 = cudaMemcpy(temp2_cpu, temp2_device,(size_t)(n_ham*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);
// printf("temp2_cpu[0].x=%f temp2_cpu[0].y=%f\n", temp2_cpu[0].x, temp2_cpu[0].y); 
// printf("temp2_cpu[1].x=%f temp2_cpu[1].y=%f\n", temp2_cpu[1].x, temp2_cpu[1].y);
// printf("temp2_cpu[2].x=%f temp2_cpu[2].y=%f\n", temp2_cpu[2].x, temp2_cpu[2].y);
// printf("temp2_cpu[3].x=%f temp2_cpu[3].y=%f\n", temp2_cpu[3].x, temp2_cpu[3].y);

   // temp = temp2;
   cudaStat1 = cudaMemcpy(temp_device, temp2_device,(size_t)(n_ham*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToDevice);
   
//printf("inside gmres pt 2\n");    

   for(int i = 0; i < k+1; i++)
    {
       // temp = temp - (dot(u(:,t),temp2)) * u(:,t);
       cubl_status_ls = cublasZdotc(cubl_handle_ls, n_ham, &Q_bar_device[i*n_ham], 1, temp2_device, 1, dot_device);
       vct1_sub_mul_vct_kernel_1_GPU<<<numBlock, threadPerBlock>>>(temp_device, &Q_bar_device[i*n_ham], n_ham, dot_device);
    }

// cudaStat1 = cudaMemcpy(temp_cpu, temp_device,(size_t)(n_ham*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);
// printf("temp_cpu[0].x=%f temp_cpu[0].y=%f\n", temp_cpu[0].x, temp_cpu[0].y); 
// printf("temp_cpu[1].x=%f temp_cpu[1].y=%f\n", temp_cpu[1].x, temp_cpu[1].y);
// printf("temp_cpu[2].x=%f temp_cpu[2].y=%f\n", temp_cpu[2].x, temp_cpu[2].y);
// printf("temp_cpu[3].x=%f temp_cpu[3].y=%f\n", temp_cpu[3].x, temp_cpu[3].y);


    //rq = b - temp2;     here b = r_device_ls
    vct1_sub_vct2_asg_vct3_kernel_1_GPU<<<numBlock, threadPerBlock>>>(r_device, temp_device, rq_device, n_ham);

// cudaStat1 = cudaMemcpy(temp2_cpu, r_device_ls,(size_t)(n_ham*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);
// printf("temp2_cpu[0].x=%f temp2_cpu[0].y=%f\n", temp2_cpu[0].x, temp2_cpu[0].y); 
// printf("temp2_cpu[1].x=%f temp2_cpu[1].y=%f\n", temp2_cpu[1].x, temp2_cpu[1].y);
// printf("temp2_cpu[2].x=%f temp2_cpu[2].y=%f\n", temp2_cpu[2].x, temp2_cpu[2].y);
// printf("temp2_cpu[3].x=%f temp2_cpu[3].y=%f\n", temp2_cpu[3].x, temp2_cpu[3].y);

    
    // beta=norm(rq);
    cubl_status_ls = cublasDznrm2(cubl_handle_ls, n_ham, rq_device, 1, beta_device);

    cudaStat1 = cudaMemcpy(beta, beta_device,(size_t)(sizeof(double)), cudaMemcpyDeviceToHost);
//printf("inside gnres norm(rq) = %f \n", *beta);

    // v(:,1)=rq/beta; 
    vct_div_slr_kernel_1_GPU<<<numBlock, threadPerBlock>>>(rq_device, v_ls_device, n_ham, *beta);


    cudaStat1 = cudaMemset(h_device, 0, sizeof(cuDoubleComplex)*restart*(restart+1));

//printf("inside gmres pt 3\n"); 

    for(int j = 0; j < restart; j++)
    {
      // temp = v(:,j);
      cudaStat1 = cudaMemcpy(temp_device, &v_ls_device[j*n_ham],(size_t)(n_ham*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToDevice);

      for(int i = 0; i < k+1; i++)
      {
       // temp = temp - (dot(Q_bar(:,i),v(:,j))) * Q_bar(:,i); norm(Q_bar) == 1
       cubl_status_ls = cublasZdotc(cubl_handle_ls, n_ham, &Q_bar_device[i*n_ham], 1, &v_ls_device[j*n_ham], 1, dot_device);
       vct1_sub_mul_vct_kernel_1_GPU<<<numBlock, threadPerBlock>>>(temp_device, &Q_bar_device[i*n_ham], n_ham, dot_device);
      }

      // temp2=A*temp;
      // w(:,j)=temp2;
      //cusp_status_ls= cusparseZcsrmv(cusp_handle_ls,CUSPARSE_OPERATION_NON_TRANSPOSE, n_ham, n_ham, size_mat_real, &scalar_1, cusp_descra_ls, valptr_device, rowptr_device, colptr_device, temp_device, &scalar_2, temp2_device);

     spmv_csr_hybrid_kernel_1_GPU<<<numBlocksMul, BLOCK_SIZE_MUL>>>(n_ham, rowptr_device_real, colptr_device_real, valptr_device_real, temp_device, temp2_device, repeat, coop, 0);

     spmv_csr_hybrid_kernel_1_GPU<<<numBlock, BLOCK_SIZE>>>(n_ham, rowptr_device_img, colptr_device_img, valptr_device_img, temp_device, temp1_device, 0);

     vct_pls_scl_mul_vct_kernel_1_GPU<<<numBlock, threadPerBlock>>>(temp2_device, temp1_device, jcmpx, n_ham);

     cudaStat1 = cudaMemcpy(&w_ls_device[j*n_ham], temp2_device,(size_t)(n_ham*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToDevice);

//printf("inside gmres pt 4\n"); 

     for(int i = 0; i < k+1; i++)
      {
       // temp = temp - (dot(Q_bar(:,i),w(:,j))) * Q_bar(:,i); norm(Q_bar) == 1
       cubl_status_ls = cublasZdotc(cubl_handle_ls, n_ham, &Q_bar_device[i*n_ham], 1, &w_ls_device[j*n_ham], 1, dot_device);
       vct1_sub_mul_vct_kernel_1_GPU<<<numBlock, threadPerBlock>>>(temp2_device, &Q_bar_device[i*n_ham], n_ham, dot_device);
      }

//printf("inside gmres pt 4.1\n"); 

      // w(:,j)=temp2;
      cudaStat1 = cudaMemcpy(&w_ls_device[j*n_ham], temp2_device,(size_t)(n_ham*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToDevice);

// cudaStat1 = cudaMemcpy(temp2_cpu, temp2_device,(size_t)(n_ham*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);
// printf("temp2_cpu[0].x=%f temp2_cpu[0].y=%f\n", temp2_cpu[0].x, temp2_cpu[0].y); 
// printf("temp2_cpu[1].x=%f temp2_cpu[1].y=%f\n", temp2_cpu[1].x, temp2_cpu[1].y);
// printf("temp2_cpu[2].x=%f temp2_cpu[2].y=%f\n", temp2_cpu[2].x, temp2_cpu[2].y);
// printf("temp2_cpu[3].x=%f temp2_cpu[3].y=%f\n", temp2_cpu[3].x, temp2_cpu[3].y);


      for(int i = 0; i <= j ; i++)
      {
        // h(i,j)=dot(w(:,j),v(:,i));
        cubl_status_ls = cublasZdotc(cubl_handle_ls, n_ham, &w_ls_device[j*n_ham], 1, &v_ls_device[i*n_ham], 1, &h_device[(j*(restart+1))+i]);

//cudaStat1 = cudaMemcpy(test, &h_device[(j*(restart+1))+i],(size_t)(sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);
//printf("h(i,j).x =%f h(i,j).y = %f\n", test[0].x, test[0].y);


        //  w(:,j)=w(:,j)-h(i,j)*v(:,i);
        vct1_sub_mul_vct_kernel_1_GPU<<<numBlock, threadPerBlock>>>(&w_ls_device[j*n_ham], &v_ls_device[i*n_ham], n_ham, &h_device[(j*(restart+1))+i]);
      }

      // h(j+1,j)=norm(w(:,j));
      cubl_status_ls = cublasDznrm2(cubl_handle_ls, n_ham, &w_ls_device[j*n_ham], 1, &h_device[(j*(restart+1))+j+1].x);

      cudaStat1 = cudaMemcpy(test, &h_device[(j*(restart+1))+j+1],(size_t)(sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);

//printf("h(j+1,j) =%f \n", test[0].x);

      //if h(j+1,j)==0
      //if(h_device[((j+1)*(restart+1))+j].x == 0.0000)
      if(test[0].x == 0.0000)
        {
//printf("inside gmres i m here if\n");
         restart=j;
        }
      else
        {
//printf("inside gmres i m here else\n"); 
           //v(:,j+1)=w(:,j)/h(j+1,j);
           vct_1_div_asg_to_vct_2_kernel_1_GPU<<<numBlock, threadPerBlock>>>(&w_ls_device[j*n_ham], &v_ls_device[(j+1)*n_ham], n_ham, &h_device[(j*(restart+1))+j+1].x);
        }

// cuDoubleComplex *v_cpu;
// v_cpu = (cuDoubleComplex *) malloc(n_ham*sizeof(cuDoubleComplex));
// cudaStat1 = cudaMemcpy(v_cpu,  &v_ls_device[(j+1)*n_ham],(size_t)(n_ham*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);
// printf("v_cpu[0].x=%f v_cpu[0].y=%f\n", v_cpu[0].x, v_cpu[0].y); 
// printf("v_cpu[1].x=%f v_cpu[1].y=%f\n", v_cpu[1].x, v_cpu[1].y);
// printf("v_cpu[2].x=%f v_cpu[2].y=%f\n", v_cpu[2].x, v_cpu[2].y);
// printf("v_cpu[3].x=%f v_cpu[3].y=%f\n", v_cpu[3].x, v_cpu[3].y);

   }

//printf("inside gmres pt 4.4\n"); 

        // g(1:m+1,:)=0;
        cudaStat1 = cudaMemset(g_device, 0, sizeof(cuDoubleComplex)*(restart+1));
        

        // g(1,:)=beta;
        cudaStat2 = cudaMemcpy(beta, beta_device,(size_t)(sizeof(double)), cudaMemcpyDeviceToHost);
        test[0].x = *beta; test[0].y = 0.0000;
        cudaStat1 = cudaMemcpy(g_device, test,(size_t)(sizeof(cuDoubleComplex)), cudaMemcpyHostToDevice);
        //if(cudaStat1 != cudaSucess

//printf("beta = %f\n", *beta);

// cudaStat2 = cudaMemcpy(g, g_device,(size_t)(sizeof(cuDoubleComplex)*(restart+1)), cudaMemcpyDeviceToHost);
// for(int i = 0; i<restart+1; i++)
// { 
// printf(" i am here g[%d].x=%f g[%d].y=%f\n", i, g[i].x, i, g[i].y); 
// }



//printf("inside gmres pt 5\n"); 

  cudaStat1 = cudaMemcpy(h, h_device,(size_t)(sizeof(cuDoubleComplex)*restart*(restart+1)), cudaMemcpyDeviceToHost);

  for(int j = 0; j < restart; j++)
    {
     cudaStat1 = cudaMalloc((void**)&P_device,sizeof(cuDoubleComplex)*(restart+1)*(restart+1));
     cudaStat2 = cudaMemset(P_device, 0, sizeof(cuDoubleComplex)*(restart+1)*(restart+1));

     // eigen_vec is kept only on CPU here
     P = (cuDoubleComplex *) calloc((restart+1)*(restart+1), sizeof(cuDoubleComplex));

     for(int i = 0; i < restart+1; i++)
     {
       // P = eye(restart+1);
       P[(i*(restart+1))+i].x = 1.0000; P[(i*(restart+1))+i].y = 0.0000;
     }

// printf("i m here\n");
// for(int i = 0; i<(restart+1)*(restart+1); i++)
// { 
// printf(" P[%d].x=%f P[%d].y=%f ", i, P[i].x, i, P[i].y); 
// }
// 
// printf("\n");


//printf("inside gmres pt 6\n"); 

     // sqrt(h(j+1,j)^2 + h(j,j)^2);
     // h(j+1,j)^2
     cudoublecomplex_temp1.x = (h[(j*(restart+1))+j+1].x * h[(j*(restart+1))+j+1].x) - (h[(j*(restart+1))+j+1].y * h[(j*(restart+1))+j+1].y);
     cudoublecomplex_temp1.y = (h[(j*(restart+1))+j+1].x * h[(j*(restart+1))+j+1].y) + (h[(j*(restart+1))+j+1].y * h[(j*(restart+1))+j+1].x);
     //h(j,j)^2
     cudoublecomplex_temp2.x = (h[j*(restart+1)+j].x * h[j*(restart+1)+j].x) - (h[j*(restart+1)+j].y * h[j*(restart+1)+j].y);
     cudoublecomplex_temp2.y = (h[j*(restart+1)+j].x * h[j*(restart+1)+j].y) + (h[j*(restart+1)+j].y * h[j*(restart+1)+j].x);

     cudoublecomplex_temp.x = cudoublecomplex_temp1.x + cudoublecomplex_temp2.x;
     cudoublecomplex_temp.y = cudoublecomplex_temp1.y + cudoublecomplex_temp2.y;

//printf("cudoublecomplex_temp1.x = %f, cudoublecomplex_temp2.x = %f\n", cudoublecomplex_temp1.x,cudoublecomplex_temp2.x);
//printf("cudoublecomplex_temp1.y = %f, cudoublecomplex_temp2.y = %f\n", cudoublecomplex_temp1.y, cudoublecomplex_temp2.y);


     // sqrt()
     sqrt_cudoublecomplex_temp.x = sqrt((cudoublecomplex_temp.x + sqrt(cudoublecomplex_temp.x*cudoublecomplex_temp.x + cudoublecomplex_temp.y*cudoublecomplex_temp.y))/2);
     if(cudoublecomplex_temp.y != 0.00000000000)
     {
     sqrt_cudoublecomplex_temp.y = (cudoublecomplex_temp.y/abs(cudoublecomplex_temp.y)) * sqrt(((-cudoublecomplex_temp.x + sqrt(cudoublecomplex_temp.x*cudoublecomplex_temp.x + cudoublecomplex_temp.y*cudoublecomplex_temp.y))/2));
     }
     else
     {
     sqrt_cudoublecomplex_temp.y = 0.0000000000;
     }


//printf("sqrt_cudoublecomplex_temp.x = %f, sqrt_cudoublecomplex_temp.y = %f\n", sqrt_cudoublecomplex_temp.x, sqrt_cudoublecomplex_temp.y);
//printf("lest check this = %f\n",(-cudoublecomplex_temp.x + sqrt(cudoublecomplex_temp.x*cudoublecomplex_temp.x + cudoublecomplex_temp.y*cudoublecomplex_temp.y))/2);
  
     // sino=h(j+1,j)/(sqrt(h(j+1,j)^2 + h(j,j)^2));
     sino.x = ((h[(j*(restart+1))+j+1].x * sqrt_cudoublecomplex_temp.x) + (h[(j*(restart+1))+j+1].y * sqrt_cudoublecomplex_temp.y)) / ((sqrt_cudoublecomplex_temp.x * sqrt_cudoublecomplex_temp.x) + (sqrt_cudoublecomplex_temp.y * sqrt_cudoublecomplex_temp.y));

     sino.y = ((h[(j*(restart+1))+j+1].y * sqrt_cudoublecomplex_temp.x) - (h[(j*(restart+1))+j+1].x * sqrt_cudoublecomplex_temp.y)) / ((sqrt_cudoublecomplex_temp.x * sqrt_cudoublecomplex_temp.x) + (sqrt_cudoublecomplex_temp.y * sqrt_cudoublecomplex_temp.y));


     // coso=h(j,j)/(sqrt(h(j+1,j)^2 + h(j,j)^2));
     coso.x = ((h[j*(restart+1)+j].x * sqrt_cudoublecomplex_temp.x) + (h[j*(restart+1)+j].y * sqrt_cudoublecomplex_temp.y)) / ((sqrt_cudoublecomplex_temp.x * sqrt_cudoublecomplex_temp.x) + (sqrt_cudoublecomplex_temp.y * sqrt_cudoublecomplex_temp.y));

     coso.y = ((h[j*(restart+1)+j].y * sqrt_cudoublecomplex_temp.x) - (h[j*(restart+1)+j].x * sqrt_cudoublecomplex_temp.y)) / ((sqrt_cudoublecomplex_temp.x * sqrt_cudoublecomplex_temp.x) + (sqrt_cudoublecomplex_temp.y * sqrt_cudoublecomplex_temp.y));

//printf("coso.x = %f, coso.y = %f\n", coso.x, coso.y);
//printf("sino.x = %f, sino.y = %f\n", sino.x, sino.y);

     // P(j,j)=conj(coso);
     P[j*(restart+1)+j] = cuConj(coso);

     // P(j+1,j+1)=coso;
     P[(j+1)*(restart+1)+j+1].x = coso.x;  P[(j+1)*(restart+1)+j+1].y = coso.y;

     // P(j,j+1)=conj(sino);
     P[(j+1)*(restart+1)+j] = cuConj(sino);

     // P(j+1,j)=-sino;
     P[(j*(restart+1))+j+1].x = -sino.x;  P[(j*(restart+1))+j+1].y = -sino.y;

     cudaStat1 = cudaMemcpy(h_device, h,(size_t)(sizeof(cuDoubleComplex)*restart*(restart+1)), cudaMemcpyHostToDevice);
     cudaStat2 = cudaMemcpy(P_device, P,(size_t)(sizeof(cuDoubleComplex)*(restart+1)*(restart+1)), cudaMemcpyHostToDevice);


// for(int i = 0; i<(restart)*(restart+1); i++)
// { 
// printf("h[%d].x=%f h[%d].y=%f\n", i, h[i].x, i, h[i].y); 
// }
// // 
// for(int i = 0; i<(restart+1)*(restart+1); i++)
// { 
// printf(" P[%d].x=%f P[%d].y=%f ", i, P[i].x, i, P[i].y); 
// }
// 
// printf("\n");


     // h=P*h;
     cubl_status_ls = cublasZgemm(cubl_handle_ls, CUBLAS_OP_N, CUBLAS_OP_N, restart+1, restart, restart+1, scalar1_device_ls, P_device, restart+1, h_device, restart+1, scalar2_device_ls, h_device, restart+1);

    // g=P*g;
    cubl_status_ls = cublasZgemv(cubl_handle_ls, CUBLAS_OP_N, restart+1, restart+1, scalar1_device_ls, P_device, restart+1, g_device, 1, scalar2_device_ls, g_device, 1);

    // update h since h = P*h
    cudaStat1 = cudaMemcpy(h, h_device,(size_t)(sizeof(cuDoubleComplex)*restart*(restart+1)), cudaMemcpyDeviceToHost);

    free(P);
    cudaFree(P_device);
   

// cudaStat2 = cudaMemcpy(g, g_device,(size_t)(sizeof(cuDoubleComplex)*(restart+1)), cudaMemcpyDeviceToHost);
// for(int i = 0; i<restart+1; i++)
// { 
// printf("g[%d].x=%f g[%d].y=%f\n", i, g[i].x, i, g[i].y); 
// }
   }

//printf("inside gmres pt 7\n"); 

   //cudaStat1 = cudaMemcpy(h, h_device,(size_t)(sizeof(cuDoubleComplex)*(restart+1)*restart), cudaMemcpyDeviceToHost);
   cudaStat2 = cudaMemcpy(g, g_device,(size_t)(sizeof(cuDoubleComplex)*(restart+1)), cudaMemcpyDeviceToHost);

   // minimizer=h(1:m,1:m)\g(1:m,1);
   for(int i = 0; i < restart; i++) { memcpy (&h_min[restart*i], &h[(restart*i)+i], sizeof(cuDoubleComplex)*restart); }
   memcpy(g_min, g, sizeof(cuDoubleComplex)*restart);
   for(int i = 0; i < restart; i++)
   {
    jpvt[i] = 0;
   //memset(jpvt, 0, sizeof(lapack_int)*restart);
   }

// for(int i = 0; i<(restart*restart); i++)
// { 
// printf("h_min[%d].x=%f h_min[%d].y=%f\n", i, h_min[i].x, i, h_min[i].y); 
// }
// 
// for(int i = 0; i<restart+1; i++)
// { 
// printf("g[%d].x=%f g[%d].y=%f\n", i, g[i].x, i, g[i].y); 
// }
// 
// 
// for(int i = 0; i<restart; i++)
// { 
// printf("g_min[%d].x=%f g_min[%d].y=%f\n", i, g_min[i].x, i, g_min[i].y); 
// }
// 
// 
// 
// printf("inside gmres pt 7.1\n"); 

   // LAPACKE_zgelsy( int matrix_layout, lapack_int m, lapack_int n, lapack_int nrhs, lapack_complex_double* a, lapack_int lda, lapack_complex_double* b, lapack_int ldb, lapack_int* jpvt, double rcond, lapack_int* rank );
   error_code = LAPACKE_zgelsy(LAPACK_COL_MAJOR, restart, restart, 1, (MKL_Complex16 *)h_min, restart, (MKL_Complex16 *)g_min, restart, jpvt, rcond, eff_rank_h_min);

// printf("error_code lapack_Zgelsy =  %d\n", error_code);
// 
// printf("g_min[0].x=%f g_min[0].y=%f\n", g_min[0].x, g_min[0].y); 
// printf("g_min[1].x=%f g_min[1].y=%f\n", g_min[1].x, g_min[1].y);
// printf("g_min[2].x=%f g_min[2].y=%f\n", g_min[2].x, g_min[2].y);
// 


   // xm=x0+v(:,1:m)*minimizer;
   cudaStat1 = cudaMemcpy(minimizer_device, g_min,(size_t)(sizeof(cuDoubleComplex)*restart), cudaMemcpyHostToDevice);
   cubl_status_ls = cublasZgemv(cubl_handle_ls, CUBLAS_OP_N, n_ham, restart, scalar1_device_ls, v_ls_device, n_ham, minimizer_device, 1, scalar2_device_ls, temp_device, 1);
   vct1_add_vct2_asg_vct3_kernel_1_GPU<<<numBlock, threadPerBlock>>>(X0_device, temp_device, temp2_device, n_ham);



/*cudaMemcpy(temp2_cpu, temp2_device,(size_t)(sizeof(cuDoubleComplex)*n_ham), cudaMemcpyDeviceToHost);
printf("temp2_cpu[0].x=%f temp2_cpu[0].y=%f\n", temp2_cpu[0].x, temp2_cpu[0].y); 
printf("temp2_cpu[1].x=%f temp2_cpu[1].y=%f\n", temp2_cpu[1].x, temp2_cpu[1].y);
printf("temp2_cpu[2].x=%f temp2_cpu[2].y=%f\n", temp2_cpu[2].x, temp2_cpu[2].y);
printf("temp2_cpu[3].x=%f temp2_cpu[3].y=%f\n", temp2_cpu[3].x, temp2_cpu[3].y);*/

// printf("abs(g(m+1,1)) = %f\n", sqrt(g[restart].x*g[restart].x+g[restart].y*g[restart].y));
// if(isnan(sqrt(g[restart].x*g[restart].x+g[restart].y*g[restart].y)))
// exit(0);
   
   // if abs(g(m+1,1))<tol 
   if(sqrt(g[restart].x*g[restart].x+g[restart].y*g[restart].y)<ls_tol)
   {
     // x = xm;
     cudaStat1 = cudaMemcpy(t_device, temp2_device,(size_t)(sizeof(cuDoubleComplex)*n_ham), cudaMemcpyDeviceToDevice);

     *ls_counter =  *ls_counter + restart_count * restart;

     cudaFree(X0_device);
     cudaFree(Xmin_device);
     cudaFree(dot_device);
     cudaFree(temp_device);
     cudaFree(temp2_device);
     cudaFree(rq_device);
     cudaFree(beta_device);
     cudaFree(v_ls_device);
     cudaFree(w_ls_device);
     cudaFree(h_device);
     cudaFree(g_device);
     cudaFree(minimizer_device);
     cudaFree(norm_r_device_ls);
     cudaFree(scalar1_device_ls);
     cudaFree(scalar2_device_ls);
     cudaFree(temp1_device);

     free(norm_r_ls);
     free(beta);
     free(h);
     free(h_min);
     free(g);
     free(g_min);
     free(test);
     free(jpvt);
     free(eff_rank_h_min);

     cudaDeviceSynchronize();

     return;
   }
   else
   {
     // x0=xm;  
     cudaStat1 = cudaMemcpy(X0_device, temp2_device,(size_t)(sizeof(cuDoubleComplex)*n_ham), cudaMemcpyDeviceToDevice);
     // restart=restart+1;
     restart_count = restart_count + 1;
     // if abs(g(m+1,1)) <= normmin
     if(sqrt(g[restart].x*g[restart].x+g[restart].y*g[restart].y)<= norm_min)
     {
       // xmin = xm;
       cudaStat1 = cudaMemcpy(Xmin_device, temp2_device,(size_t)(sizeof(cuDoubleComplex)*n_ham), cudaMemcpyDeviceToDevice);
       norm_min = sqrt(g[restart].x*g[restart].x+g[restart].y*g[restart].y);
     }
   }

  } // end of while

 // x = xmin;
 cudaStat1 = cudaMemcpy(t_device, Xmin_device,(size_t)(sizeof(cuDoubleComplex)*n_ham), cudaMemcpyDeviceToDevice);

 *ls_counter =  *ls_counter + restart_count * restart;

 cudaFree(X0_device);
 cudaFree(Xmin_device);
 cudaFree(dot_device);
 cudaFree(temp_device);
 cudaFree(temp2_device);
 cudaFree(rq_device);
 cudaFree(beta_device);
 cudaFree(v_ls_device);
 cudaFree(w_ls_device);
 cudaFree(h_device);
 cudaFree(g_device);
 cudaFree(minimizer_device);
 cudaFree(norm_r_device_ls);
 cudaFree(scalar1_device_ls);
 cudaFree(scalar2_device_ls);
 cudaFree(temp1_device);

 free(norm_r_ls);
 free(beta);
 free(h);
 free(h_min);
 free(g);
 free(g_min);
 free(test);
 free(jpvt);
 free(eff_rank_h_min);

 cudaDeviceSynchronize();

 return;
}



























__global__ void vct1_add_vct2_asg_vct3_kernel_1_GPU(cuDoubleComplex *vct1, cuDoubleComplex * vct2, cuDoubleComplex * vct3, int n_ham)
{
  int ID;
  ID = blockDim.x * blockIdx.x + threadIdx.x;
  if(ID < n_ham)
  {
   vct3[ID].x = vct1[ID].x + vct2[ID].x;
   vct3[ID].y = vct1[ID].y + vct2[ID].y;
  }
}




__global__ void vct1_sub_vct2_asg_vct3_kernel_1_GPU(cuDoubleComplex *vct1, cuDoubleComplex * vct2, cuDoubleComplex * vct3, int n_ham)
{
  int ID;
  ID = blockDim.x * blockIdx.x + threadIdx.x;
  if(ID < n_ham)
  {
   vct3[ID].x = vct1[ID].x - vct2[ID].x;
   vct3[ID].y = vct1[ID].y - vct2[ID].y;
  }
}



__global__ void vct1_neg_asg_vct2_kernel_1_GPU(cuDoubleComplex *vct1, cuDoubleComplex * vct2, int n_ham)
{
  int ID;
  ID = blockDim.x * blockIdx.x + threadIdx.x;
  if(ID < n_ham)
  {
    vct1[ID].x =  -vct2[ID].x;
    vct1[ID].y =  -vct2[ID].y;
   }
}



__global__ void vct_div_slr_minus_scl_vct_kernel_1_GPU(cuDoubleComplex * scr1, cuDoubleComplex * scr2, cuDoubleComplex * des, int n_ham, double Slr1, double Slr2)
{
 int ID;
 ID = blockDim.x * blockIdx.x + threadIdx.x;
 if(ID < n_ham)
 {
   des[ID].x = (scr1[ID].x/Slr1) - (Slr2 * scr2[ID].x);
   des[ID].y = (scr1[ID].y/Slr1) - (Slr2 * scr2[ID].y);
 }
__syncthreads();
}



__global__ void vct_div_slr_kernel_1_GPU(cuDoubleComplex * scr, cuDoubleComplex * des, int n_ham, double Slr)
{
 int ID;
 ID = blockDim.x * blockIdx.x + threadIdx.x;
 if(ID < n_ham)
 {
   des[ID].x = scr[ID].x / Slr;
   des[ID].y = scr[ID].y / Slr;
 }
}


__global__ void vct_1_div_asg_to_vct_2_kernel_1_GPU(cuDoubleComplex * vct_1, cuDoubleComplex * vct_2, int n_ham, double *beta)
{
  int ID;
  ID = blockDim.x * blockIdx.x + threadIdx.x;
  if(ID < n_ham)
  {
    vct_2[ID].x =  vct_1[ID].x / beta[0];
    vct_2[ID].y =  vct_1[ID].y / beta[0];
  }
}








__global__ void vct1_sub_mul_vct_kernel_1_GPU(cuDoubleComplex * vct1, cuDoubleComplex * vct2, int n_ham, cuDoubleComplex * dot)
{
  int  ID;
  ID = blockDim.x * blockIdx.x + threadIdx.x;
  if(ID < n_ham)
  {
    vct1[ID].x = vct1[ID].x - ((dot[0].x * vct2[ID].x) - (dot[0].y * vct2[ID].y));
    vct1[ID].y = vct1[ID].y - ((dot[0].x * vct2[ID].y) + (dot[0].y * vct2[ID].x));
  }	
}







__global__ void vct_sub_scl_mul_vct_kernel_1_GPU(cuDoubleComplex *vct1, cuDoubleComplex * vct2, double scalar, int n_ham)
{
  int ID;
  ID = blockDim.x * blockIdx.x + threadIdx.x;
  if(ID < n_ham)
  {
    vct1[ID].x =  vct1[ID].x - (scalar * vct2[ID].x);
    vct1[ID].y =  vct1[ID].y - (scalar * vct2[ID].y);
   }
}





__global__ void shift_A_kernel(float * val, int * row, int * col, int n_ham, double shift)
{
  int ID;
  ID = blockDim.x * blockIdx.x + threadIdx.x;
  if(ID < n_ham)
  {
      for(int j=row[ID]-1; j <= row[ID+1]-1; j++)
      {
       if(col[j] == ID+1)
       {
        val[j] = val[j] + shift; 
        break;
       }
      }
  }
}




__global__ void cpy_vct_1_to_vct_2_kernel_1_GPU(cuDoubleComplex * vct_scr, cuDoubleComplex * vct_des, int n_ham)
{
  int ID;
  ID = blockDim.x * blockIdx.x + threadIdx.x;
  if(ID < n_ham)
  {
    vct_des[ID].x = vct_scr[ID].x;
    vct_des[ID].y = vct_scr[ID].y;
  }
}


__global__ void mv_kernel_1_GPU(int n_ham, int col, cuDoubleComplex * xVal_kr, cuDoubleComplex * y_kr, cuDoubleComplex * Finalans_kr){
			//kernel func varaibles
			int ID, row_start, row_end, jj;
			double dot_x, dot_img;
			ID = blockDim.x*blockIdx.x+threadIdx.x;
				if(ID < n_ham){
					dot_x=0.0;
					dot_img=0.0; 
					//row_start = csrRowPtr_kr[row]-1;
					//row_end = csrRowPtr_kr[row+1]-1;
					//(x + yi)(u + vi) = (xu ? yv) + (xv + yu)i. 
					for(jj = 0; jj < col; jj++ ){
						  dot_x += ((xVal_kr[(jj*n_ham)+ID].x * y_kr[jj].x) - (xVal_kr[(jj*n_ham)+ID].y * y_kr[jj].y)); 
						  dot_img  += ((xVal_kr[(jj*n_ham)+ID].x * y_kr[jj].y) + (xVal_kr[(jj*n_ham)+ID].y * y_kr[jj].x)); 
						}
		
					Finalans_kr[ID].x = dot_x;  
					Finalans_kr[ID].y = dot_img;  
					
				}
		}





void setdeviceinit_() {
char * localRankStr = NULL;
int rank = 0, devCount = 0;
cudaError_t cudaStat1;
// We extract the local rank initialization using an environment variable
if ((localRankStr = getenv(ENV_LOCAL_RANK)) != NULL)
{
rank = atoi(localRankStr);	
}

if ((localRankStr = getenv("CUDA_SELECT_GPU")) != NULL)
{
  rank = atoi(localRankStr);
}

cudaDeviceReset();
cudaThreadExit();
cudaGetDeviceCount(&devCount);
printf("device count %d %d %d\n", devCount, rank, rank%devCount);

cudaStat1 = cudaSetDevice(rank);
if(cudaStat1 != cudaSuccess)
printf("ERROR DEVICE SET FAILED\n");
}



__global__ void spmv_csr_hybrid_kernel_1_GPU(int dimRow, const int* rowPtrs, const int* colIdxs, const float* values, 
                                       const cuDoubleComplex* x, cuDoubleComplex* y, int repeat, int coop, int offset)
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
         int rowPtr = rowPtrs[offset+i]-1;
         int stop = rowPtrs[offset+i+1]-1-rowPtr;
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
         if (coopIdx == 0) y[offset+i] = sdata[tid];
         i += blockDim.x/coop;
      }
    }

}



__global__ void spmv_csr_hybrid_kernel_1_GPU(int num_rows, const int* rowPtrs, const int* colIdxs, const float* values, 
                                       const cuDoubleComplex* x, cuDoubleComplex* y, const int offset)
{
   int i = (blockIdx.x * blockDim.x + threadIdx.x);
   if (i<num_rows)
   {
     cuDoubleComplex rowSum;
     rowSum.x = 0.0; rowSum.y = 0.0;
     int row_start = rowPtrs[offset+i]-1;
     int row_end = rowPtrs[offset+i+1]-1;
     for (int j=row_start; j<row_end; j++)
     {
        rowSum.x += values[j] * x[colIdxs[j]-1].x;
        rowSum.y += values[j] * x[colIdxs[j]-1].y;
     }
     y[offset+i] = rowSum;
   }
}


__global__ void vct_pls_scl_mul_vct_kernel_1_GPU(cuDoubleComplex *vct1, cuDoubleComplex * vct2, cuDoubleComplex scalar, int n_ham)
{
  int ID;
  ID = blockDim.x * blockIdx.x + threadIdx.x;
  if(ID < n_ham)
  {
    vct1[ID].x =  vct1[ID].x + scalar.x * vct2[ID].x - scalar.y * vct2[ID].y;
    vct1[ID].y =  vct1[ID].y + scalar.x * vct2[ID].y + scalar.y * vct2[ID].x;
  }
}
