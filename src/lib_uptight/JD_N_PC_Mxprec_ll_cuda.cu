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
#include <complex.h>



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



void gmres(float *valptr_device_real, int *rowptr_device_real, int *colptr_device_real, float *valptr_device_img, int *rowptr_device_img, int *colptr_device_img, int n_ham, int size_mat_real, cuDoubleComplex * r_device, double ls_tol, int restart, int maxit, cuDoubleComplex *Q_bar_device, cuDoubleComplex *t_device, int numBlock, int threadPerBlock, int coop, int repeat, int numBlocksMul, int k, cusparseHandle_t cusp_handle_ls, cusparseMatDescr_t cusp_descra_ls, cublasHandle_t cubl_handle_ls, int *ls_counter, int * shift_init_Mi, int * shift_end_Mi, int overlap_high, int overlap_low, int num_procs, int id,  MPI_Comm upt_comm);

static cudaError_t set_gpu_id(int gpu);

__global__ void vct_div_slr_kernel(cuDoubleComplex * scr, cuDoubleComplex * des, int n_ham, double Slr);
__global__ void vct1_sub_mul_vct_kernel(cuDoubleComplex * vct1, cuDoubleComplex * vct2, int n_ham, cuDoubleComplex * dot);
__global__ void vct_sub_scl_mul_vct_kernel(cuDoubleComplex *vct1, cuDoubleComplex * vct2, double scalar, int n_ham);
__global__ void shift_A_kernel(float * val, int * row, int * col, int n_ham, float shift, int offset);
__global__ void cpy_vct_1_to_vct_2_kernel(cuDoubleComplex * vct_scr, cuDoubleComplex * vct_des, int n_ham);
__global__ void vct_1_div_asg_to_vct_2_kernel(cuDoubleComplex * vct_1, cuDoubleComplex * vct_2, int n_ham, double *beta);
__global__ void vct_div_slr_minus_scl_vct_kernel(cuDoubleComplex * scr1, cuDoubleComplex * scr2, cuDoubleComplex * des, int n_ham, double Slr1, double Slr2);
__global__ void vct1_add_vct2_asg_vct3_kernel(cuDoubleComplex *vct1, cuDoubleComplex * vct2, cuDoubleComplex * vct3, int n_ham);
__global__ void vct1_neg_asg_vct2_kernel(cuDoubleComplex *vct1, cuDoubleComplex * vct2, int n_ham);
__global__ void vct1_sub_vct2_asg_vct3_kernel(cuDoubleComplex *vct1, cuDoubleComplex * vct2, cuDoubleComplex * vct3, int n_ham);
__global__ void mv_kernel(int n_ham, int col, cuDoubleComplex * xVal_kr, cuDoubleComplex * y_kr, cuDoubleComplex * Finalans_kr);
__global__ void vct_pls_scl_mul_vct_kernel_jd(cuDoubleComplex *vct1, cuDoubleComplex * vct2, cuDoubleComplex scalar, int n_ham);
__global__ void spmv_csr_hybrid_kernel(int num_rows, const int* rowPtrs, const int* colIdxs, const float* values, const cuDoubleComplex* x, cuDoubleComplex* y, const int offset);
__global__ void spmv_csr_hybrid_kernel(int dimRow, const int* rowPtrs, const int* colIdxs, const float* values, const cuDoubleComplex* x, cuDoubleComplex* y, int repeat, int coop, int offset);

extern "C" {
void jd_single_gpu_no_pc_split_mxprec_pal_(int * N_ham, int * Size_Mat_real, int * Size_Mat_img, float * valptr_real, int * rowptr_real, int * colptr_real, float * valptr_img, int * rowptr_img, int * colptr_img, char * sparse_fmt, int * Band_type, double * JD_tol, double * Shift,  int * JD_Min_step, int * JD_Max_step, int * Num_ev, double * lambda_out, cuDoubleComplex * eigen_vec_out, double * LS_tol, int * LS_restart, int *LS_maxit, int *col_ind_low, int *col_ind_high, int *shift_init, int *shift_end, int * NUM_PROCS, int * ID,  int * UPT_COMM);

void setdevicebeforeinit_();

} 



cudaError_t set_gpu_id(int gpu)
{

  //if (gpu < 0)
  {
    // read it from the environment
    char* gpu_env = getenv("CUDA_SELECT_GPU");
    if (gpu_env != NULL)
      gpu = atoi(gpu_env);
printf("gpu_env = %s, gpu = %d\n", gpu_env, gpu);
  }

  return cudaSetDevice(gpu);
}




void jd_single_gpu_no_pc_split_mxprec_pal_(int * N_ham, int * Size_Mat_real, int * Size_Mat_img, float * valptr_real, int * rowptr_real, int * colptr_real, float * valptr_img, int * rowptr_img, int * colptr_img, char * sparse_fmt, int * Band_type, double * JD_tol, double * Shift,  int *JD_Min_step, int * JD_Max_step, int * Num_ev, double * lambda_out, cuDoubleComplex * eigen_vec_out, double * LS_tol, int * LS_restart, int *LS_maxit, int *col_ind_low, int *col_ind_high, int *shift_init, int *shift_end, int * NUM_PROCS, int * ID,  int * UPT_COMM)
{



int num_procs = * NUM_PROCS;
int id = * ID;
MPI_Comm upt_comm;

upt_comm = MPI_Comm_f2c( * UPT_COMM);


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

//printf("JD_tol %f\n", *JD_tol);

setdevicebeforeinit_();

cuDoubleComplex * V0;
cuDoubleComplex * X_bar;
float * valptr_device_real;
int * rowptr_device_real, * colptr_device_real;
float * valptr_device_img;
int * rowptr_device_img, * colptr_device_img;
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
cuDoubleComplex * temp1_device, * temp2_device;
cuDoubleComplex * w_bar_device;
cuDoubleComplex * r_device;
cuDoubleComplex * u;
double * s_v;
double * s_v_device;
cuDoubleComplex * ubar_device;
double * norm_V0;
const cuDoubleComplex jcmpx=make_cuDoubleComplex(0.0,1.0);
cuDoubleComplex * mxv_temp_device;
cuDoubleComplex * dot_temp, *y, * dot_temp_1, * dot_temp_2, * dot_temp_device;
double * norm_temp;
cuDoubleComplex * temp_host_v, * temp_host_w;


cuDoubleComplex * v, * w;

cudaMallocHost((void**)&v, n_ham*jd_max_step*sizeof(cuDoubleComplex));
cudaMallocHost((void**)&w, n_ham*jd_max_step*sizeof(cuDoubleComplex));
cudaMallocHost((void**)&temp_host_v,n_ham*jd_max_step*sizeof(cuDoubleComplex));
cudaMallocHost((void**)&temp_host_w,n_ham*jd_max_step*sizeof(cuDoubleComplex));

cudaStream_t stream1, stream2, stream3, stream4;
cudaStreamCreate(&stream1);
cudaStreamCreate(&stream2);
cudaStreamCreate(&stream3);
cudaStreamCreate(&stream4);


int *ls_counter;
ls_counter = (int *) malloc(sizeof(int));
*ls_counter = 0;

lapack_int IL=0, IU=0, *M_out, *isuppz;
double  ABSTOL = 0.001;
double VL=0, VU=0;


MPI_Request * reqs;
MPI_Status status;

reqs=(MPI_Request *) malloc(sizeof(MPI_Request)*(num_procs+1));


int * shift_init_Mi;
int * shift_end_Mi;

shift_init_Mi=(int *) malloc(num_procs*sizeof(int));
shift_end_Mi=(int *) malloc(num_procs*sizeof(int));

shift_init_Mi[id] = *shift_init;
shift_end_Mi[id] = *shift_end;


MPI_Barrier(upt_comm);
MPI_Gather(&shift_init_Mi[id], 1, MPI_INT, &shift_init_Mi[id], 1, MPI_INT, 0, upt_comm);
MPI_Gather(&shift_end_Mi[id], 1, MPI_INT, &shift_end_Mi[id], 1, MPI_INT, 0, upt_comm);
MPI_Barrier(upt_comm);


MPI_Bcast(shift_init_Mi,num_procs, MPI_INT, 0, upt_comm);
MPI_Bcast(shift_end_Mi,num_procs, MPI_INT, 0, upt_comm);


int overlap_high = 0;
int overlap_low = 0;


if(id == 0){
    
    for (int f_loop = 0 ; f_loop <= num_procs-2; f_loop++){
    	if(col_ind_low[f_loop+1] < shift_init_Mi[f_loop]){
    	printf("use less number of nodes in MPI\n");
    	exit(0);
    	}
    }
}

MPI_Barrier(upt_comm);



 for(int f_loop = 0; f_loop <=num_procs-1; f_loop++){
 	if (overlap_high < col_ind_high[f_loop]-shift_end_Mi[f_loop]) 
                                      overlap_high = col_ind_high[f_loop]-shift_end_Mi[f_loop];
 }


 for(int f_loop = 0; f_loop <=num_procs-1; f_loop++){
 	if (overlap_low < shift_init_Mi[f_loop]-col_ind_low[f_loop]) 
                                      overlap_low = shift_init_Mi[f_loop]-col_ind_low[f_loop]; 
}

//if(id == 0)
//printf(" overlap_high, overlap_low  id %d, %d, %d\n", overlap_high, overlap_low, id);

int overlap_tol = overlap_high + overlap_low;


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
norm_temp = (double *) malloc(sizeof(double));
y = (cuDoubleComplex *) malloc(sizeof(cuDoubleComplex));
dot_temp = (cuDoubleComplex *) malloc(sizeof(cuDoubleComplex));
dot_temp_1 = (cuDoubleComplex *) malloc(sizeof(cuDoubleComplex));
dot_temp_2 = (cuDoubleComplex *) malloc(sizeof(cuDoubleComplex));

cudaError_t cudaStat1,cudaStat2,cudaStat3,cudaStat4,cudaStat5,cudaStat6,cudaStat7,cudaStat8,cudaStat9,cudaStat10,cudaStat11,cudaStat12,cudaStat13,cudaStat14,cudaStat15,cudaStat16, cudaStat17,cudaStat18,cudaStat19,cudaStat20,cudaStat21,cudaStat22,cudaStat23,cudaStat24,cudaStat25,cudaStat26,cudaStat27;

lapack_int error_code=0;

cusparseStatus_t cusp_status;
cusparseHandle_t cusp_handle=0;
cusparseMatDescr_t cusp_descra=0;
cublasStatus_t cubl_status;
cublasHandle_t cubl_handle=0;

threadPerBlock = BLOCK_SIZE;                      // blockSize
numBlock=(n_ham/threadPerBlock)+1;                // gridSize

int coop= 16  ;    
int repeat= 2 ;   // repeat = nrows * coop/ BLOCK_SIZE_MUL 
int numBlocksMul =(n_ham*coop-1)/(repeat*BLOCK_SIZE_MUL)+1;  //gridSize
 
  
if(id == 0) printf(" (CUDA) JD, threadsPerBlock %d \n", threadPerBlock);
if(id == 0) printf(" (CUDA) gridSize  %d \n",numBlock);
if(id == 0) printf(" (CUDA) JD, threadsPerBlock %d \n", BLOCK_SIZE_MUL);
if(id == 0) printf(" (CUDA) gridSize for Mxv %d \n",numBlocksMul);



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

if(id == 0) {

for(int counter=0; counter < n_ham; counter++)
{
	V0[counter].x = counter; //rand();
	V0[counter].y = counter*3; //rand();
}

}

MPI_Barrier(upt_comm);

MPI_Bcast(V0,n_ham, MPI_DOUBLE_COMPLEX, 0, upt_comm);


X_bar = (cuDoubleComplex *) malloc(sizeof(cuDoubleComplex)*n_ham*num_ev);



double memory_device = (sizeof(float)*(size_mat_real+size_mat_img)+sizeof(int)*(size_mat_real+n_ham+size_mat_img+n_ham)+sizeof(double)*(jd_max_step)+sizeof(cuDoubleComplex)*(n_ham*12+n_ham*num_ev+n_ham*jd_max_step*3)+(n_ham*((ls_restart+1)+ls_restart))+(jd_max_step*jd_max_step*2)+((ls_restart+1)*(ls_restart+1))+((ls_restart+1)*ls_restart))/1000000;

if(id == 0)
printf("memory needed on GPU is ~ %f MB\n", memory_device); 


cudaStat1 = cudaMalloc((void**)&colptr_device_real,size_mat_real*sizeof(int));
cudaStat2 = cudaMalloc((void**)&rowptr_device_real,(n_ham+1)*sizeof(int));
cudaStat3 = cudaMalloc((void**)&valptr_device_real,size_mat_real*sizeof(float));
cudaStat6 = cudaMalloc((void**)&t_device,n_ham*sizeof(cuDoubleComplex));
            cudaMemset(t_device, 0, n_ham*sizeof(cuDoubleComplex));  
cudaStat9 = cudaMalloc((void**)&y_device,sizeof(cuDoubleComplex));
cudaStat11 = cudaMalloc((void**)&colptr_device_img,size_mat_img*sizeof(int));
cudaStat14 = cudaMalloc((void**)&rowptr_device_img,(n_ham+1)*sizeof(int));
cudaStat16 = cudaMalloc((void**)&valptr_device_img,size_mat_img*sizeof(float));
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
             cudaMemset(r_device, 0, n_ham*sizeof(cuDoubleComplex));
             cudaMalloc((void**)&dot_temp_device,sizeof(cuDoubleComplex));
if( 
(cudaStat1 != cudaSuccess) ||
(cudaStat2 != cudaSuccess) ||
(cudaStat3 != cudaSuccess) ||
(cudaStat5 != cudaSuccess) ||
(cudaStat6 != cudaSuccess) ||
(cudaStat9 != cudaSuccess) ||
(cudaStat11 != cudaSuccess) ||
(cudaStat14 != cudaSuccess) ||
(cudaStat12 != cudaSuccess) ||
(cudaStat13 != cudaSuccess) ||
(cudaStat15 != cudaSuccess) ||
(cudaStat16 != cudaSuccess) ||
(cudaStat21 != cudaSuccess) ||
(cudaStat23 != cudaSuccess) ||
(cudaStat24 != cudaSuccess) ||
(cudaStat25 != cudaSuccess) ||
(cudaStat26 != cudaSuccess) ||
(cudaStat27 != cudaSuccess))
{
  CLEANUP("Device malloc failed\n");
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
  CLEANUP("copy to device failed loc\n");
}


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

//if(id == 0)
//printf("cublas and cusparse initilized\n");


// V0 = V0/norm(V0); t_device = V0
cubl_status = cublasDznrm2(cubl_handle, n_ham, t_device, 1, norm_V0_device);
cudaStat1 = cudaMemcpy(norm_V0, norm_V0_device,(size_t)(sizeof(double)), cudaMemcpyDeviceToHost);

//printf("Norm of V0 done %f\n", *norm_V0);
if(id == 0)
{
  vct_div_slr_kernel<<<numBlock, threadPerBlock, 0, stream1>>>( &t_device[shift_init_Mi[id]-1], &t_device[shift_init_Mi[id]-1], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high, *norm_V0);
}
else if(id != 0 && id != num_procs-1)
{ 
  vct_div_slr_kernel<<<numBlock, threadPerBlock, 0, stream1>>>( &t_device[shift_init_Mi[id]-1-overlap_low], &t_device[shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol, *norm_V0);
}
else
{
  vct_div_slr_kernel<<<numBlock, threadPerBlock, 0, stream1>>>( &t_device[shift_init_Mi[id]-1-overlap_low], &t_device[shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low, *norm_V0);
}


// t=v0; k = 0; m = 0; X_bar = [];
int k = 0;
int m = 0;

//cpy_vct_1_to_vct_2_kernel<<<numBlock, threadPerBlock>>>( V0_device, t_device, n_ham);

*s_v = 0;
int count = 0;
double teta = 0;

     // A = (A - shift * I)

    tol_shift = teta-shift;
    shift_A_kernel<<<numBlock, threadPerBlock, 0, stream2>>>(valptr_device_real, rowptr_device_real, colptr_device_real, (shift_end_Mi[id]-shift_init_Mi[id])+1, tol_shift, shift_init_Mi[id]-1);



while(k < num_ev)
{

    cudaStat1 = cudaMalloc((void**)&w_device,n_ham*sizeof(cuDoubleComplex));
                cudaMemset(w_device, 0, n_ham*sizeof(cuDoubleComplex));

    cudaStat2 = cudaMalloc((void**)&v_device,n_ham*sizeof(cuDoubleComplex));
                cudaMemset(v_device, 0, n_ham*sizeof(cuDoubleComplex));

    if(cudaStat1 != cudaSuccess || cudaStat1 != cudaSuccess)
      {
         CLEANUP("Device malloc failed loc 5\n");
      }

      cudaStat1 = cudaMalloc((void**)&vct_device,n_ham*sizeof(cuDoubleComplex));
                  cudaMemset(vct_device, 0, sizeof(cuDoubleComplex)*n_ham);
      cudaStat2 = cudaMalloc((void**)&mxv_temp_device,n_ham*sizeof(cuDoubleComplex));
                  cudaMemset(mxv_temp_device, 0, sizeof(cuDoubleComplex)*n_ham);
      if(cudaStat1 != cudaSuccess || cudaStat2 != cudaSuccess)
      {
         CLEANUP("Device malloc failed loc 5\n");
      }

     cudaStat1=cudaDeviceSynchronize();
     MPI_Barrier(upt_comm);

      // vct = A * t;
      //cusp_status= cusparseZcsrmv(cusp_handle,CUSPARSE_OPERATION_NON_TRANSPOSE, n_ham, n_ham, size_mat_real, &scalar_1, cusp_descra, valptr_device, rowptr_device, colptr_device, t_device, &scalar_2, vct_device);

     spmv_csr_hybrid_kernel<<<numBlocksMul, BLOCK_SIZE_MUL, 0, stream3>>>((shift_end_Mi[id]-shift_init_Mi[id])+1, rowptr_device_real, colptr_device_real, valptr_device_real, t_device, vct_device, repeat, coop, shift_init_Mi[id]-1);

     spmv_csr_hybrid_kernel<<<numBlock, BLOCK_SIZE, 0, stream4>>>((shift_end_Mi[id]-shift_init_Mi[id])+1, rowptr_device_img, colptr_device_img, valptr_device_img, t_device, mxv_temp_device, shift_init_Mi[id]-1);

     cudaStat1=cudaDeviceSynchronize();

     vct_pls_scl_mul_vct_kernel_jd<<<numBlock, threadPerBlock>>>(&vct_device[shift_init_Mi[id]-1], &mxv_temp_device[shift_init_Mi[id]-1], jcmpx, (shift_end_Mi[id]-shift_init_Mi[id])+1);

     cudaStat1=cudaDeviceSynchronize();

     cudaFree(mxv_temp_device);


     if (num_procs > 1){
     for(int rank =0; rank <=num_procs-2; rank++){
     if(id == rank+1)    MPI_Irecv(&vct_device[shift_end_Mi[rank]-overlap_low], overlap_low,MPI_DOUBLE_COMPLEX,rank, 111+rank,upt_comm,&reqs[rank+1]);
     if(id == rank)      MPI_Isend(&vct_device[shift_end_Mi[rank]-overlap_low], overlap_low,MPI_DOUBLE_COMPLEX,rank+1, 111+rank,upt_comm, &reqs[rank+1]);

     if (id == rank+1) MPI_Wait(&reqs[rank+1], &status);
     if (id == rank)   MPI_Wait(&reqs[rank+1], &status);
     }
     }

     //MPI_Barrier(upt_comm);


     if (num_procs > 1){
     for(int rank =0; rank <=num_procs-2; rank++){
     if(id == rank)    MPI_Irecv(&vct_device[shift_init_Mi[rank+1]-1],overlap_high, MPI_DOUBLE_COMPLEX,rank+1,8+rank,upt_comm, &reqs[rank+1]);
     if(id == rank+1)  MPI_Isend(&vct_device[shift_init_Mi[rank+1]-1],overlap_high, MPI_DOUBLE_COMPLEX,rank,8+rank,upt_comm, &reqs[rank+1]);

     if (id == rank+1) MPI_Wait(&reqs[rank+1], &status);
     if (id == rank)   MPI_Wait(&reqs[rank+1], &status);
     }
     }

     MPI_Barrier(upt_comm);

     //cudaStat1=cudaDeviceSynchronize();

    for(int i = 0; i < m; i++)
    {

     
     //if((i+1) != m){
     //cudaStat3 = cudaMemcpy(v_device, &v[i*n_ham], (size_t)(sizeof(cuDoubleComplex)*n_ham), cudaMemcpyHostToDevice);
     if(id ==0)
     cudaStat1 = cudaMemcpy(&v_device[shift_init_Mi[id]-1], &v[(i*n_ham)+shift_init_Mi[id]-1],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high)*sizeof(cuDoubleComplex)), cudaMemcpyHostToDevice);
     else if(id != 0 && id != num_procs-1)
     cudaStat1 = cudaMemcpy(&v_device[shift_init_Mi[id]-1-overlap_low], &v[(i*n_ham)+shift_init_Mi[id]-1-overlap_low],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol)*sizeof(cuDoubleComplex)), cudaMemcpyHostToDevice);
     else
     cudaStat1 = cudaMemcpy(&v_device[shift_init_Mi[id]-1-overlap_low], &v[(i*n_ham)+shift_init_Mi[id]-1-overlap_low],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low)*sizeof(cuDoubleComplex)), cudaMemcpyHostToDevice);

     //cudaStat4 = cudaMemcpy(w_device, &w[i*n_ham], (size_t)(sizeof(cuDoubleComplex)*n_ham), cudaMemcpyHostToDevice);
     if(id ==0)
     cudaStat1 = cudaMemcpy(&w_device[shift_init_Mi[id]-1], &w[(i*n_ham)+shift_init_Mi[id]-1],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high)*sizeof(cuDoubleComplex)), cudaMemcpyHostToDevice);
     else if(id != 0 && id != num_procs-1)
     cudaStat1 = cudaMemcpy(&w_device[shift_init_Mi[id]-1-overlap_low], &w[(i*n_ham)+shift_init_Mi[id]-1-overlap_low],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol)*sizeof(cuDoubleComplex)), cudaMemcpyHostToDevice);
     else
     cudaStat1 = cudaMemcpy(&w_device[shift_init_Mi[id]-1-overlap_low], &w[(i*n_ham)+shift_init_Mi[id]-1-overlap_low],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low)*sizeof(cuDoubleComplex)), cudaMemcpyHostToDevice);

     //}

     // y = w(:,i)' * vct;
     cubl_status = cublasZdotc(cubl_handle, (shift_end_Mi[id]-shift_init_Mi[id])+1, &w_device[shift_init_Mi[id]-1], 1, &vct_device[shift_init_Mi[id]-1], 1, y_device);

     //MPI_Barrier(upt_comm); 

     cudaStat1 = cudaMemcpy(dot_temp, y_device,(size_t)(sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);

     //MPI_Barrier(upt_comm);

     MPI_Allreduce(dot_temp, y, 1 , MPI_DOUBLE_COMPLEX, MPI_SUM, upt_comm);

     cudaStat1 = cudaMemcpy(y_device, y,(size_t)(sizeof(cuDoubleComplex)), cudaMemcpyHostToDevice);

     //cudaDeviceSynchronize();


     // vct = vct - y * w(:,i);

     if(id == 0)
     {
     vct1_sub_mul_vct_kernel<<<numBlock, threadPerBlock, 0, stream1>>>(&vct_device[shift_init_Mi[id]-1], &w_device[shift_init_Mi[id]-1], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high, y_device);
     }
     else if(id != 0 && id != num_procs-1)
     {
     vct1_sub_mul_vct_kernel<<<numBlock, threadPerBlock, 0, stream1>>>(&vct_device[shift_init_Mi[id]-1-overlap_low], &w_device[shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol, y_device);
     }
     else
     {
     vct1_sub_mul_vct_kernel<<<numBlock, threadPerBlock, 0, stream1>>>(&vct_device[shift_init_Mi[id]-1-overlap_low], &w_device[shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low, y_device);
     }


     // t = t - y * v(:,i);
     if(id == 0)
     {
     vct1_sub_mul_vct_kernel<<<numBlock, threadPerBlock, 0, stream2>>>(&t_device[shift_init_Mi[id]-1], &v_device[shift_init_Mi[id]-1], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high, y_device);
     }
     else if(id != 0 && id != num_procs-1)
     {
     vct1_sub_mul_vct_kernel<<<numBlock, threadPerBlock, 0, stream2>>>(&t_device[shift_init_Mi[id]-1-overlap_low], &v_device[shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol, y_device);
     }
     else
     {
     vct1_sub_mul_vct_kernel<<<numBlock, threadPerBlock, 0, stream2>>>(&t_device[shift_init_Mi[id]-1-overlap_low], &v_device[shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low, y_device);
     }

    cudaDeviceSynchronize();
     }

     m = m+1;

     //cubl_status = cublasDznrm2(cubl_handle, n_ham, vct_device, 1, norm_vct_device);
     cubl_status = cublasZdotc(cubl_handle, (shift_end_Mi[id]-shift_init_Mi[id])+1, &vct_device[shift_init_Mi[id]-1], 1, &vct_device[shift_init_Mi[id]-1], 1, dot_temp_device);

     //MPI_Barrier(upt_comm);

     cudaStat1 = cudaMemcpy(dot_temp_1, dot_temp_device,(size_t)(sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);

     //MPI_Barrier(upt_comm);

     MPI_Allreduce(dot_temp_1, dot_temp_2, 1 , MPI_DOUBLE_COMPLEX, MPI_SUM, upt_comm);

     //MPI_Barrier(upt_comm);

     *norm_vct = sqrt(dot_temp_2[0].x);

     cudaStat1 = cudaMemcpy(norm_vct_device, norm_vct,(size_t)(sizeof(double)), cudaMemcpyHostToDevice);



    // w(:,m) = vct/norm(vct);vct_div_slr_kernel(cuDoubleComplex * scr, cuDoubleComplex * des, int n_ham, double * Slr)
    // v(:,m) = t/norm(vct);
 
    if(id == 0)
    {
     vct_div_slr_kernel<<<numBlock, threadPerBlock, 0, stream1>>>(&vct_device[shift_init_Mi[id]-1], &w_device[shift_init_Mi[id]-1], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high, *norm_vct);
    }
    else if(id != 0 && id != num_procs-1)
    {
    vct_div_slr_kernel<<<numBlock, threadPerBlock, 0, stream1>>>(&vct_device[shift_init_Mi[id]-1-overlap_low], &w_device[shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol, *norm_vct);
    }
    else
    {
    vct_div_slr_kernel<<<numBlock, threadPerBlock, 0, stream1>>>(&vct_device[shift_init_Mi[id]-1-overlap_low], &w_device[shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low, *norm_vct);
    }

    cudaStat1 = cudaMemcpyAsync(&w[(m-1)*n_ham], w_device,(size_t)(sizeof(cuDoubleComplex)*n_ham), cudaMemcpyDeviceToHost, stream1);

    //cudaStat1 = cudaMemcpy(v_device, &v[(m-1)*n_ham],(size_t)(sizeof(cuDoubleComplex)*n_ham), cudaMemcpyHostToDevice);
    if(id == 0)
    {
     vct_div_slr_kernel<<<numBlock, threadPerBlock, 0, stream2>>>(&t_device[shift_init_Mi[id]-1], &v_device[shift_init_Mi[id]-1], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high, *norm_vct);
    }
    else if(id != 0 && id != num_procs-1)
    {
    vct_div_slr_kernel<<<numBlock, threadPerBlock, 0, stream2>>>(&t_device[shift_init_Mi[id]-1-overlap_low], &v_device[shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol, *norm_vct);
    }
    else
    {
    vct_div_slr_kernel<<<numBlock, threadPerBlock, 0, stream2>>>(&t_device[shift_init_Mi[id]-1-overlap_low], &v_device[shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low, *norm_vct);
    }

    cudaStat1 = cudaMemcpyAsync(&v[(m-1)*n_ham], v_device,(size_t)(sizeof(cuDoubleComplex)*n_ham), cudaMemcpyDeviceToHost, stream2);


    // alloc M(m,m);
    cudaStat1 = cudaMalloc((void**)&M_device,m*m*sizeof(cuDoubleComplex));
                cudaMemset(M_device, 0, m*m*sizeof(cuDoubleComplex));

    cudaStat2 = cudaMallocHost((void**)&M, sizeof(cuDoubleComplex)*m*m);

    cudaDeviceSynchronize();
   
    cudaFree(vct_device);

    for(int j=0; j < m; j++)
    { 


    for(int i = 0; i < m; i++)
    {

     // M(i,m) = w(:,i)' * v(:,m);
     cubl_status = cublasZdotc(cubl_handle, (shift_end_Mi[id]-shift_init_Mi[id])+1, &w[(i*n_ham)+shift_init_Mi[id]-1], 1, &v[(j*n_ham)+shift_init_Mi[id]-1], 1, &M_device[j*m+i]);

     //MPI_Barrier(upt_comm); 

     cudaStat1 = cudaMemcpy(dot_temp, &M_device[j*m+i],(size_t)(sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);

     //MPI_Barrier(upt_comm);

     MPI_Allreduce(dot_temp, &M[j*m+i], 1 , MPI_DOUBLE_COMPLEX, MPI_SUM, upt_comm);


     // M(m,i) = w(:,m)' * v(:,i);
     M[m*i+j] = cuConj(M[j*m+i]);

    }

    // M(m,m) = w(:,m)' * v(:,m);

    }

    cudaStat1 = cudaMemcpy(M_device, M,(size_t)(sizeof(cuDoubleComplex)*m*m), cudaMemcpyHostToDevice);
    if(cudaStat1 != cudaSuccess)
    {
     CLEANUP("cudaMemcpy failed loc 3\n");
    }
    
    eigen_vec = (cuDoubleComplex *) malloc(sizeof(cuDoubleComplex)*m*m);
    eigen_val = (double *) malloc(sizeof(double)*m);
    cudaStat1 = cudaMalloc((void**)&eigen_vec_device,sizeof(cuDoubleComplex)*m*m);
    cudaStat2 = cudaMalloc((void**)&eigen_val_device,sizeof(double)*m);
    if(cudaStat1 != cudaSuccess || cudaStat2 != cudaSuccess){
      CLEANUP("cudaMemcpy failed loc 4\n");}

       error_code = LAPACKE_zheevr( LAPACK_COL_MAJOR, 'V', 'A', 'L', m, (MKL_Complex16 *)M, m, VL, VU, IL, IU, ABSTOL, M_out, eigen_val, (MKL_Complex16 *)eigen_vec, m, isuppz );



     if(m != 1)
     {
     int index[m], temp_index;
     double temp;
     temp_ev = (cuDoubleComplex *) malloc(sizeof(cuDoubleComplex)*m);

     for(int i =0; i < m; i++)
     {
       index[i] = i;
     }

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

      // sort eigen_vector
      for(int i =0; i < m; i++)
      {
         cudaStat1 = cudaMemcpy(&eigen_vec_device[i*m], &eigen_vec[index[i]*m],(size_t)(sizeof(cuDoubleComplex)*m), cudaMemcpyHostToDevice);
         if(cudaStat1 != cudaSuccess)
         {
         CLEANUP("cudaMemcpy failed loc 5\n");
         }

      }
      }


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


      cudaStat1 = cudaMalloc((void**)&ubar_device,n_ham*sizeof(cuDoubleComplex));
                  cudaMemset(ubar_device, 0, n_ham*sizeof(cuDoubleComplex));
      if(cudaStat1 != cudaSuccess)
      {
      CLEANUP("Device malloc failed\n");
      }
     
      cudaStat1 = cudaMalloc((void**)&w_bar_device,n_ham*sizeof(cuDoubleComplex));
                  cudaMemset(w_bar_device, 0, n_ham*sizeof(cuDoubleComplex));
      if(cudaStat1 != cudaSuccess)
      {
         CLEANUP("Device malloc failed loc 11\n");
      }

      cublasSetStream(cubl_handle, stream1);
      // u_bar = v * eig_mat(:,1);
      cubl_status = cublasZgemv(cubl_handle, CUBLAS_OP_N, n_ham, m, scalar1_device, v, n_ham, eigen_vec_device, 1, scalar2_device, ubar_device, 1);

      cublasSetStream(cubl_handle, stream2);
      // w_bar = w * eig_mat(:,1);
      cubl_status = cublasZgemv(cubl_handle, transaction, n_ham, m, scalar1_device, w, n_ham, eigen_vec_device, 1, scalar2_device, w_bar_device, 1); 


        cudaDeviceSynchronize();
        cublasSetStream(cubl_handle, NULL);
         

      // s_u = norm(u_bar);
      //cubl_status = cublasDznrm2(cubl_handle, n_ham, ubar_device, 1, s_u_device);
      cubl_status = cublasZdotc(cubl_handle, (shift_end_Mi[id]-shift_init_Mi[id])+1, &ubar_device[shift_init_Mi[id]-1], 1, &ubar_device[shift_init_Mi[id]-1], 1, dot_temp_device);

      //MPI_Barrier(upt_comm);

      cudaStat1 = cudaMemcpy(dot_temp_1, dot_temp_device,(size_t)(sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);

      //MPI_Barrier(upt_comm);

      MPI_Allreduce(dot_temp_1, dot_temp_2, 1 , MPI_DOUBLE_COMPLEX, MPI_SUM, upt_comm);

      //MPI_Barrier(upt_comm);

      *s_u = sqrt(dot_temp_2[0].x);

      cudaStat1 = cudaMemcpy(s_u_device, s_u,(size_t)(sizeof(double)), cudaMemcpyHostToDevice);
      if(cudaStat1 != cudaSuccess)
      {
         CLEANUP("cudaMemcpy failed loc 10\n");
      }

      cudaStat1 = cudaMalloc((void**)&u_device,n_ham*sizeof(cuDoubleComplex));
                  cudaMemset(u_device, 0, n_ham*sizeof(cuDoubleComplex));
      if(cudaStat1 != cudaSuccess)
      {
         CLEANUP("Device malloc failed loc 10\n");
      }

      // u = u_bar/s_u; 
      if(id == 0)
      {
       vct_div_slr_kernel<<<numBlock, threadPerBlock>>>(&ubar_device[shift_init_Mi[id]-1], &u_device[shift_init_Mi[id]-1], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high, *s_u);
      }
      else if(id != 0 && id != num_procs-1)
      {
      vct_div_slr_kernel<<<numBlock, threadPerBlock>>>(&ubar_device[shift_init_Mi[id]-1-overlap_low], &u_device[shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol, *s_u);
      }
      else
      {
      vct_div_slr_kernel<<<numBlock, threadPerBlock>>>(&ubar_device[shift_init_Mi[id]-1-overlap_low], &u_device[shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low, *s_u);
      }

      cudaFree(ubar_device);

      // s_v = eig_val(1)/ (s_u * s_u);
      *s_v = eigen_val[0]/( (*s_u) * (*s_u));

//if(id == 0)
//printf("s_v %f\n", *s_v);

      cudaStat1 = cudaMemcpy(s_v_device, s_v,(size_t)(sizeof(double)), cudaMemcpyHostToDevice);
      if(cudaStat1 != cudaSuccess)
      {
         CLEANUP("cudaMemcpy failed loc 11\n");
      }


      // r = (w_bar / s_u) - (s_v * u);

      if(id == 0)
      {
       vct_div_slr_minus_scl_vct_kernel<<<numBlock, threadPerBlock>>>(&w_bar_device[shift_init_Mi[id]-1], &u_device[shift_init_Mi[id]-1], &r_device[shift_init_Mi[id]-1], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high, *s_u, *s_v);
      }
      else if(id != 0 && id != num_procs-1)
      {
      vct_div_slr_minus_scl_vct_kernel<<<numBlock, threadPerBlock>>>(&w_bar_device[shift_init_Mi[id]-1-overlap_low], &u_device[shift_init_Mi[id]-1-overlap_low], &r_device[shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol, *s_u, *s_v);
      }
      else
      {
      vct_div_slr_minus_scl_vct_kernel<<<numBlock, threadPerBlock>>>(&w_bar_device[shift_init_Mi[id]-1-overlap_low], &u_device[shift_init_Mi[id]-1-overlap_low], &r_device[shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low, *s_u, *s_v);
      }

      cudaFree(w_bar_device);


      // norm(r_bar);
      //cubl_status = cublasDznrm2(cubl_handle, n_ham, r_device, 1, norm_r_device);
      cubl_status = cublasZdotc(cubl_handle, (shift_end_Mi[id]-shift_init_Mi[id])+1, &r_device[shift_init_Mi[id]-1], 1, &r_device[shift_init_Mi[id]-1], 1, dot_temp_device);

      //MPI_Barrier(upt_comm);

      cudaStat1 = cudaMemcpy(dot_temp_1, dot_temp_device,(size_t)(sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);

      //MPI_Barrier(upt_comm);

      MPI_Allreduce(dot_temp_1, dot_temp_2, 1 , MPI_DOUBLE_COMPLEX, MPI_SUM, upt_comm);

      //MPI_Barrier(upt_comm);

      *norm_r = sqrt(dot_temp_2[0].x);

      cudaStat1 = cudaMemcpy(norm_r_device, norm_r,(size_t)(sizeof(double)), cudaMemcpyHostToDevice);



      cudaStat1 = cudaMemcpy(norm_r, norm_r_device,(size_t)(sizeof(double)), cudaMemcpyDeviceToHost);
      if(cudaStat1 != cudaSuccess)
          {
          CLEANUP("cudaMemcpy failed loc 14\n");
          }

      //while(norm(r) < tol)
      while(*norm_r < *JD_tol)
      {
        //printf("inside harmonic part\n");        

        k = k + 1;

        // X_bar = [X_bar, u];
        //cudaStat1 = cudaMemcpy(u, u_device,(size_t)(sizeof(cuDoubleComplex)*n_ham), cudaMemcpyDeviceToHost);
        if(id ==0)
        cudaStat1 = cudaMemcpy(&u[shift_init_Mi[id]-1], &u_device[shift_init_Mi[id]-1],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high)*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);
        else if(id != 0 && id != num_procs-1)
        cudaStat1 = cudaMemcpy(&u[shift_init_Mi[id]-1-overlap_low], &u_device[shift_init_Mi[id]-1-overlap_low],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol)*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);
        else
        cudaStat1 = cudaMemcpy(&u[shift_init_Mi[id]-1-overlap_low], &u_device[shift_init_Mi[id]-1-overlap_low],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low)*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);

        // gather u from all IDs
        if (num_procs > 1) 
        {
          for (int rank =0; rank <=  num_procs-1; rank++) 
             {
             if(id == 0 && rank != 0)  MPI_Irecv(&u[shift_init_Mi[rank]-1],(shift_end_Mi[rank]-shift_init_Mi[rank])+1, MPI_DOUBLE_COMPLEX, rank, 111+rank, upt_comm, &reqs[rank+1]);
             if(id != 0 && id == rank) MPI_Isend(&u[shift_init_Mi[rank]-1], (shift_end_Mi[rank]-shift_init_Mi[rank])+1, MPI_DOUBLE_COMPLEX, 0, 111+rank, upt_comm, &reqs[rank+1]);

             if (id != 0 && id == rank) MPI_Wait(&reqs[rank+1], &status);
             }

        }

        // brodcast u
        MPI_Barrier(upt_comm);
        MPI_Bcast(u,n_ham,MPI_DOUBLE_COMPLEX, 0, upt_comm);

        memcpy ( &X_bar[(k-1)*n_ham], u, sizeof(cuDoubleComplex)*n_ham );

        // lambda(k) = s_v + shift
        lambda[k-1] = *s_v + shift;
        //if(id == 0) printf("lambda[%d] = %f \n", k, lambda[k-1]);
        //printf("JD count = %d\n", count);

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
         
            for(int i =0; i < k; i++)
            {
              memcpy(&eigen_vec_out[i*n_ham], &X_bar[index[i]*n_ham], sizeof(cuDoubleComplex)*n_ham);
            }

            memcpy ( lambda_out, lambda, k*sizeof(double));

          //cudaDeviceReset();
          for(int i =0; i < k; i++){
            if(id == 0) printf("lambda[%d] = %f \n", i, lambda[i]);
            }
          printf("Total JD count = %d\n", count);
          //printf("Total GMRES iterations = %d\n", * ls_counter); 
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
        cudaStat1 = cudaMalloc((void**)&temp1_device,n_ham*sizeof(cuDoubleComplex));
                    cudaMemset(temp1_device, 0, n_ham*sizeof(cuDoubleComplex));
        cudaStat2 = cudaMalloc((void**)&temp2_device,n_ham*sizeof(cuDoubleComplex));
                    cudaMemset(temp2_device, 0, n_ham*sizeof(cuDoubleComplex));
        if(cudaStat1 != cudaSuccess)
        {
           CLEANUP("Device malloc failed loc 14.1\n");
        }

        for(int i = 0; i < m; i++)
        {
          cublasSetStream(cubl_handle, stream1);
          //v_temp(:,i) = v * eig_mat(:,i+1);
          cubl_status = cublasZgemv(cubl_handle, CUBLAS_OP_N, n_ham, m+1, scalar1_device, v, n_ham, &eigen_vec_device[(i+1)*(m+1)], 1, scalar2_device, temp1_device, 1);
if(cubl_status != CUBLAS_STATUS_SUCCESS)
{
 CLEANUP("cublasZgemv failed loc 15\n");
}

          cublasSetStream(cubl_handle, stream2);
          // w_temp(:,i) = w * eig_mat(:,i+1);
          cubl_status = cublasZgemv(cubl_handle, CUBLAS_OP_N, n_ham, m+1, scalar1_device, w, n_ham, &eigen_vec_device[(i+1)*(m+1)], 1, scalar2_device, temp2_device, 1);
if(cubl_status != CUBLAS_STATUS_SUCCESS)
{
 CLEANUP("cublasZgemv failed loc 16\n");
}
                   //cudaStat1 = cudaMemcpy(&temp_host[i*n_ham], temp_device,(size_t)(sizeof(cuDoubleComplex)*n_ham), cudaMemcpyDeviceToHost);
          if(id ==0)
          cudaStat1 = cudaMemcpyAsync(&temp_host_v[(i*n_ham)+shift_init_Mi[id]-1], &temp1_device[shift_init_Mi[id]-1],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high)*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost, stream1);
          else if(id != 0 && id != num_procs-1)
          cudaStat1 = cudaMemcpyAsync(&temp_host_v[(i*n_ham)+shift_init_Mi[id]-1-overlap_low], &temp1_device[shift_init_Mi[id]-1-overlap_low],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol)*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost, stream1);
          else
          cudaStat1 = cudaMemcpyAsync(&temp_host_v[(i*n_ham)+shift_init_Mi[id]-1-overlap_low], &temp1_device[shift_init_Mi[id]-1-overlap_low], (size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low)*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost, stream1);

          //cudaStat1 = cudaMemcpy(&temp_host[i*n_ham], temp_device,(size_t)(sizeof(cuDoubleComplex)*n_ham), cudaMemcpyDeviceToHost);
          if(id ==0)
          cudaStat1 = cudaMemcpyAsync(&temp_host_w[(i*n_ham)+shift_init_Mi[id]-1], &temp2_device[shift_init_Mi[id]-1],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high)*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost, stream2);
          else if(id != 0 && id != num_procs-1)
          cudaStat1 = cudaMemcpyAsync(&temp_host_w[(i*n_ham)+shift_init_Mi[id]-1-overlap_low], &temp2_device[shift_init_Mi[id]-1-overlap_low],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol)*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost, stream2);
          else
          cudaStat1 = cudaMemcpyAsync(&temp_host_w[(i*n_ham)+shift_init_Mi[id]-1-overlap_low], &temp2_device[shift_init_Mi[id]-1-overlap_low], (size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low)*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost, stream2);
         }

         // v = v_temp(:,(1:m));
        cudaStat2 = cudaMemcpy(v, temp_host_v,(size_t)(sizeof(cuDoubleComplex)*n_ham*m), cudaMemcpyHostToHost);

        //cudaStat2 = cudaMemcpyAsync(v_device, v,(size_t)(sizeof(cuDoubleComplex)*n_ham), cudaMemcpyHostToDevice, stream1);
         if(id ==0)
         cudaStat1 = cudaMemcpyAsync(&v_device[shift_init_Mi[id]-1], &v[shift_init_Mi[id]-1],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high)*sizeof(cuDoubleComplex)), cudaMemcpyHostToDevice, stream1);
         else if(id != 0 && id != num_procs-1)
         cudaStat1 = cudaMemcpyAsync(&v_device[shift_init_Mi[id]-1-overlap_low], &v[shift_init_Mi[id]-1-overlap_low],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol)*sizeof(cuDoubleComplex)), cudaMemcpyHostToDevice, stream1);
         else
         cudaStat1 = cudaMemcpyAsync(&v_device[shift_init_Mi[id]-1-overlap_low], &v[shift_init_Mi[id]-1-overlap_low], (size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low)*sizeof(cuDoubleComplex)), cudaMemcpyHostToDevice, stream1);

         // w = w_temp(:,(1:m)); 
         cudaStat2 = cudaMemcpy(w, temp_host_w,(size_t)(sizeof(cuDoubleComplex)*n_ham*m), cudaMemcpyHostToHost);

         //cudaStat1 = cudaMemcpyAsync(w_device, w,(size_t)(sizeof(cuDoubleComplex)*n_ham), cudaMemcpyHostToDevice, stream2);
         if(id ==0)
         cudaStat1 = cudaMemcpyAsync(&w_device[shift_init_Mi[id]-1], &w[shift_init_Mi[id]-1],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high)*sizeof(cuDoubleComplex)), cudaMemcpyHostToDevice, stream2);
         else if(id != 0 && id != num_procs-1)
         cudaStat1 = cudaMemcpyAsync(&w_device[shift_init_Mi[id]-1-overlap_low], &w[shift_init_Mi[id]-1-overlap_low],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol)*sizeof(cuDoubleComplex)), cudaMemcpyHostToDevice, stream2);
         else
         cudaStat1 = cudaMemcpyAsync(&w_device[shift_init_Mi[id]-1-overlap_low], &w[shift_init_Mi[id]-1-overlap_low], (size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low)*sizeof(cuDoubleComplex)), cudaMemcpyHostToDevice, stream2);

         cudaFree(temp1_device);
         cudaFree(temp2_device);

         cublasSetStream(cubl_handle, NULL);

          // M(i,i) = eig_val(i+1);
          //M[(i*m)+i].x = eigen_val[i+1];

         for(int i = 0; i < m; i++)
         {
          // eig_val(i) = eig_val(i+1);
          eigen_val[i] = eigen_val[i+1];

         }

         
        cudaFree(eigen_vec_device);
        free(eigen_vec);

        cudaStat1=cudaDeviceSynchronize();

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
        //cubl_status = cublasDznrm2(cubl_handle, n_ham, v_device, 1, s_u_device);
        cubl_status = cublasZdotc(cubl_handle, (shift_end_Mi[id]-shift_init_Mi[id])+1, &v_device[shift_init_Mi[id]-1], 1, &v_device[shift_init_Mi[id]-1], 1, dot_temp_device);

        //MPI_Barrier(upt_comm);

        cudaStat1 = cudaMemcpy(dot_temp_1, dot_temp_device,(size_t)(sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);

        //MPI_Barrier(upt_comm);

        MPI_Allreduce(dot_temp_1, dot_temp_2, 1 , MPI_DOUBLE_COMPLEX, MPI_SUM, upt_comm);

        //MPI_Barrier(upt_comm);

        *s_u = sqrt(dot_temp_2[0].x);

        cudaStat1 = cudaMemcpy(s_u_device, s_u,(size_t)(sizeof(double)), cudaMemcpyHostToDevice);

        // s_v = eig_val(1) / (s_u * s_u); 
        *s_v = eigen_val[0]/((*s_u) * (*s_u));
        cudaStat1 = cudaMemcpy(s_v_device, s_v,(size_t)(sizeof(double)), cudaMemcpyHostToDevice);
        
        // u = v(:,1)/s_u;
        if(id == 0)
        {
        vct_div_slr_kernel<<<numBlock, threadPerBlock>>>(&v_device[shift_init_Mi[id]-1], &u_device[shift_init_Mi[id]-1], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high, *s_u);
        }
        else if(id != 0 && id != num_procs-1)
        {
        vct_div_slr_kernel<<<numBlock, threadPerBlock>>>(&v_device[shift_init_Mi[id]-1-overlap_low], &u_device[shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol, *s_u);
        }
        else
        {
        vct_div_slr_kernel<<<numBlock, threadPerBlock>>>(&v_device[shift_init_Mi[id]-1-overlap_low], &u_device[shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low, *s_u);
        }

        // r = (w(:,1) / s_u) - (s_v * u); 
        if(id == 0)
        {
        vct_div_slr_minus_scl_vct_kernel<<<numBlock, threadPerBlock>>>(&w_device[shift_init_Mi[id]-1], &u_device[shift_init_Mi[id]-1], &r_device[shift_init_Mi[id]-1], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high, *s_u, *s_v);
        }
        else if(id != 0 && id != num_procs-1)
        {
        vct_div_slr_minus_scl_vct_kernel<<<numBlock, threadPerBlock>>>(&w_device[shift_init_Mi[id]-1-overlap_low], &u_device[shift_init_Mi[id]-1-overlap_low], &r_device[shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol, *s_u, *s_v);
        }
        else
        {
        vct_div_slr_minus_scl_vct_kernel<<<numBlock, threadPerBlock>>>(&w_device[shift_init_Mi[id]-1-overlap_low], &u_device[shift_init_Mi[id]-1-overlap_low], &r_device[shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low, *s_u, *s_v);
        }

        // norm(r_bar);
        //cubl_status = cublasDznrm2(cubl_handle, n_ham, r_device, 1, norm_r_device);
        cubl_status = cublasZdotc(cubl_handle, (shift_end_Mi[id]-shift_init_Mi[id])+1, &r_device[shift_init_Mi[id]-1], 1, &r_device[shift_init_Mi[id]-1], 1, dot_temp_device);

        //MPI_Barrier(upt_comm);

        cudaStat1 = cudaMemcpy(dot_temp_1, dot_temp_device,(size_t)(sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);

        //MPI_Barrier(upt_comm);

        MPI_Allreduce(dot_temp_1, dot_temp_2, 1 , MPI_DOUBLE_COMPLEX, MPI_SUM, upt_comm);

        //MPI_Barrier(upt_comm);

        *norm_r = sqrt(dot_temp_2[0].x);

        cudaStat1 = cudaMemcpy(norm_r_device, norm_r,(size_t)(sizeof(double)), cudaMemcpyHostToDevice);
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

        cudaStat1 = cudaMalloc((void**)&temp1_device,n_ham*sizeof(cuDoubleComplex));
                    cudaMemset(temp1_device, 0, n_ham*sizeof(cuDoubleComplex));
        cudaStat2 = cudaMalloc((void**)&temp2_device,n_ham*sizeof(cuDoubleComplex));
                    cudaMemset(temp2_device, 0, n_ham*sizeof(cuDoubleComplex));
        if(cudaStat1 != cudaSuccess)
        {
           CLEANUP("Device malloc failed loc 14.1\n");
        }

        // for i = 1 : m_min
        for(int i = 0; i < jd_min_step; i++)
        {
          cublasSetStream(cubl_handle, stream1);
          //v_temp(:,i) = v * eig_mat(:,i);
          cubl_status = cublasZgemv(cubl_handle, CUBLAS_OP_N, n_ham, m, scalar1_device, v, n_ham, &eigen_vec_device[i*m], 1, scalar2_device, temp1_device, 1);

          cublasSetStream(cubl_handle, stream2);
          // w_temp(:,i) = w * eig_mat(:,i);
          cubl_status = cublasZgemv(cubl_handle, CUBLAS_OP_N, n_ham, m, scalar1_device, w, n_ham, &eigen_vec_device[i*m], 1, scalar2_device, temp2_device, 1);

          //cudaStat1 = cudaMemcpy(&temp_host[i*n_ham], temp_device,(size_t)(sizeof(cuDoubleComplex)*n_ham), cudaMemcpyDeviceToHost);
          if(id ==0)
          cudaStat1 = cudaMemcpyAsync(&temp_host_v[(i*n_ham)+shift_init_Mi[id]-1], &temp1_device[shift_init_Mi[id]-1],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high)*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost, stream1);
          else if(id != 0 && id != num_procs-1)
          cudaStat1 = cudaMemcpyAsync(&temp_host_v[(i*n_ham)+shift_init_Mi[id]-1-overlap_low], &temp1_device[shift_init_Mi[id]-1-overlap_low],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol)*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost, stream1);
          else
          cudaStat1 = cudaMemcpyAsync(&temp_host_v[(i*n_ham)+shift_init_Mi[id]-1-overlap_low], &temp1_device[shift_init_Mi[id]-1-overlap_low], (size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low)*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost, stream1);

         //cudaStat1 = cudaMemcpy(&temp_host[i*n_ham], temp_device,(size_t)(sizeof(cuDoubleComplex)*n_ham), cudaMemcpyDeviceToHost); 
         if(id ==0)
         cudaStat1 = cudaMemcpyAsync(&temp_host_w[(i*n_ham)+shift_init_Mi[id]-1], &temp2_device[shift_init_Mi[id]-1],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high)*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost, stream2);
         else if(id != 0 && id != num_procs-1)
         cudaStat1 = cudaMemcpyAsync(&temp_host_w[(i*n_ham)+shift_init_Mi[id]-1-overlap_low], &temp2_device[shift_init_Mi[id]-1-overlap_low],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol)*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost, stream2);
         else
         cudaStat1 = cudaMemcpyAsync(&temp_host_w[(i*n_ham)+shift_init_Mi[id]-1-overlap_low], &temp2_device[shift_init_Mi[id]-1-overlap_low], (size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low)*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost, stream2);
          // M(i,i) = eig_val(i);
          //M[(i*jd_min_step)+i].x = eigen_val[i];

         }

        // v = v_temp(:,(1:m_min));
       cudaStat1 = cudaMemcpy(v, temp_host_v,(size_t)(n_ham*jd_min_step*sizeof(cuDoubleComplex)), cudaMemcpyHostToHost);

       //cudaStat1 = cudaMemcpyAsync(v_device, v,(size_t)(sizeof(cuDoubleComplex)*n_ham), cudaMemcpyHostToDevice, stream1);
       if(id ==0)
       cudaStat1 = cudaMemcpyAsync(&v_device[shift_init_Mi[id]-1], &v[shift_init_Mi[id]-1],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high)*sizeof(cuDoubleComplex)), cudaMemcpyHostToDevice, stream1);
       else if(id != 0 && id != num_procs-1)
       cudaStat1 = cudaMemcpyAsync(&v_device[shift_init_Mi[id]-1-overlap_low], &v[shift_init_Mi[id]-1-overlap_low],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol)*sizeof(cuDoubleComplex)), cudaMemcpyHostToDevice, stream1);
       else
       cudaStat1 = cudaMemcpyAsync(&v_device[shift_init_Mi[id]-1-overlap_low], &v[shift_init_Mi[id]-1-overlap_low], (size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low)*sizeof(cuDoubleComplex)), cudaMemcpyHostToDevice, stream1);

         // w = w_temp(:,(1:m_min)); 
         cudaStat2 = cudaMemcpy(w, temp_host_w,(size_t)(n_ham*jd_min_step*sizeof(cuDoubleComplex)), cudaMemcpyHostToHost);

       cudaFree(temp1_device);
       cudaFree(temp2_device);

       cublasSetStream(cubl_handle, NULL);

       m = jd_min_step;

      } // end if(restart)

      cudaFree(M_device);
      cudaFree(eigen_vec_device);
      cudaFree(eigen_val_device);
      cudaFreeHost(M);
      free(eigen_val);
      free(eigen_vec);


      //cubl_status = cublasDznrm2(cubl_handle, n_ham, r_device, 1, norm_r_device);
      cubl_status = cublasZdotc(cubl_handle, (shift_end_Mi[id]-shift_init_Mi[id])+1, &r_device[shift_init_Mi[id]-1], 1, &r_device[shift_init_Mi[id]-1], 1, dot_temp_device);

      //MPI_Barrier(upt_comm);

      cudaStat1 = cudaMemcpy(dot_temp_1, dot_temp_device,(size_t)(sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);

      //MPI_Barrier(upt_comm);

      MPI_Allreduce(dot_temp_1, dot_temp_2, 1 , MPI_DOUBLE_COMPLEX, MPI_SUM, upt_comm);

      //MPI_Barrier(upt_comm);

      *norm_r = sqrt(dot_temp_2[0].x);

      cudaStat1 = cudaMemcpy(norm_r_device, norm_r,(size_t)(sizeof(double)), cudaMemcpyHostToDevice);
      if(cudaStat1 != cudaSuccess)
      {
         CLEANUP("cudaMemcpy failed loc 19\n");
      }
      
      if(*norm_r > 0.1 && flag ==1)
      {
        teta = shift;
      }
      else
      {
        teta = *s_v+shift;
        flag =0;
      }

     // A = (A - shift * I)
     tol_shift = shift-teta;
     shift_A_kernel<<<numBlock, threadPerBlock, 0, stream1>>>(valptr_device_real, rowptr_device_real, colptr_device_real, (shift_end_Mi[id]-shift_init_Mi[id])+1, tol_shift, shift_init_Mi[id]-1);

      //cudaDeviceSynchronize(); 
      //cudaStat1 = cudaMemcpy(u, u_device,(size_t)(sizeof(cuDoubleComplex)*n_ham), cudaMemcpyDeviceToHost);
      if(id ==0)
      cudaStat1 = cudaMemcpy(&u[shift_init_Mi[id]-1], &u_device[shift_init_Mi[id]-1],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high)*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);
      else if(id != 0 && id != num_procs-1)
      cudaStat1 = cudaMemcpy(&u[shift_init_Mi[id]-1-overlap_low], &u_device[shift_init_Mi[id]-1-overlap_low],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol)*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);
      else
      cudaStat1 = cudaMemcpy(&u[shift_init_Mi[id]-1-overlap_low], &u_device[shift_init_Mi[id]-1-overlap_low],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low)*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);

      cudaFree(u_device);


      if(k !=0)
      {
         Q_bar = (cuDoubleComplex *) malloc(sizeof(cuDoubleComplex)*n_ham*(k+1));
         cudaStat1 = cudaMalloc((void**)&Q_bar_device,sizeof(cuDoubleComplex)*n_ham*(k+1)); 
         memcpy (Q_bar, X_bar, sizeof(cuDoubleComplex)*n_ham*k);
         memcpy (&Q_bar[k*n_ham], u, sizeof(cuDoubleComplex)*n_ham);
         cudaStat2 = cudaMemcpy(Q_bar_device, Q_bar,(size_t)(sizeof(cuDoubleComplex)*n_ham*(k+1)), cudaMemcpyHostToDevice);
      }
      else
      {
        Q_bar = (cuDoubleComplex *) malloc(n_ham*sizeof(cuDoubleComplex));
        cudaStat1 = cudaMalloc((void**)&Q_bar_device,n_ham*sizeof(cuDoubleComplex)); 
        memcpy (Q_bar, u, n_ham*sizeof(cuDoubleComplex));
        cudaStat2 = cudaMemcpy(Q_bar_device, Q_bar,(size_t)(n_ham*sizeof(cuDoubleComplex)), cudaMemcpyHostToDevice);
      }


     cudaFree(v_device);
     cudaFree(w_device);
     //cudaFree(v_buffer_device);
     //cudaFree(w_buffer_device);

     cudaStat1=cudaDeviceSynchronize();
     MPI_Barrier(upt_comm); 

     gmres(valptr_device_real, rowptr_device_real, colptr_device_real, valptr_device_img, rowptr_device_img, colptr_device_img, n_ham, size_mat_real, r_device, ls_tol, ls_restart, ls_maxit, Q_bar_device, t_device, numBlock, threadPerBlock, coop, repeat, numBlocksMul, k, cusp_handle, cusp_descra, cubl_handle, ls_counter, shift_init_Mi, shift_end_Mi, overlap_high, overlap_low, num_procs, id,  upt_comm);

     cudaStat1=cudaDeviceSynchronize();
     MPI_Barrier(upt_comm); 

     free(Q_bar);
     cudaFree(Q_bar_device);

     if(k!=0) // MGS after linear solver
     {
       cudaStat1 = cudaMalloc((void**)&X_bar_device,n_ham*k*sizeof(cuDoubleComplex));
       if(cudaStat1 != cudaSuccess)
       {
         CLEANUP("cudaMemcpy device malloc falied loc 21\n");
       }

       cudaStat1   = cudaMemcpy(X_bar_device, X_bar,(size_t)(sizeof(cuDoubleComplex)*n_ham*k), cudaMemcpyHostToDevice);

       for(int i = 0; i < k; i++)
       {
         // y = X_bar(:,i)' * t;
         cubl_status = cublasZdotc(cubl_handle, (shift_end_Mi[id]-shift_init_Mi[id])+1, &X_bar_device[i*n_ham+shift_init_Mi[id]-1], 1, &t_device[shift_init_Mi[id]-1], 1, y_device); 

         //MPI_Barrier(upt_comm); 

         cudaStat1 = cudaMemcpy(dot_temp, y_device,(size_t)(sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);

         //MPI_Barrier(upt_comm);

         MPI_Allreduce(dot_temp, y, 1 , MPI_DOUBLE_COMPLEX, MPI_SUM, upt_comm);

         cudaStat1 = cudaMemcpy(y_device, y,(size_t)(sizeof(cuDoubleComplex)), cudaMemcpyHostToDevice);
 

         // t = t - y * X_bar(:,i);
         if(id == 0)
         {
         vct1_sub_mul_vct_kernel<<<numBlock, threadPerBlock>>>( &t_device[shift_init_Mi[id]-1], &X_bar_device[i*n_ham+shift_init_Mi[id]-1], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high, y_device);
         }
         else if(id != 0 && id != num_procs-1)
         { 
         vct1_sub_mul_vct_kernel<<<numBlock, threadPerBlock>>>( &t_device[shift_init_Mi[id]-1-overlap_low], &X_bar_device[i*n_ham+shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol, y_device);
         }
         else
         {
         vct1_sub_mul_vct_kernel<<<numBlock, threadPerBlock>>>( &t_device[shift_init_Mi[id]-1-overlap_low], &X_bar_device[i*n_ham+shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low, y_device);
         }

       }

         cudaFree(X_bar_device);
      }

     tol_shift = teta-shift;
     shift_A_kernel<<<numBlock, threadPerBlock, 0, stream2>>>(valptr_device_real, rowptr_device_real, colptr_device_real, (shift_end_Mi[id]-shift_init_Mi[id])+1, tol_shift, shift_init_Mi[id]-1);

     if (num_procs > 1){
     for(int rank =0; rank <=num_procs-2; rank++){
     if(id == rank+1)    MPI_Irecv(&t_device[shift_end_Mi[rank]-overlap_low], overlap_low,MPI_DOUBLE_COMPLEX,rank, 111+rank,upt_comm,&reqs[rank+1]);
     if(id == rank)      MPI_Isend(&t_device[shift_end_Mi[rank]-overlap_low], overlap_low,MPI_DOUBLE_COMPLEX,rank+1, 111+rank,upt_comm, &reqs[rank+1]);

     if (id == rank+1) MPI_Wait(&reqs[rank+1], &status);
     if (id == rank)   MPI_Wait(&reqs[rank+1], &status);
     }
     }

     //MPI_Barrier(upt_comm);
//printf("first t_device send done\n");

     if (num_procs > 1){
     for(int rank =0; rank <=num_procs-2; rank++){
     if(id == rank)    MPI_Irecv(&t_device[shift_init_Mi[rank+1]-1],overlap_high, MPI_DOUBLE_COMPLEX,rank+1,8+rank,upt_comm, &reqs[rank+1]);
     if(id == rank+1)  MPI_Isend(&t_device[shift_init_Mi[rank+1]-1],overlap_high, MPI_DOUBLE_COMPLEX,rank,8+rank,upt_comm, &reqs[rank+1]);

     if (id == rank+1) MPI_Wait(&reqs[rank+1], &status);
     if (id == rank)   MPI_Wait(&reqs[rank+1], &status);
     }
     }

count = count+1;
//printf("count = %d\n", count);
//if(count == 30)
//exit(0);
  //MPI_Barrier(upt_comm); 
  //cudaStat1=cudaDeviceSynchronize();
} // end while(k < num_ev)

printf("Total JD count = %d\n", count);
printf("Total MxV iterations = %d\n", *ls_counter+count+(*ls_counter/ls_restart)); 


}// end of JD



















void gmres(float *valptr_device_real, int *rowptr_device_real, int *colptr_device_real, float *valptr_device_img, int *rowptr_device_img, int *colptr_device_img, int n_ham, int size_mat_real, cuDoubleComplex * r_device, double ls_tol, int restart, int maxit, cuDoubleComplex *Q_bar_device, cuDoubleComplex *t_device, int numBlock, int threadPerBlock, int coop, int repeat, int numBlocksMul, int k, cusparseHandle_t cusp_handle_ls, cusparseMatDescr_t cusp_descra_ls, cublasHandle_t cubl_handle_ls, int *ls_counter, int * shift_init_Mi, int * shift_end_Mi, int overlap_high, int overlap_low, int num_procs, int id,  MPI_Comm upt_comm)
{

  //cuDoubleComplex *r_device_ls;
  cuDoubleComplex *X0_device, *Xmin_device, *dot_device, *P, *P_device, *v_ls, *buffer_device;
  cuDoubleComplex *rq_device, *h_min, *g, cudoublecomplex_temp1, cudoublecomplex_temp2, cudoublecomplex_temp, sqrt_cudoublecomplex_temp;
  cuDoubleComplex *temp_device, *temp2_device, *test, *temp1_device;
  double *beta_device, *norm_r_device_ls, *norm_r_ls, *beta, *norm_temp;
  cuDoubleComplex *v_ls_device, *w_ls_device, *g_device, *g_min, *minimizer_device, *h_device, *h; 
  cuDoubleComplex *scalar1_device_ls, *scalar2_device_ls;
  cuDoubleComplex *scalar1_ls, *scalar2_ls, *dot, *dot_temp_1, *dot_temp_2, *dot_temp, *dot_temp_device_ls;
  lapack_int *jpvt;
  lapack_int error_code=0;
  cuDoubleComplex coso, sino, scalar_1, scalar_2;

  MPI_Request * reqs_ls;
  MPI_Status status_ls;

  reqs_ls=(MPI_Request *) malloc(sizeof(MPI_Request)*(num_procs+1));

  cudaStream_t stream5, stream6, stream7, stream8;
  cudaStreamCreate(&stream5);
  cudaStreamCreate(&stream6);
  cudaStreamCreate(&stream7);
  cudaStreamCreate(&stream8);

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
  cublasStatus_t cubl_status_ls;

  int overlap_tol = overlap_high + overlap_low;

  cudaStat1 = cudaMalloc((void**)&X0_device,sizeof(cuDoubleComplex)*n_ham);
  cudaStat2 = cudaMalloc((void**)&Xmin_device,sizeof(cuDoubleComplex)*n_ham); 
  cudaStat3 = cudaMemset(X0_device, 0, n_ham*sizeof(cuDoubleComplex));
  cudaStat4 = cudaMemset(Xmin_device, 0, n_ham*sizeof(cuDoubleComplex));

  cudaStat5 = cudaMalloc((void**)&dot_device,sizeof(cuDoubleComplex)); 

  cudaStat6 = cudaMalloc((void**)&temp_device,sizeof(cuDoubleComplex)*n_ham);
              cudaMemset(temp_device, 0, n_ham*sizeof(cuDoubleComplex));  

  cudaStat7 = cudaMalloc((void**)&temp2_device,sizeof(cuDoubleComplex)*n_ham);
              cudaMemset(temp2_device, 0, n_ham*sizeof(cuDoubleComplex));

  cudaStat8 = cudaMalloc((void**)&rq_device,sizeof(cuDoubleComplex)*n_ham);
              cudaMemset(rq_device, 0, n_ham*sizeof(cuDoubleComplex));

  cudaStat9 = cudaMalloc((void**)&beta_device,sizeof(double));

  cudaStat10 = cudaMalloc((void**)&v_ls_device,sizeof(cuDoubleComplex)*n_ham);
               cudaMemset(v_ls_device, 0, sizeof(cuDoubleComplex)*n_ham);

  cudaStat11 = cudaMalloc((void**)&w_ls_device,sizeof(cuDoubleComplex)*n_ham);
               cudaMemset(w_ls_device, 0, sizeof(cuDoubleComplex)*n_ham);

  cudaStat12 = cudaMalloc((void**)&h_device,sizeof(cuDoubleComplex)*(restart+1)*restart);
               cudaMemset(h_device, 0, sizeof(cuDoubleComplex)*(restart+1)*restart);

  cudaStat13 = cudaMalloc((void**)&g_device,sizeof(cuDoubleComplex)*(restart+1));
               cudaMemset(g_device, 0, sizeof(cuDoubleComplex)*(restart+1));

  cudaStat14 = cudaMalloc((void**)&minimizer_device,sizeof(cuDoubleComplex)*restart);
               cudaMemset(minimizer_device, 0, sizeof(cuDoubleComplex)*restart);

  cudaStat15 = cudaMalloc((void**)&norm_r_device_ls,sizeof(double));
  cudaStat16 = cudaMalloc((void**)&scalar1_device_ls,sizeof(cuDoubleComplex)); 
  cudaStat17 = cudaMalloc((void**)&scalar2_device_ls,sizeof(cuDoubleComplex));  
  cudaStat18 = cudaMalloc((void**)&temp1_device,n_ham*sizeof(cuDoubleComplex));
               cudaMemset(temp1_device, 0, n_ham*sizeof(cuDoubleComplex));
               cudaMalloc((void**)&dot_temp_device_ls,sizeof(cuDoubleComplex));
               
               cudaMalloc((void**)&buffer_device,n_ham*sizeof(cuDoubleComplex));
               cudaMemset(buffer_device, 0, n_ham*sizeof(cuDoubleComplex));
               cudaMalloc((void**)&dot_temp_device_ls,sizeof(cuDoubleComplex));

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
    CLEANUP_LS("cudaMemcpy failed\n");
  }


  // r_ls = -r
  //vct1_neg_asg_vct2_kernel<<<numBlock, threadPerBlock>>>(r_device, r_device, n_ham);
  if(id == 0)
  {
  vct1_neg_asg_vct2_kernel<<<numBlock, threadPerBlock>>>( &r_device[shift_init_Mi[id]-1], &r_device[shift_init_Mi[id]-1], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high);
  }
  else if(id != 0 && id != num_procs-1)
  { 
  vct1_neg_asg_vct2_kernel<<<numBlock, threadPerBlock>>>( &r_device[shift_init_Mi[id]-1-overlap_low], &r_device[shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol);
  }
  else
  {
  vct1_neg_asg_vct2_kernel<<<numBlock, threadPerBlock>>>( &r_device[shift_init_Mi[id]-1-overlap_low], &r_device[shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low);
  }

  double norm_min =10000000000;
  double rcond;
  lapack_int *eff_rank_h_min;


  dot = (cuDoubleComplex *) malloc(sizeof(cuDoubleComplex));
  dot_temp = (cuDoubleComplex *) malloc(sizeof(cuDoubleComplex));
  dot_temp_1 = (cuDoubleComplex *) malloc(sizeof(cuDoubleComplex));
  dot_temp_2 = (cuDoubleComplex *) malloc(sizeof(cuDoubleComplex));
  norm_temp = (double *) malloc(sizeof(double));
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
  cudaMallocHost((void**)&v_ls, (restart+1)*n_ham*sizeof(cuDoubleComplex)); 
  cudaMemset(v_ls, 0, n_ham*(restart+1)*sizeof(cuDoubleComplex));

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
  //cubl_status_ls = cublasDznrm2(cubl_handle_ls, n_ham, r_device, 1, norm_r_device_ls);
  cubl_status_ls = cublasZdotc(cubl_handle_ls, (shift_end_Mi[id]-shift_init_Mi[id])+1, &r_device[shift_init_Mi[id]-1], 1, &r_device[shift_init_Mi[id]-1], 1, dot_temp_device_ls);

  //MPI_Barrier(upt_comm);

  cudaStat1 = cudaMemcpy(dot_temp_1, dot_temp_device_ls,(size_t)(sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);

  //MPI_Barrier(upt_comm);

  MPI_Allreduce(dot_temp_1, dot_temp_2, 1 , MPI_DOUBLE_COMPLEX, MPI_SUM, upt_comm);

  //MPI_Barrier(upt_comm);

  *norm_r_ls = sqrt(dot_temp_2[0].x);

  cudaStat1 = cudaMemcpy(norm_r_device_ls, norm_r_ls,(size_t)(sizeof(double)), cudaMemcpyHostToDevice); 

  ls_tol = (*norm_r_ls) * ls_tol;

  int restart_count = 0;

  while(restart_count < maxit)
  {
    // temp = x0;
    //cudaStat1 = cudaMemcpy(temp_device, X0_device,(size_t)(n_ham*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToDevice);
    if(id ==0)
    cudaStat1 = cudaMemcpy(&temp_device[shift_init_Mi[id]-1], &X0_device[shift_init_Mi[id]-1],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high)*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToDevice);
    else if(id != 0 && id != num_procs-1)
    cudaStat1 = cudaMemcpy(&temp_device[shift_init_Mi[id]-1-overlap_low], &X0_device[shift_init_Mi[id]-1-overlap_low],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol)*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToDevice);
    else
    cudaStat1 = cudaMemcpy(&temp_device[shift_init_Mi[id]-1-overlap_low], &X0_device[shift_init_Mi[id]-1-overlap_low],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low)*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToDevice);


    for(int i = 0; i < k+1; i++)
    {
      cubl_status_ls = cublasZdotc(cubl_handle_ls, (shift_end_Mi[id]-shift_init_Mi[id])+1, &Q_bar_device[i*n_ham+shift_init_Mi[id]-1], 1, &X0_device[shift_init_Mi[id]-1], 1, dot_device);

       //MPI_Barrier(upt_comm); 

       cudaStat1 = cudaMemcpy(dot_temp, dot_device,(size_t)(sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);

       //MPI_Barrier(upt_comm);

       MPI_Allreduce(dot_temp, dot, 1 , MPI_DOUBLE_COMPLEX, MPI_SUM, upt_comm);

       cudaStat1 = cudaMemcpy(dot_device, dot,(size_t)(sizeof(cuDoubleComplex)), cudaMemcpyHostToDevice);

       //MPI_Barrier(upt_comm);


       if(id == 0)
       {
       vct1_sub_mul_vct_kernel<<<numBlock, threadPerBlock>>>( &temp_device[shift_init_Mi[id]-1], &Q_bar_device[i*n_ham+shift_init_Mi[id]-1], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high, dot_device);
       }
       else if(id != 0 && id != num_procs-1)
       { 
       vct1_sub_mul_vct_kernel<<<numBlock, threadPerBlock>>>( &temp_device[shift_init_Mi[id]-1-overlap_low], &Q_bar_device[i*n_ham+shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol, dot_device);
       }
       else
       {
       vct1_sub_mul_vct_kernel<<<numBlock, threadPerBlock>>>( &temp_device[shift_init_Mi[id]-1-overlap_low], &Q_bar_device[i*n_ham+shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low, dot_device);
       }
    }

     spmv_csr_hybrid_kernel<<<numBlocksMul, BLOCK_SIZE_MUL, 0, stream5>>>((shift_end_Mi[id]-shift_init_Mi[id])+1, rowptr_device_real, colptr_device_real, valptr_device_real, temp_device, temp2_device, repeat, coop, shift_init_Mi[id]-1);

     spmv_csr_hybrid_kernel<<<numBlock, BLOCK_SIZE, 0, stream6>>>((shift_end_Mi[id]-shift_init_Mi[id])+1, rowptr_device_img, colptr_device_img, valptr_device_img, temp_device, temp1_device, shift_init_Mi[id]-1);

     cudaStat1=cudaDeviceSynchronize();

     vct_pls_scl_mul_vct_kernel_jd<<<numBlock, threadPerBlock>>>(&temp2_device[shift_init_Mi[id]-1], &temp1_device[shift_init_Mi[id]-1], jcmpx, (shift_end_Mi[id]-shift_init_Mi[id])+1);

     cudaStat1=cudaDeviceSynchronize();
     //MPI_Barrier(upt_comm);

     if (num_procs > 1){
     for(int rank =0; rank <=num_procs-2; rank++){
     if(id == rank+1)    MPI_Irecv(&temp2_device[shift_end_Mi[rank]-overlap_low], overlap_low, MPI_DOUBLE_COMPLEX, rank, 411+rank, upt_comm, &reqs_ls[rank+1]);
     if(id == rank)      MPI_Isend(&temp2_device[shift_end_Mi[rank]-overlap_low], overlap_low, MPI_DOUBLE_COMPLEX, rank+1, 411+rank, upt_comm, &reqs_ls[rank+1]);

     if (id == rank+1) MPI_Wait(&reqs_ls[rank+1], &status_ls);
     if (id == rank)   MPI_Wait(&reqs_ls[rank+1], &status_ls);
     }
     }


     if (num_procs > 1){
     for(int rank =0; rank <=num_procs-2; rank++){
     if(id == rank)    MPI_Irecv(&temp2_device[shift_init_Mi[rank+1]-1],overlap_high, MPI_DOUBLE_COMPLEX,rank+1,450+rank,upt_comm, &reqs_ls[rank+1]);
     if(id == rank+1)  MPI_Isend(&temp2_device[shift_init_Mi[rank+1]-1],overlap_high, MPI_DOUBLE_COMPLEX,rank,450+rank,upt_comm, &reqs_ls[rank+1]);

     if (id == rank+1) MPI_Wait(&reqs_ls[rank+1], &status_ls);
     if (id == rank)   MPI_Wait(&reqs_ls[rank+1], &status_ls);
     }
     }

     MPI_Barrier(upt_comm);

   // temp = temp2;
   //cudaStat1 = cudaMemcpyAsync(temp_device, temp2_device,(size_t)(n_ham*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToDevice, stream5);
   if(id ==0)
    cudaStat1 = cudaMemcpy(&temp_device[shift_init_Mi[id]-1], &temp2_device[shift_init_Mi[id]-1],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high)*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToDevice);
    else if(id != 0 && id != num_procs-1)
    cudaStat1 = cudaMemcpy(&temp_device[shift_init_Mi[id]-1-overlap_low], &temp2_device[shift_init_Mi[id]-1-overlap_low],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol)*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToDevice);
    else
    cudaStat1 = cudaMemcpy(&temp_device[shift_init_Mi[id]-1-overlap_low], &temp2_device[shift_init_Mi[id]-1-overlap_low],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low)*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToDevice);
   

   for(int i = 0; i < k+1; i++)
    {
       // temp = temp - (dot(u(:,t),temp2)) * u(:,t);
       cubl_status_ls = cublasZdotc(cubl_handle_ls, (shift_end_Mi[id]-shift_init_Mi[id])+1, &Q_bar_device[i*n_ham+shift_init_Mi[id]-1], 1, &temp2_device[shift_init_Mi[id]-1], 1, dot_device);

       //MPI_Barrier(upt_comm); 

       cudaStat1 = cudaMemcpy(dot_temp, dot_device,(size_t)(sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);

       //MPI_Barrier(upt_comm);

       MPI_Allreduce(dot_temp, dot, 1 , MPI_DOUBLE_COMPLEX, MPI_SUM, upt_comm);

       cudaStat1 = cudaMemcpy(dot_device, dot,(size_t)(sizeof(cuDoubleComplex)), cudaMemcpyHostToDevice);

       //MPI_Barrier(upt_comm);


       if(id == 0)
       {
       vct1_sub_mul_vct_kernel<<<numBlock, threadPerBlock>>>( &temp_device[shift_init_Mi[id]-1], &Q_bar_device[i*n_ham+shift_init_Mi[id]-1], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high, dot_device);
       }
       else if(id != 0 && id != num_procs-1)
       { 
       vct1_sub_mul_vct_kernel<<<numBlock, threadPerBlock>>>( &temp_device[shift_init_Mi[id]-1-overlap_low], &Q_bar_device[i*n_ham+shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol, dot_device);
       }
       else
       {
       vct1_sub_mul_vct_kernel<<<numBlock, threadPerBlock>>>( &temp_device[shift_init_Mi[id]-1-overlap_low], &Q_bar_device[i*n_ham+shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low, dot_device);
       }
    }

    //rq = b - temp2;     here b = r_device_ls
    if(id == 0)
    {
    vct1_sub_vct2_asg_vct3_kernel<<<numBlock, threadPerBlock>>>( &r_device[shift_init_Mi[id]-1], &temp_device[shift_init_Mi[id]-1], &rq_device[shift_init_Mi[id]-1], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high);
    }
    else if(id != 0 && id != num_procs-1)
    { 
    vct1_sub_vct2_asg_vct3_kernel<<<numBlock, threadPerBlock>>>( &r_device[shift_init_Mi[id]-1-overlap_low], &temp_device[shift_init_Mi[id]-1-overlap_low], &rq_device[shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol);
    }
    else
    {
    vct1_sub_vct2_asg_vct3_kernel<<<numBlock, threadPerBlock>>>( &r_device[shift_init_Mi[id]-1-overlap_low], &temp_device[shift_init_Mi[id]-1-overlap_low], &rq_device[shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low);
    }
   
    // beta=norm(rq);
    //cubl_status_ls = cublasDznrm2(cubl_handle_ls, n_ham, rq_device, 1, beta_device);
    cubl_status_ls = cublasZdotc(cubl_handle_ls, (shift_end_Mi[id]-shift_init_Mi[id])+1, &rq_device[shift_init_Mi[id]-1], 1, &rq_device[shift_init_Mi[id]-1], 1, dot_temp_device_ls);

    //MPI_Barrier(upt_comm);

    cudaStat1 = cudaMemcpy(dot_temp_1, dot_temp_device_ls,(size_t)(sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);

    //MPI_Barrier(upt_comm);

    MPI_Allreduce(dot_temp_1, dot_temp_2, 1 , MPI_DOUBLE_COMPLEX, MPI_SUM, upt_comm);

    //MPI_Barrier(upt_comm);

    *beta = sqrt(dot_temp_2[0].x);

    cudaStat1 = cudaMemcpy(beta_device, beta,(size_t)(sizeof(double)), cudaMemcpyHostToDevice); 

    // v(:,1)=rq/beta; 
    if(id == 0)
    {
    vct_div_slr_kernel<<<numBlock, threadPerBlock>>>( &rq_device[shift_init_Mi[id]-1], &v_ls_device[shift_init_Mi[id]-1], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high, *beta);
    }
    else if(id != 0 && id != num_procs-1)
    { 
    vct_div_slr_kernel<<<numBlock, threadPerBlock>>>( &rq_device[shift_init_Mi[id]-1-overlap_low], &v_ls_device[shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol, *beta);
    }
    else
    {
    vct_div_slr_kernel<<<numBlock, threadPerBlock>>>( &rq_device[shift_init_Mi[id]-1-overlap_low], &v_ls_device[shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low, *beta);
    }
    
    if(id ==0)
    cudaStat1 = cudaMemcpyAsync(&v_ls[shift_init_Mi[id]-1], &v_ls_device[shift_init_Mi[id]-1],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high)*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost, stream6);
    else if(id != 0 && id != num_procs-1)
    cudaStat1 = cudaMemcpyAsync(&v_ls[shift_init_Mi[id]-1-overlap_low], &v_ls_device[shift_init_Mi[id]-1-overlap_low],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol)*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost, stream6);
    else
    cudaStat1 = cudaMemcpyAsync(&v_ls[shift_init_Mi[id]-1-overlap_low], &v_ls_device[shift_init_Mi[id]-1-overlap_low],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low)*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost, stream6);

    cudaStat1 = cudaMemset(h_device, 0, sizeof(cuDoubleComplex)*restart*(restart+1));
    

//printf("inside gmres pt 3\n"); 

    for(int j = 0; j < restart; j++)
    {
      // temp = v(:,j);
      //cudaStat1 = cudaMemcpy(temp_device, &v_ls[j*n_ham],(size_t)(n_ham*sizeof(cuDoubleComplex)), cudaMemcpyHostToDevice);
      if(id ==0)
      cudaStat1 = cudaMemcpyAsync(&temp_device[shift_init_Mi[id]-1], &v_ls_device[(shift_init_Mi[id]-1)],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high)*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToDevice, stream7);
      else if(id != 0 && id != num_procs-1)
      cudaStat1 = cudaMemcpyAsync(&temp_device[shift_init_Mi[id]-1-overlap_low], &v_ls_device[(shift_init_Mi[id]-1-overlap_low)],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol)*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToDevice, stream7);
      else
      cudaStat1 = cudaMemcpyAsync(&temp_device[shift_init_Mi[id]-1-overlap_low], &v_ls_device[(shift_init_Mi[id]-1-overlap_low)],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low)*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToDevice, stream7);

      //cudaStat1=cudaDeviceSynchronize();

      for(int i = 0; i < k+1; i++)
      {
       cublasSetStream(cubl_handle_ls, stream7);

       // temp = temp - (dot(Q_bar(:,i),v(:,j))) * Q_bar(:,i); norm(Q_bar) == 1
       cubl_status_ls = cublasZdotc(cubl_handle_ls, (shift_end_Mi[id]-shift_init_Mi[id])+1, &Q_bar_device[i*n_ham+shift_init_Mi[id]-1], 1, &v_ls_device[shift_init_Mi[id]-1], 1, dot_device);

       //MPI_Barrier(upt_comm); 

       cudaStat1 = cudaMemcpyAsync(dot_temp, dot_device,(size_t)(sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost, stream7);

       //MPI_Barrier(upt_comm);
       //cudaStat1=cudaDeviceSynchronize();

       MPI_Allreduce(dot_temp, dot, 1 , MPI_DOUBLE_COMPLEX, MPI_SUM, upt_comm);

       //MPI_Barrier(upt_comm);

       cudaStat1 = cudaMemcpyAsync(dot_device, dot,(size_t)(sizeof(cuDoubleComplex)), cudaMemcpyHostToDevice, stream7);


       if(id == 0)
       {
       vct1_sub_mul_vct_kernel<<<numBlock, threadPerBlock, 0, stream7>>>( &temp_device[shift_init_Mi[id]-1], &Q_bar_device[i*n_ham+shift_init_Mi[id]-1], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high, dot_device);
       }
       else if(id != 0 && id != num_procs-1)
       { 
       vct1_sub_mul_vct_kernel<<<numBlock, threadPerBlock, 0, stream7>>>( &temp_device[shift_init_Mi[id]-1-overlap_low], &Q_bar_device[i*n_ham+shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol, dot_device);
       }
       else
       {
       vct1_sub_mul_vct_kernel<<<numBlock, threadPerBlock, 0, stream7>>>( &temp_device[shift_init_Mi[id]-1-overlap_low], &Q_bar_device[i*n_ham+shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low, dot_device);
       }
      }


      // temp2=A*temp;
      // w(:,j)=temp2;
      //cusp_status_ls= cusparseZcsrmv(cusp_handle_ls,CUSPARSE_OPERATION_NON_TRANSPOSE, n_ham, n_ham, size_mat_real, &scalar_1, cusp_descra_ls, valptr_device, rowptr_device, colptr_device, temp_device, &scalar_2, temp2_device);
     cudaStat1=cudaDeviceSynchronize();

     spmv_csr_hybrid_kernel<<<numBlocksMul, BLOCK_SIZE_MUL, 0, stream5>>>((shift_end_Mi[id]-shift_init_Mi[id])+1, rowptr_device_real, colptr_device_real, valptr_device_real, temp_device, temp2_device, repeat, coop, shift_init_Mi[id]-1);

     //cudaStat1=cudaDeviceSynchronize();

     spmv_csr_hybrid_kernel<<<numBlock, BLOCK_SIZE, 0, stream6>>>((shift_end_Mi[id]-shift_init_Mi[id])+1, rowptr_device_img, colptr_device_img, valptr_device_img, temp_device, temp1_device, shift_init_Mi[id]-1);

     cudaStat1=cudaDeviceSynchronize();

     vct_pls_scl_mul_vct_kernel_jd<<<numBlock, threadPerBlock>>>(&temp2_device[shift_init_Mi[id]-1], &temp1_device[shift_init_Mi[id]-1], jcmpx, (shift_end_Mi[id]-shift_init_Mi[id])+1);

     cudaStat1=cudaDeviceSynchronize();

     if (num_procs > 1){
     for(int rank =0; rank <=num_procs-2; rank++){
     if(id == rank+1)    MPI_Irecv(&temp2_device[shift_end_Mi[rank]-overlap_low], overlap_low,MPI_DOUBLE_COMPLEX,rank, 311+rank,upt_comm,&reqs_ls[rank+1]);
     if(id == rank)      MPI_Isend(&temp2_device[shift_end_Mi[rank]-overlap_low], overlap_low,MPI_DOUBLE_COMPLEX,rank+1, 311+rank,upt_comm, &reqs_ls[rank+1]);

     if (id == rank+1) MPI_Wait(&reqs_ls[rank+1], &status_ls);
     if (id == rank)   MPI_Wait(&reqs_ls[rank+1], &status_ls);
     }
     }

     //MPI_Barrier(upt_comm);


     if (num_procs > 1){
     for(int rank =0; rank <=num_procs-2; rank++){
     if(id == rank)    MPI_Irecv(&temp2_device[shift_init_Mi[rank+1]-1],overlap_high, MPI_DOUBLE_COMPLEX,rank+1,350+rank,upt_comm, &reqs_ls[rank+1]);
     if(id == rank+1)  MPI_Isend(&temp2_device[shift_init_Mi[rank+1]-1],overlap_high, MPI_DOUBLE_COMPLEX,rank,350+rank,upt_comm, &reqs_ls[rank+1]);

     if (id == rank+1) MPI_Wait(&reqs_ls[rank+1], &status_ls);
     if (id == rank)   MPI_Wait(&reqs_ls[rank+1], &status_ls);
     }
     }

     MPI_Barrier(upt_comm);

     //cudaStat1=cudaDeviceSynchronize();

     //cudaStat1 = cudaMemcpy(w_ls_device, temp2_device,(size_t)(n_ham*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToDevice);
     if(id ==0)
     cudaStat1 = cudaMemcpy(&w_ls_device[shift_init_Mi[id]-1], &temp2_device[shift_init_Mi[id]-1],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high)*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToDevice);
     else if(id != 0 && id != num_procs-1)
     cudaStat1 = cudaMemcpy(&w_ls_device[shift_init_Mi[id]-1-overlap_low], &temp2_device[shift_init_Mi[id]-1-overlap_low],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol)*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToDevice);
     else
     cudaStat1 = cudaMemcpy(&w_ls_device[shift_init_Mi[id]-1-overlap_low], &temp2_device[shift_init_Mi[id]-1-overlap_low],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low)*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToDevice);

//printf("inside gmres pt 4\n");

     for(int i = 0; i < k+1; i++)
      {
       // temp = temp - (dot(Q_bar(:,i),w(:,j))) * Q_bar(:,i); norm(Q_bar) == 1
       cubl_status_ls = cublasZdotc(cubl_handle_ls, (shift_end_Mi[id]-shift_init_Mi[id])+1, &Q_bar_device[i*n_ham+shift_init_Mi[id]-1], 1, &temp2_device[shift_init_Mi[id]-1], 1, dot_device);

       //MPI_Barrier(upt_comm); 

       cudaStat1 = cudaMemcpy(dot_temp, dot_device,(size_t)(sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);

       //MPI_Barrier(upt_comm);

       MPI_Allreduce(dot_temp, dot, 1 , MPI_DOUBLE_COMPLEX, MPI_SUM, upt_comm);

       cudaStat1 = cudaMemcpy(dot_device, dot,(size_t)(sizeof(cuDoubleComplex)), cudaMemcpyHostToDevice);

       //MPI_Barrier(upt_comm);

       if(id == 0)
       {
       vct1_sub_mul_vct_kernel<<<numBlock, threadPerBlock>>>( &w_ls_device[shift_init_Mi[id]-1], &Q_bar_device[i*n_ham+shift_init_Mi[id]-1], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high, dot_device);
       }
       else if(id != 0 && id != num_procs-1)
       { 
       vct1_sub_mul_vct_kernel<<<numBlock, threadPerBlock>>>( &w_ls_device[shift_init_Mi[id]-1-overlap_low], &Q_bar_device[i*n_ham+shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol, dot_device);
       }
       else
       {
       vct1_sub_mul_vct_kernel<<<numBlock, threadPerBlock>>>( &w_ls_device[shift_init_Mi[id]-1-overlap_low], &Q_bar_device[i*n_ham+shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low, dot_device);
       }
      }

      // w(:,j)=temp2;
      //cudaStat1 = cudaMemcpy(w_ls_device, temp2_device,(size_t)(n_ham*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToDevice);


    if(id ==0)
    cudaStat1 = cudaMemcpy(&buffer_device[shift_init_Mi[id]-1], &v_ls[shift_init_Mi[id]-1],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high)*sizeof(cuDoubleComplex)), cudaMemcpyHostToDevice);
    else if(id != 0 && id != num_procs-1)
    cudaStat1 = cudaMemcpy(&buffer_device[shift_init_Mi[id]-1-overlap_low], &v_ls[shift_init_Mi[id]-1-overlap_low],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol)*sizeof(cuDoubleComplex)), cudaMemcpyHostToDevice);
    else
    cudaStat1 = cudaMemcpy(&buffer_device[shift_init_Mi[id]-1-overlap_low], &v_ls[shift_init_Mi[id]-1-overlap_low],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low)*sizeof(cuDoubleComplex)), cudaMemcpyHostToDevice);


      for(int i = 0; i <= j ; i++)
      {
        //cudaStat1 = cudaMemcpy(v_ls_device, &v_ls[i*n_ham],(size_t)(n_ham*sizeof(cuDoubleComplex)), cudaMemcpyHostToDevice);
        if(id ==0)
        cudaStat1 = cudaMemcpy(&v_ls_device[shift_init_Mi[id]-1], &buffer_device[shift_init_Mi[id]-1],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high)*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToDevice);
        else if(id != 0 && id != num_procs-1)
        cudaStat1 = cudaMemcpy(&v_ls_device[shift_init_Mi[id]-1-overlap_low], &buffer_device[shift_init_Mi[id]-1-overlap_low],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol)*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToDevice);
        else
        cudaStat1 = cudaMemcpy(&v_ls_device[shift_init_Mi[id]-1-overlap_low], &buffer_device[shift_init_Mi[id]-1-overlap_low],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low)*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToDevice);

        cudaStat1=cudaDeviceSynchronize();

        if(i+1 != j+1){
        if(id ==0)
        cudaStat2 = cudaMemcpyAsync(&buffer_device[shift_init_Mi[id]-1], &v_ls[(i+1)*n_ham+shift_init_Mi[id]-1],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high)*sizeof(cuDoubleComplex)), cudaMemcpyHostToDevice, stream5);
        else if(id != 0 && id != num_procs-1)
        cudaStat2 = cudaMemcpyAsync(&buffer_device[shift_init_Mi[id]-1-overlap_low], &v_ls[(i+1)*n_ham+shift_init_Mi[id]-1-overlap_low],  (size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol)*sizeof(cuDoubleComplex)), cudaMemcpyHostToDevice, stream5);
        else
        cudaStat2 = cudaMemcpyAsync(&buffer_device[shift_init_Mi[id]-1-overlap_low], &v_ls[(i+1)*n_ham+shift_init_Mi[id]-1-overlap_low],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low)*sizeof(cuDoubleComplex)), cudaMemcpyHostToDevice, stream5);
        }

        cublasSetStream(cubl_handle_ls, stream6);
        // h(i,j)=dot(w(:,j),v(:,i));
        cubl_status_ls = cublasZdotc(cubl_handle_ls, (shift_end_Mi[id]-shift_init_Mi[id])+1, &w_ls_device[shift_init_Mi[id]-1], 1, &v_ls_device[shift_init_Mi[id]-1], 1, &h_device[(j*(restart+1))+i]);

        //cudaStat1=cudaDeviceSynchronize();
        //MPI_Barrier(upt_comm); 

        cudaStat1 = cudaMemcpyAsync(dot_temp, &h_device[(j*(restart+1))+i],(size_t)(sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost, stream6);

        //cudaStat1=cudaDeviceSynchronize();
        //MPI_Barrier(upt_comm);

        MPI_Allreduce(dot_temp, dot, 1 , MPI_DOUBLE_COMPLEX, MPI_SUM, upt_comm);

        //MPI_Barrier(upt_comm);

        cudaStat1 = cudaMemcpyAsync(&h_device[(j*(restart+1))+i], dot,(size_t)(sizeof(cuDoubleComplex)), cudaMemcpyHostToDevice, stream6);

        //MPI_Barrier(upt_comm);


        //  w(:,j)=w(:,j)-h(i,j)*v(:,i);
        if(id == 0)
        {
        vct1_sub_mul_vct_kernel<<<numBlock, threadPerBlock, 0, stream6>>>( &w_ls_device[shift_init_Mi[id]-1], &v_ls_device[shift_init_Mi[id]-1], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high, &h_device[(j*(restart+1))+i]);
        }
        else if(id != 0 && id != num_procs-1)
        { 
        vct1_sub_mul_vct_kernel<<<numBlock, threadPerBlock, 0, stream6>>>( &w_ls_device[shift_init_Mi[id]-1-overlap_low], &v_ls_device[shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol, &h_device[(j*(restart+1))+i]);
        }
        else
        {
        vct1_sub_mul_vct_kernel<<<numBlock, threadPerBlock, 0, stream6>>>( &w_ls_device[shift_init_Mi[id]-1-overlap_low], &v_ls_device[shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low, &h_device[(j*(restart+1))+i]);
        }
      }

      cudaDeviceSynchronize();
      cublasSetStream(cubl_handle_ls, stream7);

      // h(j+1,j)=norm(w(:,j));
      //cubl_status_ls = cublasDznrm2(cubl_handle_ls, n_ham, &w_ls_device[j*n_ham], 1, &h_device[(j*(restart+1))+j+1].x);
      cubl_status_ls = cublasZdotc(cubl_handle_ls, (shift_end_Mi[id]-shift_init_Mi[id])+1, &w_ls_device[shift_init_Mi[id]-1], 1, &w_ls_device[shift_init_Mi[id]-1], 1, dot_temp_device_ls);

      //MPI_Barrier(upt_comm);

      cudaStat1 = cudaMemcpyAsync(dot_temp_1, dot_temp_device_ls,(size_t)(sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost, stream7);

      //MPI_Barrier(upt_comm);

      MPI_Allreduce(dot_temp_1, dot_temp_2, 1 , MPI_DOUBLE_COMPLEX, MPI_SUM, upt_comm);

      //MPI_Barrier(upt_comm);

      *norm_temp = sqrt(dot_temp_2[0].x);

      cudaStat1 = cudaMemcpyAsync(&h_device[(j*(restart+1))+j+1].x, norm_temp,(size_t)(sizeof(double)), cudaMemcpyHostToDevice, stream7); 

      cudaStat1 = cudaMemcpyAsync(test, &h_device[(j*(restart+1))+j+1],(size_t)(sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost, stream7);

      //if h(j+1,j)==0
      //if(h_device[((j+1)*(restart+1))+j].x == 0.0000)
      if(test[0].x == 0.0000)
        {
         restart=j;
        }
      else
        {
           //v(:,j+1)=w(:,j)/h(j+1,j);
           if(id == 0)
           {
           vct_1_div_asg_to_vct_2_kernel<<<numBlock, threadPerBlock, 0, stream7>>>( &w_ls_device[shift_init_Mi[id]-1], &v_ls_device[shift_init_Mi[id]-1], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high, &h_device[(j*(restart+1))+j+1].x);
           }
           else if(id != 0 && id != num_procs-1)
           { 
           vct_1_div_asg_to_vct_2_kernel<<<numBlock, threadPerBlock, 0, stream7>>>( &w_ls_device[shift_init_Mi[id]-1-overlap_low], &v_ls_device[shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol, &h_device[(j*(restart+1))+j+1].x);
           }
           else
          {
          vct_1_div_asg_to_vct_2_kernel<<<numBlock, threadPerBlock, 0, stream7>>>( &w_ls_device[shift_init_Mi[id]-1-overlap_low], &v_ls_device[shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low, &h_device[(j*(restart+1))+j+1].x);
          }
          cudaDeviceSynchronize();
          //cudaStat1 = cudaMemcpy(&v_ls[(j+1)*n_ham], v_ls_device,(size_t)(n_ham*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost);
          if(id ==0)
          cudaStat1 = cudaMemcpyAsync(&v_ls[(j+1)*n_ham+(shift_init_Mi[id]-1)], &v_ls_device[shift_init_Mi[id]-1], (size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high)*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost, stream6);
          else if(id != 0 && id != num_procs-1)
          cudaStat1 = cudaMemcpyAsync(&v_ls[(j+1)*n_ham+(shift_init_Mi[id]-1-overlap_low)], &v_ls_device[shift_init_Mi[id]-1-overlap_low], (size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol)*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost, stream6);
          else
          cudaStat1 = cudaMemcpyAsync(&v_ls[(j+1)*n_ham+(shift_init_Mi[id]-1-overlap_low)], &v_ls_device[shift_init_Mi[id]-1-overlap_low], (size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low)*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToHost, stream6);
        }

   }

        cublasSetStream(cubl_handle_ls, NULL);


        // g(1:m+1,:)=0;
        cudaStat1 = cudaMemset(g_device, 0, sizeof(cuDoubleComplex)*(restart+1));
        

        // g(1,:)=beta;
        cudaStat2 = cudaMemcpy(beta, beta_device,(size_t)(sizeof(double)), cudaMemcpyDeviceToHost);
        test[0].x = *beta; test[0].y = 0.0000;
        cudaStat1 = cudaMemcpy(g_device, test,(size_t)(sizeof(cuDoubleComplex)), cudaMemcpyHostToDevice);
        //if(cudaStat1 != cudaSucess

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

     // sqrt(h(j+1,j)^2 + h(j,j)^2);
     // h(j+1,j)^2
     cudoublecomplex_temp1.x = (h[(j*(restart+1))+j+1].x * h[(j*(restart+1))+j+1].x) - (h[(j*(restart+1))+j+1].y * h[(j*(restart+1))+j+1].y);
     cudoublecomplex_temp1.y = (h[(j*(restart+1))+j+1].x * h[(j*(restart+1))+j+1].y) + (h[(j*(restart+1))+j+1].y * h[(j*(restart+1))+j+1].x);
     //h(j,j)^2
     cudoublecomplex_temp2.x = (h[j*(restart+1)+j].x * h[j*(restart+1)+j].x) - (h[j*(restart+1)+j].y * h[j*(restart+1)+j].y);
     cudoublecomplex_temp2.y = (h[j*(restart+1)+j].x * h[j*(restart+1)+j].y) + (h[j*(restart+1)+j].y * h[j*(restart+1)+j].x);

     cudoublecomplex_temp.x = cudoublecomplex_temp1.x + cudoublecomplex_temp2.x;
     cudoublecomplex_temp.y = cudoublecomplex_temp1.y + cudoublecomplex_temp2.y;

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

     // sino=h(j+1,j)/(sqrt(h(j+1,j)^2 + h(j,j)^2));
     sino.x = ((h[(j*(restart+1))+j+1].x * sqrt_cudoublecomplex_temp.x) + (h[(j*(restart+1))+j+1].y * sqrt_cudoublecomplex_temp.y)) / ((sqrt_cudoublecomplex_temp.x * sqrt_cudoublecomplex_temp.x) + (sqrt_cudoublecomplex_temp.y * sqrt_cudoublecomplex_temp.y));

     sino.y = ((h[(j*(restart+1))+j+1].y * sqrt_cudoublecomplex_temp.x) - (h[(j*(restart+1))+j+1].x * sqrt_cudoublecomplex_temp.y)) / ((sqrt_cudoublecomplex_temp.x * sqrt_cudoublecomplex_temp.x) + (sqrt_cudoublecomplex_temp.y * sqrt_cudoublecomplex_temp.y));


     // coso=h(j,j)/(sqrt(h(j+1,j)^2 + h(j,j)^2));
     coso.x = ((h[j*(restart+1)+j].x * sqrt_cudoublecomplex_temp.x) + (h[j*(restart+1)+j].y * sqrt_cudoublecomplex_temp.y)) / ((sqrt_cudoublecomplex_temp.x * sqrt_cudoublecomplex_temp.x) + (sqrt_cudoublecomplex_temp.y * sqrt_cudoublecomplex_temp.y));

     coso.y = ((h[j*(restart+1)+j].y * sqrt_cudoublecomplex_temp.x) - (h[j*(restart+1)+j].x * sqrt_cudoublecomplex_temp.y)) / ((sqrt_cudoublecomplex_temp.x * sqrt_cudoublecomplex_temp.x) + (sqrt_cudoublecomplex_temp.y * sqrt_cudoublecomplex_temp.y));

     // P(j,j)=conj(coso);
     P[j*(restart+1)+j] = cuConj(coso);

     // P(j+1,j+1)=coso;
     P[(j+1)*(restart+1)+j+1].x = coso.x;  P[(j+1)*(restart+1)+j+1].y = coso.y;

     // P(j,j+1)=conj(sino);
     P[(j+1)*(restart+1)+j] = cuConj(sino);

     // P(j+1,j)=-sino;
     P[(j*(restart+1))+j+1].x = -sino.x;  P[(j*(restart+1))+j+1].y = -sino.y;


     cudaStat2 = cudaMemcpy(P_device, P,(size_t)(sizeof(cuDoubleComplex)*(restart+1)*(restart+1)), cudaMemcpyHostToDevice);
     cudaStat1 = cudaMemcpy(h_device, h,(size_t)(sizeof(cuDoubleComplex)*restart*(restart+1)), cudaMemcpyHostToDevice);
     
     cudaDeviceSynchronize();

     cublasSetStream(cubl_handle_ls, stream5);
     // h=P*h;
     cubl_status_ls = cublasZgemm(cubl_handle_ls, CUBLAS_OP_N, CUBLAS_OP_N, restart+1, restart, restart+1, scalar1_device_ls, P_device, restart+1, h_device, restart+1, scalar2_device_ls, h_device, restart+1);

    cublasSetStream(cubl_handle_ls, stream6);
    // g=P*g;
    cubl_status_ls = cublasZgemv(cubl_handle_ls, CUBLAS_OP_N, restart+1, restart+1, scalar1_device_ls, P_device, restart+1, g_device, 1, scalar2_device_ls, g_device, 1);

    // update h since h = P*h
    cudaStat1 = cudaMemcpyAsync(h, h_device,(size_t)(sizeof(cuDoubleComplex)*restart*(restart+1)), cudaMemcpyDeviceToHost, stream5);

    free(P);
    cudaFree(P_device);
   }

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

   // LAPACKE_zgelsy( int matrix_layout, lapack_int m, lapack_int n, lapack_int nrhs, lapack_complex_double* a, lapack_int lda, lapack_complex_double* b, lapack_int ldb, lapack_int* jpvt, double rcond, lapack_int* rank );

   error_code = LAPACKE_zgelsy(LAPACK_COL_MAJOR, restart, restart, 1, (MKL_Complex16 *)h_min, restart, (MKL_Complex16 *)g_min, restart, jpvt, rcond, eff_rank_h_min);


   cublasSetStream(cubl_handle_ls, NULL);

   // xm=x0+v(:,1:m)*minimizer;
   cudaStat1 = cudaMemcpy(minimizer_device, g_min,(size_t)(sizeof(cuDoubleComplex)*restart), cudaMemcpyHostToDevice);
   //cudaStat1 = cudaMemcpy(v_ls_device, v_ls,(size_t)(n_ham*sizeof(cuDoubleComplex)), cudaMemcpyHostToDevice);

   cubl_status_ls = cublasZgemv(cubl_handle_ls, CUBLAS_OP_N, n_ham, restart, scalar1_device_ls, v_ls, n_ham, minimizer_device, 1, scalar2_device_ls, temp_device, 1);
   //vct1_add_vct2_asg_vct3_kernel<<<numBlock, threadPerBlock>>>(X0_device, temp_device, temp2_device, n_ham);
   if(id == 0)
   {
   vct1_add_vct2_asg_vct3_kernel<<<numBlock, threadPerBlock>>>(&X0_device[shift_init_Mi[id]-1], &temp_device[shift_init_Mi[id]-1], &temp2_device[shift_init_Mi[id]-1], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high);
   }
   else if(id != 0 && id != num_procs-1)
   { 
   vct1_add_vct2_asg_vct3_kernel<<<numBlock, threadPerBlock>>>(&X0_device[shift_init_Mi[id]-1-overlap_low], &temp_device[shift_init_Mi[id]-1-overlap_low], &temp2_device[shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol);
   }
   else
   {
   vct1_add_vct2_asg_vct3_kernel<<<numBlock, threadPerBlock>>>(&X0_device[shift_init_Mi[id]-1-overlap_low], &temp_device[shift_init_Mi[id]-1-overlap_low], &temp2_device[shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low);
   }

   // if abs(g(m+1,1))<tol 
   if(sqrt(g[restart].x*g[restart].x+g[restart].y*g[restart].y)<ls_tol)
   {   
     // x = xm;
     cudaStat1 = cudaMemcpy(t_device, temp2_device,(size_t)(sizeof(cuDoubleComplex)*n_ham), cudaMemcpyDeviceToDevice);

     *ls_counter =  *ls_counter + restart_count * restart;

     //printf("i m out from here, gmres converged\n");

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
     cudaFree(dot_temp_device_ls);
     cudaStreamDestroy(stream5);
     cudaStreamDestroy(stream6);
     cudaStreamDestroy(stream7);
     cudaStreamDestroy(stream8);
     cudaFreeHost(v_ls);
     cudaFree(buffer_device);

     free(norm_r_ls);
     free(beta);
     free(h);
     free(h_min);
     free(g);
     free(g_min);
     free(test);
     free(jpvt);
     free(eff_rank_h_min);
     free(reqs_ls);
     free(dot_temp);
     free(dot_temp_1);
     free(dot_temp_2);
     free(norm_temp);

     cudaDeviceSynchronize();

     return;
   }
   else
   {
     // x0=xm;  
     //cudaStat1 = cudaMemcpy(X0_device, temp2_device,(size_t)(sizeof(cuDoubleComplex)*n_ham), cudaMemcpyDeviceToDevice);
     if(id ==0)
     cudaStat1 = cudaMemcpy(&X0_device[shift_init_Mi[id]-1], &temp2_device[shift_init_Mi[id]-1],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high)*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToDevice);
     else if(id != 0 && id != num_procs-1)
     cudaStat1 = cudaMemcpy(&X0_device[shift_init_Mi[id]-1-overlap_low], &temp2_device[shift_init_Mi[id]-1-overlap_low],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol)*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToDevice);
     else
     cudaStat1 = cudaMemcpy(&X0_device[shift_init_Mi[id]-1-overlap_low], &temp2_device[shift_init_Mi[id]-1-overlap_low],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low)*sizeof(cuDoubleComplex)), cudaMemcpyDeviceToDevice);
     // restart=restart+1;
     restart_count = restart_count + 1;
     // if abs(g(m+1,1)) <= normmin
     if(sqrt(g[restart].x*g[restart].x+g[restart].y*g[restart].y)<= norm_min)
     {
       // xmin = xm;
       cudaStat1 = cudaMemcpyAsync(Xmin_device, temp2_device,(size_t)(sizeof(cuDoubleComplex)*n_ham), cudaMemcpyDeviceToDevice, stream6);
       norm_min = sqrt(g[restart].x*g[restart].x+g[restart].y*g[restart].y);
     }
   }

  } // end of while

     cudaStat1=cudaDeviceSynchronize();


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
 cudaFree(dot_temp_device_ls);
 cudaStreamDestroy(stream5);
 cudaStreamDestroy(stream6);
 cudaStreamDestroy(stream7);
 cudaStreamDestroy(stream8);
 cudaFreeHost(v_ls);
 cudaFree(buffer_device);

 free(norm_r_ls);
 free(beta);
 free(h);
 free(h_min);
 free(g);
 free(g_min);
 free(test);
 free(jpvt);
 free(eff_rank_h_min);
 free(dot_temp);
 free(dot_temp_1);
 free(dot_temp_2);
 free(norm_temp);
 free(reqs_ls);

 cudaDeviceSynchronize();

 return;
}



























__global__ void vct1_add_vct2_asg_vct3_kernel(cuDoubleComplex *vct1, cuDoubleComplex * vct2, cuDoubleComplex * vct3, int n_ham)
{
  int ID;
  ID = blockDim.x * blockIdx.x + threadIdx.x;
  if(ID < n_ham)
  {
   vct3[ID].x = vct1[ID].x + vct2[ID].x;
   vct3[ID].y = vct1[ID].y + vct2[ID].y;
  }
  __syncthreads();
}




__global__ void vct1_sub_vct2_asg_vct3_kernel(cuDoubleComplex *vct1, cuDoubleComplex * vct2, cuDoubleComplex * vct3, int n_ham)
{
  int ID;
  ID = blockDim.x * blockIdx.x + threadIdx.x;
  if(ID < n_ham)
  {
   vct3[ID].x = vct1[ID].x - vct2[ID].x;
   vct3[ID].y = vct1[ID].y - vct2[ID].y;
  }
  __syncthreads();
}



__global__ void vct1_neg_asg_vct2_kernel(cuDoubleComplex *vct1, cuDoubleComplex * vct2, int n_ham)
{
  int ID;
  ID = blockDim.x * blockIdx.x + threadIdx.x;
  if(ID < n_ham)
  {
    vct1[ID].x =  -vct2[ID].x;
    vct1[ID].y =  -vct2[ID].y;
   }
   __syncthreads();
}



__global__ void vct_div_slr_minus_scl_vct_kernel(cuDoubleComplex * scr1, cuDoubleComplex * scr2, cuDoubleComplex * des, int n_ham, double Slr1, double Slr2)
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



__global__ void vct_div_slr_kernel(cuDoubleComplex * scr, cuDoubleComplex * des, int n_ham, double Slr)
{
 int ID;
 ID = blockDim.x * blockIdx.x + threadIdx.x;
 if(ID < n_ham)
 {
   des[ID].x = scr[ID].x / Slr;
   des[ID].y = scr[ID].y / Slr;
 }
 __syncthreads();
}


__global__ void vct_1_div_asg_to_vct_2_kernel(cuDoubleComplex * vct_1, cuDoubleComplex * vct_2, int n_ham, double *beta)
{
  int ID;
  ID = blockDim.x * blockIdx.x + threadIdx.x;
  if(ID < n_ham)
  {
    vct_2[ID].x =  vct_1[ID].x / beta[0];
    vct_2[ID].y =  vct_1[ID].y / beta[0];
  }
  __syncthreads();
}








__global__ void vct1_sub_mul_vct_kernel(cuDoubleComplex * vct1, cuDoubleComplex * vct2, int n_ham, cuDoubleComplex * dot)
{
  int  ID;
  ID = blockDim.x * blockIdx.x + threadIdx.x;
  if(ID < n_ham)
  {
    vct1[ID].x = vct1[ID].x - ((dot[0].x * vct2[ID].x) - (dot[0].y * vct2[ID].y));
    vct1[ID].y = vct1[ID].y - ((dot[0].x * vct2[ID].y) + (dot[0].y * vct2[ID].x));
  }
  __syncthreads();	
}







__global__ void vct_sub_scl_mul_vct_kernel(cuDoubleComplex *vct1, cuDoubleComplex * vct2, double scalar, int n_ham)
{
  int ID;
  ID = blockDim.x * blockIdx.x + threadIdx.x;
  if(ID < n_ham)
  {
    vct1[ID].x =  vct1[ID].x - (scalar * vct2[ID].x);
    vct1[ID].y =  vct1[ID].y - (scalar * vct2[ID].y);
   }
   __syncthreads();
}





__global__ void shift_A_kernel(float * val, int * row, int * col, int n_ham, float shift, int offset)
{
  int ID;
  ID = blockDim.x * blockIdx.x + threadIdx.x;
  if(ID < n_ham)
  {
      for(int j=row[ID+offset]-1; j <= row[ID+1+offset]-1; j++)
      {
       if(col[j] == ID+1+offset)
       {
        val[j] = val[j] + shift; 
        break;
       }
      }
  }
  __syncthreads();
}




__global__ void cpy_vct_1_to_vct_2_kernel(cuDoubleComplex * vct_scr, cuDoubleComplex * vct_des, int n_ham)
{
  int ID;
  ID = blockDim.x * blockIdx.x + threadIdx.x;
  if(ID < n_ham)
  {
    vct_des[ID].x = vct_scr[ID].x;
    vct_des[ID].y = vct_scr[ID].y;
  }
  __syncthreads();
}


__global__ void mv_kernel(int n_ham, int col, cuDoubleComplex * xVal_kr, cuDoubleComplex * y_kr, cuDoubleComplex * Finalans_kr){
			//kernel func varaibles
			int ID, row_start, row_end, jj;
			double dot_x, dot_img;
			ID = blockDim.x*blockIdx.x+threadIdx.x;
				if(ID < n_ham){
					dot_x=0.0;
					dot_img=0.0; 
					for(jj = 0; jj < col; jj++ ){
						  dot_x += ((xVal_kr[(jj*n_ham)+ID].x * y_kr[jj].x) - (xVal_kr[(jj*n_ham)+ID].y * y_kr[jj].y)); 
						  dot_img  += ((xVal_kr[(jj*n_ham)+ID].x * y_kr[jj].y) + (xVal_kr[(jj*n_ham)+ID].y * y_kr[jj].x)); 
						}
		
					Finalans_kr[ID].x = dot_x;  
					Finalans_kr[ID].y = dot_img;  
					
				}
                                __syncthreads();
		}





void setdevicebeforeinit_() {
  char * localRankStr = NULL;
  int rank = 0, devCount = 0;
  cudaError_t cudaStat1;
  // We extract the local rank initialization using an environment variable
  if ((localRankStr = getenv(ENV_LOCAL_RANK)) != NULL)
  {
     rank = atoi(localRankStr);	
  }
  cudaDeviceReset();
  cudaThreadExit();
  cudaGetDeviceCount(&devCount);
  
  printf("device count = %d\n", devCount);

  if (devCount > 0)
  {	  
     cudaStat1 = set_gpu_id(rank % devCount);
     if(cudaStat1 != cudaSuccess)
     printf("ERROR DEVICE SET FAILED\n");
  }
  else
  {
     printf("NO CUDA DEVICE FOUND: PLEASE CHANGE SOLVER\n");
  }
}




__global__ void spmv_csr_hybrid_kernel(int dimRow, const int* rowPtrs, const int* colIdxs, const float* values, 
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
    __syncthreads();

}



__global__ void spmv_csr_hybrid_kernel(int num_rows, const int* rowPtrs, const int* colIdxs, const float* values, 
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
   __syncthreads();
}


__global__ void vct_pls_scl_mul_vct_kernel_jd(cuDoubleComplex *vct1, cuDoubleComplex * vct2, cuDoubleComplex scalar, int n_ham)
{
  int ID;
  ID = blockDim.x * blockIdx.x + threadIdx.x;
  if(ID < n_ham)
  {
    vct1[ID].x =  vct1[ID].x + scalar.x * vct2[ID].x - scalar.y * vct2[ID].y;
    vct1[ID].y =  vct1[ID].y + scalar.x * vct2[ID].y + scalar.y * vct2[ID].x;
  }
  __syncthreads();
}

