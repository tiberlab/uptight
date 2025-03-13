    // !=============================================================================
    // !
    // !                              JD SOLVER
    // !
    // !=============================================================================
    // !
    // ! Walter Rodrigues
    // ! Dipartimento di Ingegneria Elettronica
    // ! Universita` di Roma "Tor Vergata"
    // ! 17-04-2014
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
    
    #include <mkl.h>
    #include <mkl_lapacke.h>
    
    
    #include <unistd.h>
    #include <sys/unistd.h>
    #include <ctype.h>
    
    #include <time.h>
    #include <sys/time.h>
    
    
   
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
    
    

    
     void GMRES_CPU(float *valptr_real, int *rowptr_real, int *colptr_real, float *valptr_img, int *rowptr_img, int *colptr_img, int n_ham, int size_mat_real, double complex * r, double ls_tol, int restart, int maxit, double complex *Q_bar, double complex *t, int k, int *ls_counter, int * shift_init_Mi, int * shift_end_Mi, int overlap_high, int overlap_low, int num_procs, int id,  MPI_Comm upt_comm);

    void c_zdotu(int size, double complex * array_a, double complex * array_b, double complex * output);
    
    
    void vct_div_slr(double complex * scr, double complex * des, int n_ham, double Slr);
    void vct1_sub_mul_vct(double complex * vct1, double complex * vct2, int n_ham, double complex * dot);
    void vct_sub_scl_mul_vct(double complex *vct1, double complex * vct2, double scalar, int n_ham);
    void shift_A(float * val, int * row, int * col, int n_ham, float shift, int offset);
    void cpy_vct_1_to_vct_2(double complex * vct_scr, double complex * vct_des, int n_ham);
    void vct_1_div_asg_to_vct_2(double complex * vct_1, double complex * vct_2, int n_ham, double beta);
    void vct_div_slr_minus_scl_vct(double complex * scr1, double complex * scr2, double complex * des, int n_ham, double Slr1, double Slr2);
    void vct1_add_vct2_asg_vct3(double complex *vct1, double complex * vct2, double complex * vct3, int n_ham);
    void vct1_neg_asg_vct2(double complex *vct1, double complex * vct2, int n_ham);
    void vct1_sub_vct2_asg_vct3(double complex *vct1, double complex * vct2, double complex * vct3, int n_ham);
    void mv(int n_ham, int col, double complex * xVal_kr, double complex * y_kr, double complex * Finalans_kr);
    void vct_pls_scl_mul_vct(double complex *vct1, double complex * vct2, double complex scalar, int n_ham);
    void spmv_csr_hybrid(int num_rows, const int* rowPtrs, const int* colIdxs, const float* values, const double complex* x, double complex* y, const int offset);


#ifdef __cplusplus
extern "C" {
#endif

    //void *__gxx_personality_v0;

    extern void jd_cpu_no_pc_split_mxprec_pal_(int * N_ham, int * Size_Mat_real, int * Size_Mat_img, float * valptr_real, int * rowptr_real, int * colptr_real, float * valptr_img, int * rowptr_img, int * colptr_img, char * sparse_fmt, int * Band_type, double * JD_tol, double * Shift,  int * JD_Min_step, int * JD_Max_step, int * Num_ev, double * lambda_out, double complex * eigen_vec_out, double * LS_tol, int * LS_restart, int *LS_maxit, int *col_ind_low, int *col_ind_high, int *shift_init, int *shift_end, int * NUM_PROCS, int * ID,  int * UPT_COMM);
 
    //void setdevicebeforeinit_();

#ifdef __cplusplus
}
#endif
    
    
    
    void jd_cpu_no_pc_split_mxprec_pal_(int * N_ham, int * Size_Mat_real, int * Size_Mat_img, float * valptr_real, int * rowptr_real, int * colptr_real, float * valptr_img, int * rowptr_img, int * colptr_img, char * sparse_fmt, int * Band_type, double * JD_tol, double * Shift,  int *JD_Min_step, int * JD_Max_step, int * Num_ev, double * lambda_out, double complex * eigen_vec_out, double * LS_tol, int * LS_restart, int *LS_maxit, int *col_ind_low, int *col_ind_high, int *shift_init, int *shift_end, int * NUM_PROCS, int * ID,  int * UPT_COMM)
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
    double complex * V0;
    double complex * X_bar;
    double complex * t; 
    double complex * vct; 
    double * norm_V0;
    double * norm_vct;
    double * norm_t;
    //double * norm_t_in;
    double tol_shift;
    double complex  * scalar_1, * scalar_2;
    double complex * y;
    double complex * temp_ev;
    double complex * eigen_vec;
    double complex * u;
    double * eigen_val;
    double complex * M;
    double * s_u;
    double lambda[num_ev];
    double * norm_r;
    double complex * Q_bar;;
    double complex * temp1, * temp2;
    double complex * w_bar;
    double complex * r;
    double * s_v;
    double complex * ubar;
    const double complex jcmpx= 0.000000 + 1.0000000 * I;
    double complex * mxv_temp;
    double complex * dot_temp, * dot_temp_1, * dot_temp_2;
    //double * norm_temp;
    double complex * temp_host_v, * temp_host_w;
    double complex * v, * w;
    int f_loop;
    double complex * eigen_vec_temp;
    int i, j;
    int rank;

    
    v = (double complex *) malloc(n_ham*jd_max_step*sizeof(double complex));
    memset(v, 0, n_ham*jd_max_step*sizeof(double complex));
    w = (double complex *) malloc(n_ham*jd_max_step*sizeof(double complex));
    memset(w, 0, n_ham*jd_max_step*sizeof(double complex));
//     temp_host_v = (double complex *) malloc(n_ham*jd_max_step*sizeof(double complex));
//     memset(temp_host_v, 0, n_ham*jd_max_step*sizeof(double complex));
//     temp_host_w = (double complex *) malloc(n_ham*jd_max_step*sizeof(double complex));
//     memset(temp_host_w, 0, n_ham*jd_max_step*sizeof(double complex));
    
    
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
    int * shift_init_M;
    int * shift_end_M;
    
    shift_init_Mi=(int *) malloc(num_procs*sizeof(int));
    shift_end_Mi=(int *) malloc(num_procs*sizeof(int));
    shift_init_M=(int *) malloc(num_procs*sizeof(int));
    shift_end_M=(int *) malloc(num_procs*sizeof(int));

    shift_init_M[id]= * shift_init;
    shift_end_M[id]= * shift_end;
    
//     if ( id == 0 && id == num_procs-1 )
//     { 
//             shift_init_Mi[id] = 1; 
//             shift_end_Mi[id]  = n_ham;
//     }
//     else if ( id == 0)
//     {
//             shift_init_Mi[id] = 1;
//             shift_end_Mi[id] =  row_offset[id+1];
//     }
//     else if (  id != num_procs-1 )
//     {
//             shift_init_Mi[id] = row_offset[id]+1;
//             shift_end_Mi[id] =  row_offset[id+1];
//     }
//     else if (  id == num_procs-1 )
//     {
//             shift_init_Mi[id] = row_offset[id]+1;
//             shift_end_Mi[id] = n_ham;
//     }
    
    
    MPI_Barrier(upt_comm);
    MPI_Gather(&shift_init_M[id], 1, MPI_INT, &shift_init_Mi[id], 1, MPI_INT, 0, upt_comm);
    MPI_Gather(&shift_end_M[id], 1, MPI_INT, &shift_end_Mi[id], 1, MPI_INT, 0, upt_comm);
    MPI_Barrier(upt_comm);
    
    
    MPI_Bcast(shift_init_Mi,num_procs, MPI_INT, 0, upt_comm);
    MPI_Bcast(shift_end_Mi,num_procs, MPI_INT, 0, upt_comm);
    
    
    int overlap_high = 0;
    int overlap_low = 0;
    
    
    if(id == 0){
        
        for (f_loop = 0 ; f_loop <= num_procs-2; f_loop++)
        {
            if(col_ind_low[f_loop+1] < shift_init_Mi[f_loop])
            {
            printf("use less number of nodes in MPI\n");
            exit(0);
            }
        }
    }
    
    MPI_Barrier(upt_comm);
    
    
    
    for(f_loop = 0; f_loop <=num_procs-1; f_loop++)
    {
            if (overlap_high < col_ind_high[f_loop]-shift_end_Mi[f_loop]) 
                                          overlap_high = col_ind_high[f_loop]-shift_end_Mi[f_loop];
    }
    
    
    for(f_loop = 0; f_loop <=num_procs-1; f_loop++)
    {
            if (overlap_low < shift_init_Mi[f_loop]-col_ind_low[f_loop]) 
                                          overlap_low = shift_init_Mi[f_loop]-col_ind_low[f_loop]; 
    }
    
    int overlap_tol = overlap_high + overlap_low;
    
    //printf("overlap_high = %d, overlap_low = %d\n", overlap_high, overlap_low);
    
    
    norm_r = (double *) malloc(sizeof(double));
    norm_t = (double *) malloc(sizeof(double));
    u = (double complex *) malloc(sizeof(double complex)*n_ham);
    memset(u, 0, n_ham*sizeof(double complex));
    norm_vct = (double *) malloc(sizeof(double));
    //norm_t_in = (double *) malloc(sizeof(double));
    //norm_r = (double *) malloc(sizeof(double));
    M_out = (lapack_int *) malloc(sizeof(lapack_int));
    isuppz = (lapack_int *) malloc(2*jd_max_step*sizeof(lapack_int));
    s_u = (double *) malloc(sizeof(double));
    s_v = (double *) malloc(sizeof(double));
    norm_V0 = (double *) malloc(sizeof(double));
    //norm_temp = (double *) malloc(sizeof(double));
    y = (double complex *) malloc(sizeof(double complex));
    dot_temp = (double complex *) malloc(sizeof(double complex));
    dot_temp_1 = (double complex *) malloc(sizeof(double complex));
    dot_temp_2 = (double complex *) malloc(sizeof(double complex));
    scalar_1 = (double complex *) malloc(sizeof(double complex));
    scalar_2 = (double complex *) malloc(sizeof(double complex));

    lapack_int error_code=0;


    int counter;
    
    *scalar_1= 1.00 + 0.00 * I;
    //cimag(scalar_1)= 0.00;
    *scalar_2= 0.00 + 0.00 * I;
    //cimag(scalar_2)= 0.00;
    
    
    int flag =1;
    V0 = (double complex *) malloc(sizeof(double complex)*n_ham);
    
    if(id == 0) {

    for(counter=0; counter < n_ham; counter++)
    {
            V0[counter]= (counter+2.0000000)+ 0.0000000*I; //(1000*rand()) + 0.00I;
            V0[counter]= creal(V0[counter]) + (counter*3.000000)*I; //(1000*rand()) * I;
    }
    
    }
    
    MPI_Barrier(upt_comm);
    
    MPI_Bcast(V0, n_ham, MPI_DOUBLE_COMPLEX, 0, upt_comm);
    
    X_bar = (double complex *) malloc(sizeof(double complex)*n_ham*num_ev);
    memset(X_bar, 0, n_ham*num_ev*sizeof(double complex));
    
    double memory_device = (sizeof(float)*(size_mat_real+size_mat_img)+sizeof(int)*(size_mat_real+n_ham+size_mat_img+n_ham) + sizeof(double complex)*(n_ham*(10+10)) + sizeof(double complex)*(n_ham*(ls_restart+1)) + sizeof(double complex)*(((ls_restart+1)*ls_restart) + ls_restart+1 + ls_restart)) /1000000;

//(sizeof(float)*(size_mat_real+size_mat_img)+sizeof(int)*(size_mat_real+n_ham+size_mat_img+n_ham)+sizeof(double)*(jd_max_step)+sizeof(double complex)*(n_ham*12+n_ham*num_ev+n_ham*jd_max_step*3)+(n_ham*((ls_restart+1)+ls_restart))+(jd_max_step*jd_max_step*2)+((ls_restart+1)*(ls_restart+1))+((ls_restart+1)*ls_restart))/1000000;

    //printf("memory needed per device = %f MB\n", memory_device);

    t = (double complex *) malloc(n_ham*sizeof(double complex));
    memset(t, 0, n_ham*sizeof(double complex));  
    r = (double complex *) malloc(n_ham*sizeof(double complex));
    memset(r, 0, n_ham*sizeof(double complex));
   
    memcpy(t, V0, n_ham*sizeof(double complex));
   
    //printf(" i m here 1\n");


    double teta = 0.000000;
    
    // A = (A - shift * I)
    tol_shift = teta-shift;
    shift_A(valptr_real, rowptr_real, colptr_real, (shift_end_Mi[id]-shift_init_Mi[id])+1, tol_shift, shift_init_Mi[id]-1);

    MPI_Barrier(upt_comm);

    // V0 = V0/norm(V0); t_device = V0
    *norm_V0 = cblas_dznrm2(n_ham, t, 1);
    
    if(id == 0)
    {
      vct_div_slr( &t[shift_init_Mi[id]-1], &t[shift_init_Mi[id]-1], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high, *norm_V0);
    }
    else if(id != 0 && id != num_procs-1)
    { 
      vct_div_slr( &t[shift_init_Mi[id]-1-overlap_low], &t[shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol, *norm_V0);
    }
    else
    {
      vct_div_slr( &t[shift_init_Mi[id]-1-overlap_low], &t[shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low, *norm_V0);
    }
    
    // t=v0; k = 0; m = 0; X_bar = [];
    int k = 0;
    int m = 0;
    
    *s_v = 0;
    int count = 0;
    
    
    while(k < num_ev)
    {
        vct = (double complex *) malloc(n_ham*sizeof(double complex));
                    memset(vct, 0, sizeof(double complex)*n_ham);
        mxv_temp = (double complex *) malloc(n_ham*sizeof(double complex));
                    memset(mxv_temp, 0, sizeof(double complex)*n_ham);
                    
    
        MPI_Barrier(upt_comm);
    
       // vct = A * t;
       spmv_csr_hybrid((shift_end_Mi[id]-shift_init_Mi[id])+1, rowptr_real, colptr_real, valptr_real, t, vct, shift_init_Mi[id]-1);
    
       spmv_csr_hybrid((shift_end_Mi[id]-shift_init_Mi[id])+1, rowptr_img, colptr_img, valptr_img, t, mxv_temp, shift_init_Mi[id]-1);

       MPI_Barrier(upt_comm);
    
       vct_pls_scl_mul_vct(&vct[shift_init_Mi[id]-1], &mxv_temp[shift_init_Mi[id]-1], jcmpx, (shift_end_Mi[id]-shift_init_Mi[id])+1);
       
       MPI_Barrier(upt_comm);

       free(mxv_temp);
    
  
      if (num_procs > 1)
      {
      for(rank =0; rank <=num_procs-2; rank++)
      {
      if(id == rank+1)    MPI_Irecv(&vct[shift_end_Mi[rank]-overlap_low], overlap_low,MPI_DOUBLE_COMPLEX,rank, 111+rank,upt_comm,&reqs[rank+1]);
      if(id == rank)      MPI_Isend(&vct[shift_end_Mi[rank]-overlap_low], overlap_low,MPI_DOUBLE_COMPLEX,rank+1, 111+rank,upt_comm, &reqs[rank+1]);
   
      if (id == rank+1) MPI_Wait(&reqs[rank+1], &status);
      if (id == rank)   MPI_Wait(&reqs[rank+1], &status);
      }
      }
      
      MPI_Barrier(upt_comm);
  
  
      if (num_procs > 1)
      {
      for(rank =0; rank <=num_procs-2; rank++)
      {
      if(id == rank)    MPI_Irecv(&vct[shift_init_Mi[rank+1]-1],overlap_high, MPI_DOUBLE_COMPLEX,rank+1,8+rank,upt_comm, &reqs[rank+1]);
      if(id == rank+1)  MPI_Isend(&vct[shift_init_Mi[rank+1]-1],overlap_high, MPI_DOUBLE_COMPLEX,rank,8+rank,upt_comm, &reqs[rank+1]);
  
      if (id == rank+1) MPI_Wait(&reqs[rank+1], &status);
      if (id == rank)   MPI_Wait(&reqs[rank+1], &status);
      }
      }
  
      MPI_Barrier(upt_comm);
      

      for(i = 0; i < m; i++)
      {
      
      // y = w(:,i)' * vct;
      //c_zdotu((shift_end_Mi[id]-shift_init_Mi[id])+1, &w[i*n_ham+shift_init_Mi[id]-1], &vct[shift_init_Mi[id]-1], dot_temp);
      cblas_zdotc_sub((shift_end_Mi[id]-shift_init_Mi[id])+1, &w[i*n_ham+shift_init_Mi[id]-1], 1, &vct[shift_init_Mi[id]-1], 1, dot_temp);
  
      MPI_Barrier(upt_comm);
  
      MPI_Allreduce(dot_temp, y, 1 , MPI_DOUBLE_COMPLEX, MPI_SUM, upt_comm);
  
      // vct = vct - y * w(:,i);
      if(id == 0)
      {
      vct1_sub_mul_vct(&vct[shift_init_Mi[id]-1], &w[i*n_ham+shift_init_Mi[id]-1], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high, y);
      }
      else if(id != 0 && id != num_procs-1)
      {
      vct1_sub_mul_vct(&vct[shift_init_Mi[id]-1-overlap_low], &w[i*n_ham+shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol, y);
      }
      else
      {
      vct1_sub_mul_vct(&vct[shift_init_Mi[id]-1-overlap_low], &w[i*n_ham+shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low, y);
      }
  
      // t = t - y * v(:,i);
      if(id == 0)
      {
      vct1_sub_mul_vct(&t[shift_init_Mi[id]-1], &v[i*n_ham+shift_init_Mi[id]-1], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high, y);
      }
      else if(id != 0 && id != num_procs-1)
      {
      vct1_sub_mul_vct(&t[shift_init_Mi[id]-1-overlap_low], &v[i*n_ham+shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol, y);
      }
      else
      {
      vct1_sub_mul_vct(&t[shift_init_Mi[id]-1-overlap_low], &v[i*n_ham+shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low, y);
      }
      
      }
    
  
      m = m+1;

         //cblas_dznrm2 (const MKL_INT N, const void *X, const MKL_INT incX);
         //*norm_vct = cblas_dznrm2(n_ham, vct, 1);
      //c_zdotu((shift_end_Mi[id]-shift_init_Mi[id])+1, &vct[shift_init_Mi[id]-1], &vct[shift_init_Mi[id]-1], dot_temp_1);
      cblas_zdotc_sub((shift_end_Mi[id]-shift_init_Mi[id])+1, &vct[shift_init_Mi[id]-1], 1, &vct[shift_init_Mi[id]-1], 1, dot_temp_1);
  
      MPI_Barrier(upt_comm);
  
      MPI_Allreduce(dot_temp_1, dot_temp_2, 1 , MPI_DOUBLE_COMPLEX, MPI_SUM, upt_comm);
  
      *norm_vct = sqrt(creal(*dot_temp_2));

      // w(:,m) = vct/norm(vct);
      // v(:,m) = t/norm(vct);
  
      if(id == 0)
      {
      vct_div_slr(&vct[shift_init_Mi[id]-1], &w[(m-1)*n_ham+shift_init_Mi[id]-1], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high, *norm_vct);
      }
      else if(id != 0 && id != num_procs-1)
      {
      vct_div_slr(&vct[shift_init_Mi[id]-1-overlap_low], &w[(m-1)*n_ham+shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol, *norm_vct);
      }
      else
      {
      vct_div_slr(&vct[shift_init_Mi[id]-1-overlap_low], &w[(m-1)*n_ham+shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low, *norm_vct);
      }
  
      if(id == 0)
      {
      vct_div_slr(&t[shift_init_Mi[id]-1], &v[(m-1)*n_ham+shift_init_Mi[id]-1], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high, *norm_vct);
      }
      else if(id != 0 && id != num_procs-1)
      {
      vct_div_slr(&t[shift_init_Mi[id]-1-overlap_low], &v[(m-1)*n_ham+shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol, *norm_vct);
      }
      else
      {
      vct_div_slr(&t[shift_init_Mi[id]-1-overlap_low], &v[(m-1)*n_ham+shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low, *norm_vct);
      }


      free(vct);
      
      // alloc M(m,m);
      M = (double complex *) malloc(m*m*sizeof(double complex));
      memset(M, 0, m*m*sizeof(double complex));  

      for(j=0; j < m; j++)
      { 
      for(i = 0; i < m; i++)
      {
      // M(i,m) = w(:,i)' * v(:,m);
      //c_zdotu((shift_end_Mi[id]-shift_init_Mi[id])+1, &w[i*n_ham+shift_init_Mi[id]-1], &v[j*n_ham+shift_init_Mi[id]-1], dot_temp);
      cblas_zdotc_sub((shift_end_Mi[id]-shift_init_Mi[id])+1, &w[i*n_ham+shift_init_Mi[id]-1], 1, &v[j*n_ham+shift_init_Mi[id]-1], 1, dot_temp);
  
      MPI_Barrier(upt_comm);
  
      MPI_Allreduce(dot_temp, &M[j*m+i], 1 , MPI_DOUBLE_COMPLEX, MPI_SUM, upt_comm);

      // M(m,i) = w(:,m)' * v(:,i);
      M[m*i+j] = conj(M[j*m+i]);
      }
      }
  
//       free(temp1);
      
      eigen_vec = (double complex *) malloc(sizeof(double complex)*m*m);
      eigen_val = (double *) malloc(sizeof(double)*m);
      eigen_vec_temp = (double complex *) malloc(sizeof(double complex)*m*m);
      memset(eigen_vec, 0, m*m*sizeof(double complex));
      memset(eigen_vec_temp, 0, m*m*sizeof(double complex));
      memset(eigen_val, 0, m*sizeof(double));

      error_code = LAPACKE_zheevr( LAPACK_COL_MAJOR, 'V', 'A', 'L', m, (MKL_Complex16 *)M, m, VL, VU, IL, IU, ABSTOL, M_out, eigen_val, (MKL_Complex16 *)eigen_vec, m, isuppz );

  
      if(m != 1)
      {
      int index[m], temp_index;
      double temp;
      //temp_ev = (double complex *) malloc(sizeof(double complex)*m);
      //memset(temp_ev, 0, m*sizeof(double complex));
  
      for(i =0; i < m; i++)
      {
        index[i] = i;
      }
  
      // sort(eigen_val - shift)  band_type = 1 for CB and band_type = 2 for VB
      if(band_type == 1)
      {
      for(i =0; i < m-1; i++)
      {
        for(j = 0; j < m-1-i; j++)
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
        for(i =0; i < m-1; i++)
      {
        for(j = 0; j < m-1-i; j++)
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
        for(i =0; i < m; i++)
        {
          memcpy(&eigen_vec_temp[i*m], &eigen_vec[index[i]*m],(sizeof(double complex)*m));
        }
      }
      
  
        if(m == 1) 
        { 
          memcpy(eigen_vec_temp, eigen_vec,(sizeof(double complex)*m*m));
        }
        else 
        {
          memcpy(eigen_vec, eigen_vec_temp,(sizeof(double complex)*m*m));
        }

        free(eigen_vec_temp);
       
        ubar = (double complex *) malloc(n_ham*sizeof(double complex));
        memset(ubar, 0, n_ham*sizeof(double complex));
       
        w_bar = (double complex *) malloc(n_ham*sizeof(double complex));
        memset(w_bar, 0, n_ham*sizeof(double complex));

  
        // u_bar = v * eig_mat(:,1);
        cblas_zgemv(CblasColMajor, CblasNoTrans, n_ham, m, scalar_1, v, n_ham, eigen_vec, 1, scalar_2, ubar, 1);
  
        // w_bar = w * eig_mat(:,1);
        cblas_zgemv(CblasColMajor, CblasNoTrans, n_ham, m, scalar_1, w, n_ham, eigen_vec, 1, scalar_2, w_bar, 1); 

        MPI_Barrier(upt_comm);
        
//        for(i = 0; i< n_ham; i++)
  
        // s_u = norm(u_bar);
        //c_zdotu((shift_end_Mi[id]-shift_init_Mi[id])+1, &ubar[shift_init_Mi[id]-1], &ubar[shift_init_Mi[id]-1], dot_temp);
        cblas_zdotc_sub((shift_end_Mi[id]-shift_init_Mi[id])+1, &ubar[shift_init_Mi[id]-1], 1, &ubar[shift_init_Mi[id]-1], 1, dot_temp);

        MPI_Barrier(upt_comm);
  
        MPI_Allreduce(dot_temp, dot_temp_2, 1 , MPI_DOUBLE_COMPLEX, MPI_SUM, upt_comm);
  
        *s_u = sqrt(creal(*dot_temp_2));
  
        //u = (double complex *) malloc(n_ham*sizeof(double complex));
        memset(u, 0, n_ham*sizeof(double complex));
  
        // u = u_bar/s_u; 
        if(id == 0)
        {
          vct_div_slr(&ubar[shift_init_Mi[id]-1], &u[shift_init_Mi[id]-1], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high, *s_u);
        }
        else if(id != 0 && id != num_procs-1)
        {
          vct_div_slr(&ubar[shift_init_Mi[id]-1-overlap_low], &u[shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol, *s_u);
        }
        else
        {
          vct_div_slr(&ubar[shift_init_Mi[id]-1-overlap_low], &u[shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low, *s_u);
        }
  
        free(ubar);
  
        // s_v = eig_val(1)/ (s_u * s_u);
        *s_v = eigen_val[0]/( (*s_u) * (*s_u));
  
        // r = (w_bar / s_u) - (s_v * u);
        if(id == 0)
        {
          vct_div_slr_minus_scl_vct(&w_bar[shift_init_Mi[id]-1], &u[shift_init_Mi[id]-1], &r[shift_init_Mi[id]-1], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high, *s_u, *s_v);
        }
        else if(id != 0 && id != num_procs-1)
        {
          vct_div_slr_minus_scl_vct(&w_bar[shift_init_Mi[id]-1-overlap_low], &u[shift_init_Mi[id]-1-overlap_low], &r[shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol, *s_u, *s_v);
        }
        else
        {
          vct_div_slr_minus_scl_vct(&w_bar[shift_init_Mi[id]-1-overlap_low], &u[shift_init_Mi[id]-1-overlap_low], &r[shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low, *s_u, *s_v);
        }
  
        free(w_bar);
  
        // norm(r_bar);
        //c_zdotu((shift_end_Mi[id]-shift_init_Mi[id])+1, &r[shift_init_Mi[id]-1], &r[shift_init_Mi[id]-1], dot_temp);
        cblas_zdotc_sub((shift_end_Mi[id]-shift_init_Mi[id])+1, &r[shift_init_Mi[id]-1], 1, &r[shift_init_Mi[id]-1], 1, dot_temp);

        MPI_Allreduce(dot_temp, dot_temp_2, 1 , MPI_DOUBLE_COMPLEX, MPI_SUM, upt_comm);

        *norm_r = sqrt(creal(dot_temp_2[0]));
          
        
        while(*norm_r < *JD_tol)
        {
  
         //"inside harmonic part"        
  
          k = k + 1;
  
          // X_bar = [X_bar, u];
   
          // gather u from all IDs
          if (num_procs > 1) 
          {
            int rank;
            for (rank =0; rank <=  num_procs-1; rank++) 
              {
              if(id == 0 && rank != 0)  MPI_Irecv(&u[shift_init_Mi[rank]-1],(shift_end_Mi[rank]-shift_init_Mi[rank])+1, MPI_DOUBLE_COMPLEX, rank, 111+rank, upt_comm, &reqs[rank+1]);
              if(id != 0 && id == rank) MPI_Isend(&u[shift_init_Mi[rank]-1], (shift_end_Mi[rank]-shift_init_Mi[rank])+1, MPI_DOUBLE_COMPLEX, 0, 111+rank, upt_comm, &reqs[rank+1]);
  
              if (id != 0 && id == rank) MPI_Wait(&reqs[rank+1], &status);
              }
  
          }
          
          //printf(" i m here 13\n");
  
          // brodcast u
          MPI_Barrier(upt_comm);
          MPI_Bcast(u, n_ham, MPI_DOUBLE_COMPLEX, 0, upt_comm);
  
          memcpy ( &X_bar[(k-1)*n_ham], u, sizeof(double complex)*n_ham );
  
          // lambda(k) = s_v + shift
          lambda[k-1] = *s_v + shift;
          //if(id == 0) printf("Before sort lambda[%d] = %f \n", k, lambda[k-1]);
          //if(id == 0) printf("JD count = %d\n", count);
  
          if(k == num_ev)
          {     
            double temp;
            int index[k], temp_index;
            for(i =0; i < k; i++)
            {
              index[i] = i;
            }
            
            if(band_type == 1)
            {
              for(i =0; i < k-1; i++)
              {
                for(j = 0; j < k-1-i; j++)
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
              for(i =0; i < k-1; i++)
              {
                for(j = 0; j < k-1-i; j++)
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
         
            for(i =0; i < k; i++)
            {
              memcpy(&eigen_vec_out[i*n_ham], &X_bar[index[i]*n_ham], sizeof(double complex)*n_ham);
            }
            //eigen_vec_out = X_bar;
            //lambda_out = lambda;
            //memcpy ( eigen_vec_out, X_bar, k*n_ham*sizeof(double complex));
            memcpy ( lambda_out, lambda, k*sizeof(double));
            //cudaProfilerStop();
            //cudaDeviceReset();
            for(i =0; i < k; i++){
            if(id == 0) printf("lambda[%d] = %f \n", i, lambda[i]);
            }
            //if(id == 0) printf("Total JD count = %d\n", count);
            //if(id == 0) printf("Total GMRES iterations = %d\n", * ls_counter);

//             cudaFree(colptr_real);
//             cudaFree(rowptr_real);
//             cudaFree(valptr_real);
//             cudaFree(colptr_img);
//             cudaFree(rowptr_img);
//             cudaFree(valptr_img);
//             cudaFree(t_device);
//             cudaFree(scalar1_device);
//             cudaFree(scalar2_device);
//             cudaFree(y_device);
//             cudaFree(s_v_device);
//             cudaFree(norm_r_device);
//             cudaFree(norm_vct_device);
//             cudaFree(norm_t_in_device);
//             cudaFree(s_u_device);
//             cudaFree(norm_t_device);
//             cudaFree(norm_V0_device);
//             cudaFree(r_device);
//             cudaFree(dot_temp_device);
//             cudaFree(w_device);
//             cudaFree(v_device);
//             //cudaFree(vct_device);
//             cudaFree(M_device);
//             cudaFree(temp1_device);
//             cudaFree(eigen_vec_device);
//             cudaFree(eigen_val_device);
//             //cudaFree(ubar_device);
//             //cudaFree(w_bar_device);
//             cudaFree(u_device);
//             //cudaFree(temp1_device);
//             //cudaFree(temp2_device);
//             cudaFree(eigen_vec_device);
//             //cudaFree(Q_bar_device);
//             //cudaFree(X_bar_device);
//             
//             cudaStreamDestroy(stream1);
//             cudaStreamDestroy(stream2);
//             cudaStreamDestroy(stream3);
//             cudaStreamDestroy(stream4);
//             
//             cudaFreeHost(v);
//             cudaFreeHost(w);
//             cudaFreeHost(temp_host_v);
//             cudaFreeHost(M);
//             
//             cublasDestroy(cubl_handle);
//             cusparseDestroy(cusp_handle);
//             cusparseDestroyMatDescr(cusp_descra);
//             
//             free(M_out);
//             free(isuppz);
//             free(s_u);
//             free(s_v);
//             free(scalar1);
//             free(scalar2);
//             free(norm_r);
//             free(norm_t);
//             free(u);
//             free(norm_vct);
//             free(norm_t_in);
//             free(ls_counter);
//             free(reqs);
//             free(shift_init_Mi);
//             free(shift_end_Mi);
//             free(norm_V0);
//             free(norm_temp);
//             free(y);
//             free(dot_temp);
//             free(dot_temp_1);
//             free(dot_temp_2);
//             free(V0);
//             free(X_bar);
//             free(eigen_vec);
//             free(eigen_val);
//             free(temp_ev);
//             
//             cudaDeviceSynchronize();
//             cudaDeviceReset();
            
            //printf(" i m here 15\n");
 
            return; // call cleanup
          }
  
          m = m -1;
  
          // no need of keeping M bcoz its calculated on every run
          //M = [];
  
          // M is kept only on CPU here
          //M = (double complex *) calloc(m*m, sizeof(double complex));
          temp1 = (double complex *) malloc(n_ham*sizeof(double complex));
          memset(temp1, 0, n_ham*sizeof(double complex));
          temp2 = (double complex *) malloc(n_ham*sizeof(double complex));
          memset(temp2, 0, n_ham*sizeof(double complex));
          temp_host_v = (double complex *) malloc(n_ham*jd_max_step*sizeof(double complex));
          memset(temp_host_v, 0, n_ham*jd_max_step*sizeof(double complex));
          temp_host_w = (double complex *) malloc(n_ham*jd_max_step*sizeof(double complex));
          memset(temp_host_w, 0, n_ham*jd_max_step*sizeof(double complex));

          for(i = 0; i < m; i++)
          {
            //v_temp(:,i) = v * eig_mat(:,i+1);
  
            cblas_zgemv(CblasColMajor, CblasNoTrans, n_ham, m+1, scalar_1, v, n_ham, &eigen_vec[(i+1)*(m+1)], 1, scalar_2, temp1, 1);

            // w_temp(:,i) = w * eig_mat(:,i+1);
    
            cblas_zgemv(CblasColMajor, CblasNoTrans, n_ham, m+1, scalar_1, w, n_ham, &eigen_vec[(i+1)*(m+1)], 1, scalar_2, temp2, 1);

            if(id ==0)
            memcpy(&temp_host_v[(i*n_ham)+shift_init_Mi[id]-1], &temp1[shift_init_Mi[id]-1],(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high)*sizeof(double complex)));
            else if(id != 0 && id != num_procs-1)
            memcpy(&temp_host_v[(i*n_ham)+shift_init_Mi[id]-1-overlap_low], &temp1[shift_init_Mi[id]-1-overlap_low],(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol)*sizeof(double complex)));
            else
            memcpy(&temp_host_v[(i*n_ham)+shift_init_Mi[id]-1-overlap_low], &temp1[shift_init_Mi[id]-1-overlap_low], (((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low)*sizeof(double complex)));
  
            if(id ==0)
            memcpy(&temp_host_w[(i*n_ham)+shift_init_Mi[id]-1], &temp2[shift_init_Mi[id]-1],(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high)*sizeof(double complex)));
            else if(id != 0 && id != num_procs-1)
            memcpy(&temp_host_w[(i*n_ham)+shift_init_Mi[id]-1-overlap_low], &temp2[shift_init_Mi[id]-1-overlap_low],(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol)*sizeof(double complex)));
            else
            memcpy(&temp_host_w[(i*n_ham)+shift_init_Mi[id]-1-overlap_low], &temp2[shift_init_Mi[id]-1-overlap_low], (((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low)*sizeof(double complex)));

          }
          
  
          // v = v_temp(:,(1:m));
          memcpy(v, temp_host_v,(sizeof(double complex)*n_ham*m));
          
//           if(id ==0)
//           memcpy(&v[shift_init_Mi[id]-1], &v[shift_init_Mi[id]-1],(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high)*sizeof(double complex)));
//           else if(id != 0 && id != num_procs-1)
//           memcpy(&v[shift_init_Mi[id]-1-overlap_low], &v[shift_init_Mi[id]-1-overlap_low],(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol)*sizeof(double complex)));
//           else
//           memcpy(&v[shift_init_Mi[id]-1-overlap_low], &v[shift_init_Mi[id]-1-overlap_low],(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low)*sizeof(double complex)));
  
          // w = w_temp(:,(1:m)); 
          memcpy(w, temp_host_w,(sizeof(double complex)*n_ham*m));
  
          free(temp1);
          free(temp2);
          free(temp_host_v);
          free(temp_host_w);
  
          MPI_Barrier(upt_comm);
  
          // M(i,i) = eig_val(i+1);
          //M[(i*m)+i].x = eigen_val[i+1];
  
          for(i = 0; i < m; i++)
          {
            eigen_val[i] = eigen_val[i+1];
  
          }
  
          free(eigen_vec);
  
          MPI_Barrier(upt_comm);
  
  
          eigen_vec = (double complex *) calloc(m*m, sizeof(double complex));
          
          for(i = 0; i < m; i++)
          {
            // eig_mat(:,i) = e(:,i);
            eigen_vec[(i*m)+i] = 1.0000 + 0.0000 * I;
          }
          
          //printf(" i m here 18\n");
          
          // s_u = norm(v(:,1));
          //c_zdotu((shift_end_Mi[id]-shift_init_Mi[id])+1, &v[shift_init_Mi[id]-1], &v[shift_init_Mi[id]-1], dot_temp);
          cblas_zdotc_sub((shift_end_Mi[id]-shift_init_Mi[id])+1, &v[shift_init_Mi[id]-1], 1, &v[shift_init_Mi[id]-1], 1, dot_temp);
  
          MPI_Barrier(upt_comm);
  
          MPI_Allreduce(dot_temp, dot_temp_2, 1 , MPI_DOUBLE_COMPLEX, MPI_SUM, upt_comm);
  
          *s_u = sqrt(creal(dot_temp_2[0]));
  
          // s_v = eig_val(1) / (s_u * s_u); 
          *s_v = eigen_val[0]/((*s_u) * (*s_u));
          
          // u = v(:,1)/s_u;
          if(id == 0)
          {
            vct_div_slr(&v[shift_init_Mi[id]-1], &u[shift_init_Mi[id]-1], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high, *s_u);
          }
          else if(id != 0 && id != num_procs-1)
          {
            vct_div_slr(&v[shift_init_Mi[id]-1-overlap_low], &u[shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol, *s_u);
          }
          else
          {
            vct_div_slr(&v[shift_init_Mi[id]-1-overlap_low], &u[shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low, *s_u);
          }
  
          // r = (w(:,1) / s_u) - (s_v * u); 
          if(id == 0)
          {
            vct_div_slr_minus_scl_vct(&w[shift_init_Mi[id]-1], &u[shift_init_Mi[id]-1], &r[shift_init_Mi[id]-1], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high, *s_u, *s_v);
          }
          else if(id != 0 && id != num_procs-1)
          {
            vct_div_slr_minus_scl_vct(&w[shift_init_Mi[id]-1-overlap_low], &u[shift_init_Mi[id]-1-overlap_low], &r[shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol, *s_u, *s_v);
          }
          else
          {
            vct_div_slr_minus_scl_vct(&w[shift_init_Mi[id]-1-overlap_low], &u[shift_init_Mi[id]-1-overlap_low], &r[shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low, *s_u, *s_v);
          }
  
          //printf(" i m here 19\n");
          
          // norm(r_bar);
          //c_zdotu((shift_end_Mi[id]-shift_init_Mi[id])+1, &r[shift_init_Mi[id]-1], &r[shift_init_Mi[id]-1], dot_temp);
          cblas_zdotc_sub((shift_end_Mi[id]-shift_init_Mi[id])+1, &r[shift_init_Mi[id]-1], 1, &r[shift_init_Mi[id]-1], 1, dot_temp);
  
          MPI_Barrier(upt_comm);
  
          MPI_Allreduce(dot_temp, dot_temp_2, 1 , MPI_DOUBLE_COMPLEX, MPI_SUM, upt_comm);
  
          *norm_r = sqrt(creal(dot_temp_2[0]));

        }// end while(norm(r) < tol)
  
        
        // restart
        if (m >= jd_max_step) 
        {
          // M is kept only on CPU here
          //M = (double complex *) calloc(jd_min_step*jd_min_step, sizeof(double complex));
  
          temp1 = (double complex *) malloc(n_ham*sizeof(double complex));
          memset(temp1, 0, n_ham*sizeof(double complex));
          temp2 = (double complex *) malloc(n_ham*sizeof(double complex));
          memset(temp2, 0, n_ham*sizeof(double complex));
          temp_host_v = (double complex *) malloc(n_ham*jd_max_step*sizeof(double complex));
          memset(temp_host_v, 0, n_ham*jd_max_step*sizeof(double complex));
          temp_host_w = (double complex *) malloc(n_ham*jd_max_step*sizeof(double complex));
          memset(temp_host_w, 0, n_ham*jd_max_step*sizeof(double complex));
         
          // for i = 1 : m_min
          for(i = 0; i < jd_min_step; i++)
          {
            //v_temp(:,i) = v * eig_mat(:,i);
	    cblas_zgemv(CblasColMajor, CblasNoTrans, n_ham, m, scalar_1, v, n_ham, &eigen_vec[i*m], 1, scalar_2, temp1, 1);
            
            // w_temp(:,i) = w * eig_mat(:,i);
            cblas_zgemv(CblasColMajor, CblasNoTrans, n_ham, m, scalar_1, w, n_ham, &eigen_vec[i*m], 1, scalar_2, temp2, 1);


            if(id ==0)
            memcpy(&temp_host_v[(i*n_ham)+shift_init_Mi[id]-1], &temp1[shift_init_Mi[id]-1],(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high)*sizeof(double complex)));
            else if(id != 0 && id != num_procs-1)
            memcpy(&temp_host_v[(i*n_ham)+shift_init_Mi[id]-1-overlap_low], &temp1[shift_init_Mi[id]-1-overlap_low],(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol)*sizeof(double complex)));
            else
            memcpy(&temp_host_v[(i*n_ham)+shift_init_Mi[id]-1-overlap_low], &temp1[shift_init_Mi[id]-1-overlap_low], (((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low)*sizeof(double complex)));
  
            if(id ==0)
            memcpy(&temp_host_w[(i*n_ham)+shift_init_Mi[id]-1], &temp2[shift_init_Mi[id]-1],(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high)*sizeof(double complex)));
            else if(id != 0 && id != num_procs-1)
            memcpy(&temp_host_w[(i*n_ham)+shift_init_Mi[id]-1-overlap_low], &temp2[shift_init_Mi[id]-1-overlap_low],(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol)*sizeof(double complex)));
            else
            memcpy(&temp_host_w[(i*n_ham)+shift_init_Mi[id]-1-overlap_low], &temp2[shift_init_Mi[id]-1-overlap_low], (((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low)*sizeof(double complex)));
          }
        
         // v = v_temp(:,(1:m_min));
        memcpy(v, temp_host_v,(n_ham*jd_min_step*sizeof(double complex)));
        memcpy(w, temp_host_w,(n_ham*jd_min_step*sizeof(double complex)));
  
        MPI_Barrier(upt_comm);

//         if(id ==0)
//         memcpy(&v[shift_init_Mi[id]-1], &v[shift_init_Mi[id]-1],(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high)*sizeof(double complex)));
//         else if(id != 0 && id != num_procs-1)
//         memcpy(&v[shift_init_Mi[id]-1-overlap_low], &v[shift_init_Mi[id]-1-overlap_low],(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol)*sizeof(double complex)));
//         else
//         memcpy(&v[shift_init_Mi[id]-1-overlap_low], &v[shift_init_Mi[id]-1-overlap_low],(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low)*sizeof(double complex)));
  
        // w = w_temp(:,(1:m_min)); // check the next line
        //cudaStat2 = cudaMemcpyAsync(w, w,(size_t)(n_ham*jd_min_step*sizeof(double complex)), cudaMemcpyHostToDevice, stream2);
  
        free(temp1);
        free(temp2);
        free(temp_host_v);
        free(temp_host_w);
  
        m = jd_min_step;

        } // end if(restart)

        MPI_Barrier(upt_comm);

        //cudaStat1 = cudaMemcpyAsync(w, w,(size_t)(sizeof(double complex)*n_ham*m), cudaMemcpyDeviceToHost, stream2);

        free(M);
        free(eigen_val);
        free(eigen_vec);

        //c_zdotu((shift_end_Mi[id]-shift_init_Mi[id])+1, &r[shift_init_Mi[id]-1], &r[shift_init_Mi[id]-1], dot_temp);
        cblas_zdotc_sub((shift_end_Mi[id]-shift_init_Mi[id])+1, &r[shift_init_Mi[id]-1], 1, &r[shift_init_Mi[id]-1], 1, dot_temp);
  
        MPI_Barrier(upt_comm);
  
        MPI_Allreduce(dot_temp, dot_temp_2, 1 , MPI_DOUBLE_COMPLEX, MPI_SUM, upt_comm);
  
        *norm_r = sqrt(creal(dot_temp_2[0]));
  
        
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
      shift_A(valptr_real, rowptr_real, colptr_real, (shift_end_Mi[id]-shift_init_Mi[id])+1, tol_shift, shift_init_Mi[id]-1);

  
        if(k !=0)
        {
          Q_bar = (double complex *) malloc(sizeof(double complex)*n_ham*(k+1));
          memset(Q_bar, 0, n_ham*(k+1)*sizeof(double complex));
          memcpy (Q_bar, X_bar, sizeof(double complex)*n_ham*k);
          memcpy (&Q_bar[k*n_ham], u, sizeof(double complex)*n_ham);
        }
        else
        {
          Q_bar = (double complex *) malloc(n_ham*sizeof(double complex));
          memset(Q_bar, 0, n_ham*sizeof(double complex));
          memcpy (Q_bar, u, n_ham*sizeof(double complex));
        }


      MPI_Barrier(upt_comm); 

  
      GMRES_CPU(valptr_real, rowptr_real, colptr_real, valptr_img, rowptr_img, colptr_img, n_ham, size_mat_real, r, ls_tol, ls_restart, ls_maxit, Q_bar, t, k, ls_counter, shift_init_Mi, shift_end_Mi, overlap_high, overlap_low, num_procs, id,  upt_comm);
  
      MPI_Barrier(upt_comm); 
      
  
      free(Q_bar);
  
  
      if(k!=0) // MGS after linear solver
      {
  
        for(i = 0; i < k; i++)
        {
          // y = X_bar(:,i)' * t;
          //c_zdotu((shift_end_Mi[id]-shift_init_Mi[id])+1, &X_bar[i*n_ham+shift_init_Mi[id]-1], &t[shift_init_Mi[id]-1], dot_temp); 
          cblas_zdotc_sub((shift_end_Mi[id]-shift_init_Mi[id])+1, &X_bar[i*n_ham+shift_init_Mi[id]-1], 1, &t[shift_init_Mi[id]-1], 1, dot_temp);
  
          MPI_Barrier(upt_comm); 
  
          MPI_Allreduce(dot_temp, y, 1 , MPI_DOUBLE_COMPLEX, MPI_SUM, upt_comm);
  
  
          // t = t - y * X_bar(:,i);
          if(id == 0)
          {
            vct1_sub_mul_vct( &t[shift_init_Mi[id]-1], &X_bar[i*n_ham+shift_init_Mi[id]-1], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high, y);
          }
          else if(id != 0 && id != num_procs-1)
          { 
            vct1_sub_mul_vct( &t[shift_init_Mi[id]-1-overlap_low], &X_bar[i*n_ham+shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol, y);
          }
          else
          {
            vct1_sub_mul_vct( &t[shift_init_Mi[id]-1-overlap_low], &X_bar[i*n_ham+shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low, y);
          }
  
        }
      }
  
      // A = (A - shift * I)
      tol_shift = teta-shift;
      shift_A(valptr_real, rowptr_real, colptr_real, (shift_end_Mi[id]-shift_init_Mi[id])+1, tol_shift, shift_init_Mi[id]-1);
  
      if (num_procs > 1)
      {
        int rank;
        for(rank =0; rank <=num_procs-2; rank++)
        {
          if(id == rank+1)    MPI_Irecv(&t[shift_end_Mi[rank]-overlap_low], overlap_low,MPI_DOUBLE_COMPLEX,rank, 111+rank,upt_comm,&reqs[rank+1]);
          if(id == rank)      MPI_Isend(&t[shift_end_Mi[rank]-overlap_low], overlap_low,MPI_DOUBLE_COMPLEX,rank+1, 111+rank,upt_comm, &reqs[rank+1]);
  
          if (id == rank+1) MPI_Wait(&reqs[rank+1], &status);
          if (id == rank)   MPI_Wait(&reqs[rank+1], &status);
        }
      }

      if (num_procs > 1)
      {
        int rank;
        for(rank =0; rank <=num_procs-2; rank++)
        {
          if(id == rank)    MPI_Irecv(&t[shift_init_Mi[rank+1]-1],overlap_high, MPI_DOUBLE_COMPLEX,rank+1,8+rank,upt_comm, &reqs[rank+1]);
          if(id == rank+1)  MPI_Isend(&t[shift_init_Mi[rank+1]-1],overlap_high, MPI_DOUBLE_COMPLEX,rank,8+rank,upt_comm, &reqs[rank+1]);
  
          if (id == rank+1) MPI_Wait(&reqs[rank+1], &status);
          if (id == rank)   MPI_Wait(&reqs[rank+1], &status);
        }
      }
  
  count = count+1;

  } // end while(k < num_ev)

  free(shift_init_Mi);
  free(shift_end_Mi);
  free(shift_init_M);
  free(shift_end_M);
  
  printf("Total JD count = %d\n", count);
  printf("Total GMRES iterations = %d\n", *ls_counter); 
    
    //cudaProfilerStop();
 
 }// end of JD
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
     void GMRES_CPU(float *valptr_real, int *rowptr_real, int *colptr_real, float *valptr_img, int *rowptr_img, int *colptr_img, int n_ham, int size_mat_real, double complex * r, double ls_tol, int restart, int maxit, double complex *Q_bar, double complex *t, int k, int *ls_counter, int * shift_init_Mi, int * shift_end_Mi, int overlap_high, int overlap_low, int num_procs, int id,  MPI_Comm upt_comm)
     {
       double complex *X0, *Xmin, *dot, *P;
       double complex *rq, *h_min, *g;
       double complex sqrt_temp, sqrt_complex_temp;
       double complex *test;
       double *norm_r_ls, *beta, *norm_temp;
       double complex *v_ls, *w_ls, *g_min, *minimizer, *h; 
       double complex *scalar1_ls, *scalar2_ls;
       double complex *dot_temp_1, *dot_temp_2, *dot_temp, *dot_temp_ls;
       lapack_int *jpvt;
       lapack_int error_code=0;
       double complex coso, sino, scalar_1, scalar_2, complex_temp1, complex_temp2, complex_temp;
       double complex * temp, * temp1, * temp2;
       int i, j;
     
       MPI_Request * reqs_ls;
       MPI_Status status_ls;
     
       reqs_ls=(MPI_Request *) malloc(sizeof(MPI_Request)*(num_procs+1));

       scalar_1 = 1.00 + 0.00 * I;

       scalar_2 = 0.00 + 0.00 * I;
//     
//       scalar1_ls = (double complex *) malloc(sizeof(double complex));
//       scalar2_ls = (double complex *) malloc(sizeof(double complex));
//     
//       scalar1_ls[0].x = 1.00;
//       scalar1_ls[0].y = 0.00;
//       scalar2_ls[0].x = 0.00;
//       scalar2_ls[0].y = 0.00;
//     
       const double complex jcmpx= 0.00000 + 1.00000 * I;

	   int overlap_tol = overlap_high + overlap_low;
 
       X0 = (double complex *) malloc(sizeof(double complex)*n_ham);
       Xmin= (double complex *) malloc(sizeof(double complex)*n_ham); 
       memset(X0, 0, n_ham*sizeof(double complex));
       memset(Xmin, 0, n_ham*sizeof(double complex));
       temp = (double complex *) malloc(sizeof(double complex)*n_ham);
       memset(temp, 0, n_ham*sizeof(double complex));  
       temp2 = (double complex *) malloc(sizeof(double complex)*n_ham);
       memset(temp2, 0, n_ham*sizeof(double complex));
       rq = (double complex *) malloc(sizeof(double complex)*n_ham);
       memset(rq, 0, n_ham*sizeof(double complex));
       v_ls = (double complex *) malloc(sizeof(double complex)*n_ham*(restart+1));
       memset(v_ls, 0, sizeof(double complex)*n_ham*(restart+1));
       w_ls = (double complex *) malloc(sizeof(double complex)*n_ham*restart);
       memset(w_ls, 0, sizeof(double complex)*n_ham*restart);
       h = (double complex *) malloc(sizeof(double complex)*(restart+1)*restart);
       memset(h, 0, sizeof(double complex)*(restart+1)*restart);
       g = (double complex *) malloc(sizeof(double complex)*(restart+1));
       memset(g, 0, sizeof(double complex)*(restart+1));
       minimizer = (double complex *) malloc(sizeof(double complex)*restart);
       memset(minimizer, 0, sizeof(double complex)*restart);
       norm_r_ls = (double *) malloc(sizeof(double));
       scalar1_ls = (double complex *) malloc(sizeof(double complex)); 
       scalar2_ls = (double complex *) malloc(sizeof(double complex));  
       temp1 = (double complex *) malloc(n_ham*sizeof(double complex));
       memset(temp1, 0, n_ham*sizeof(double complex));
       dot_temp_ls = (double complex *) malloc(sizeof(double complex));

       *scalar1_ls = 1.0 + 0.0 * I;
       *scalar2_ls = 0.0 + 0.0 * I;

       double complex * h_temp, *g_temp;
       h_temp = (double complex *) calloc((restart+1)*(restart), sizeof(double complex));
       g_temp = (double complex *) calloc((restart+1), sizeof(double complex));

       // r_ls = -r
       if(id == 0)
       {
         vct1_neg_asg_vct2( &r[shift_init_Mi[id]-1], &r[shift_init_Mi[id]-1], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high);
       }
       else if(id != 0 && id != num_procs-1)
       { 
         vct1_neg_asg_vct2( &r[shift_init_Mi[id]-1-overlap_low], &r[shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol);
       }
       else
       {
         vct1_neg_asg_vct2( &r[shift_init_Mi[id]-1-overlap_low], &r[shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low);
       }
    
       double norm_min =1000000000000;
       double rcond;
       lapack_int *eff_rank_h_min;
          
       dot = (double complex *) malloc(sizeof(double complex));
       dot_temp = (double complex *) malloc(sizeof(double complex));
       dot_temp_1 = (double complex *) malloc(sizeof(double complex));
       dot_temp_2 = (double complex *) malloc(sizeof(double complex));
       norm_temp = (double *) malloc(sizeof(double));
       //norm_r_ls = (double *) malloc(sizeof(double));
       beta = (double *) malloc(sizeof(double));
       h_min = (double complex *) malloc(restart*restart*sizeof(double complex));
       g_min = (double complex *) malloc(restart*sizeof(double complex));
       test = (double complex *) malloc(sizeof(double complex));
       jpvt = (lapack_int *) malloc(restart*sizeof(lapack_int));
       eff_rank_h_min = (lapack_int *) malloc(sizeof(lapack_int));
          
       if(maxit > n_ham)
       {
           maxit = n_ham;
       }
     
       if (restart > n_ham)
       {
           restart = n_ham;
       }
     
       //tol =tol * norm(b);
       //c_zdotu((shift_end_Mi[id]-shift_init_Mi[id])+1, &r[shift_init_Mi[id]-1], &r[shift_init_Mi[id]-1], dot_temp_ls);
       cblas_zdotc_sub((shift_end_Mi[id]-shift_init_Mi[id])+1, &r[shift_init_Mi[id]-1], 1, &r[shift_init_Mi[id]-1], 1, dot_temp_ls);
     
       MPI_Barrier(upt_comm);
     
       MPI_Allreduce(dot_temp_ls, dot_temp_2, 1 , MPI_DOUBLE_COMPLEX, MPI_SUM, upt_comm);
     
       *norm_r_ls = sqrt(creal(dot_temp_2[0]));
         
       ls_tol = (*norm_r_ls) * ls_tol;
     
       int restart_count = 0;
       
       while(restart_count < maxit)
       {
         // temp = x0;
// //         if(id ==0)
// //         cudaStat1 = cudaMemcpy(&temp_device[shift_init_Mi[id]-1], &X0_device[shift_init_Mi[id]-1],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high)*sizeof(double complex)), cudaMemcpyDeviceToDevice);
// //         else if(id != 0 && id != num_procs-1)
// //         cudaStat1 = cudaMemcpy(&temp_device[shift_init_Mi[id]-1-overlap_low], &X0_device[shift_init_Mi[id]-1-overlap_low],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol)*sizeof(double complex)), cudaMemcpyDeviceToDevice);
// //         else
// //         cudaStat1 = cudaMemcpy(&temp_device[shift_init_Mi[id]-1-overlap_low], &X0_device[shift_init_Mi[id]-1-overlap_low],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low)*sizeof(double complex)), cudaMemcpyDeviceToDevice);
// //     
// //         for(int i = 0; i < k+1; i++)
// //         {
// //           cubl_status_ls = cublasZdotc(cubl_handle_ls, (shift_end_Mi[id]-shift_init_Mi[id])+1, &Q_bar_device[i*n_ham+shift_init_Mi[id]-1], 1, &X0_device[shift_init_Mi[id]-1], 1, dot_device);
// //     
// //           cudaStat1 = cudaMemcpy(dot_temp, dot_device,(size_t)(sizeof(double complex)), cudaMemcpyDeviceToHost);
// //     
// //           MPI_Allreduce(dot_temp, dot, 1 , MPI_DOUBLE_COMPLEX, MPI_SUM, upt_comm);
// //     
// //           cudaStat1 = cudaMemcpy(dot_device, dot,(size_t)(sizeof(double complex)), cudaMemcpyHostToDevice);
// //     
// //     
// //           if(id == 0)
// //           {
// //             vct1_sub_mul_vct_kernel<<<numBlock, threadPerBlock>>>( &temp_device[shift_init_Mi[id]-1], &Q_bar_device[i*n_ham+shift_init_Mi[id]-1], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high, dot_device);
// //           }
// //           else if(id != 0 && id != num_procs-1)
// //           { 
// //             vct1_sub_mul_vct_kernel<<<numBlock, threadPerBlock>>>( &temp_device[shift_init_Mi[id]-1-overlap_low], &Q_bar_device[i*n_ham+shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol, dot_device);
// //           }
// //           else
// //           {
// //             vct1_sub_mul_vct_kernel<<<numBlock, threadPerBlock>>>( &temp_device[shift_init_Mi[id]-1-overlap_low], &Q_bar_device[i*n_ham+shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low, dot_device);
// //           }
// //         }
// 
     
         spmv_csr_hybrid((shift_end_Mi[id]-shift_init_Mi[id])+1, rowptr_real, colptr_real, valptr_real, X0, temp2, shift_init_Mi[id]-1);
     
         spmv_csr_hybrid((shift_end_Mi[id]-shift_init_Mi[id])+1, rowptr_img, colptr_img, valptr_img, X0, temp1, shift_init_Mi[id]-1);

         MPI_Barrier(upt_comm);
     
         vct_pls_scl_mul_vct(&temp2[shift_init_Mi[id]-1], &temp1[shift_init_Mi[id]-1], jcmpx, (shift_end_Mi[id]-shift_init_Mi[id])+1);

         MPI_Barrier(upt_comm);
     
         if (num_procs > 1)
         {
           int rank;
           for(rank =0; rank <=num_procs-2; rank++)
           {
             if(id == rank+1)    MPI_Irecv(&temp2[shift_end_Mi[rank]-overlap_low], overlap_low, MPI_DOUBLE_COMPLEX, rank, 411+rank, upt_comm, &reqs_ls[rank+1]);
             if(id == rank)      MPI_Isend(&temp2[shift_end_Mi[rank]-overlap_low], overlap_low, MPI_DOUBLE_COMPLEX, rank+1, 411+rank, upt_comm, &reqs_ls[rank+1]);
     
             if (id == rank+1) MPI_Wait(&reqs_ls[rank+1], &status_ls);
             if (id == rank)   MPI_Wait(&reqs_ls[rank+1], &status_ls);
           }
         }

         MPI_Barrier(upt_comm);
 
         if (num_procs > 1)
         {
	   int rank;
           for(rank =0; rank <=num_procs-2; rank++)
           {
             if(id == rank)    MPI_Irecv(&temp2[shift_init_Mi[rank+1]-1],overlap_high, MPI_DOUBLE_COMPLEX,rank+1,450+rank,upt_comm, &reqs_ls[rank+1]);
             if(id == rank+1)  MPI_Isend(&temp2[shift_init_Mi[rank+1]-1],overlap_high, MPI_DOUBLE_COMPLEX,rank,450+rank,upt_comm, &reqs_ls[rank+1]);
     
             if (id == rank+1) MPI_Wait(&reqs_ls[rank+1], &status_ls);
             if (id == rank)   MPI_Wait(&reqs_ls[rank+1], &status_ls);
           }
         }
     
         MPI_Barrier(upt_comm);
 
       // temp = temp2;
       if(id ==0)
       memcpy(&temp[shift_init_Mi[id]-1], &temp2[shift_init_Mi[id]-1],(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high)*sizeof(double complex)));
       else if(id != 0 && id != num_procs-1)
       memcpy(&temp[shift_init_Mi[id]-1-overlap_low], &temp2[shift_init_Mi[id]-1-overlap_low],(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol)*sizeof(double complex)));
       else
       memcpy(&temp[shift_init_Mi[id]-1-overlap_low], &temp2[shift_init_Mi[id]-1-overlap_low],(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low)*sizeof(double complex)));
     
	   
       for(i = 0; i < k+1; i++)
         {
           // temp = temp - (dot(u(:,t),temp2)) * u(:,t);
           //c_zdotu((shift_end_Mi[id]-shift_init_Mi[id])+1, &Q_bar[i*n_ham+shift_init_Mi[id]-1], &temp2[shift_init_Mi[id]-1], dot_temp);
           cblas_zdotc_sub((shift_end_Mi[id]-shift_init_Mi[id])+1, &Q_bar[i*n_ham+shift_init_Mi[id]-1], 1, &temp2[shift_init_Mi[id]-1], 1, dot_temp);
     
           MPI_Barrier(upt_comm);

           MPI_Allreduce(dot_temp, dot, 1 , MPI_DOUBLE_COMPLEX, MPI_SUM, upt_comm);

           //printf("dot %f %f \n", creal(dot[0]), cimag(dot[0]));
   
           if(id == 0)
           {
             vct1_sub_mul_vct( &temp[shift_init_Mi[id]-1], &Q_bar[i*n_ham+shift_init_Mi[id]-1], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high, dot);
           }
           else if(id != 0 && id != num_procs-1)
           { 
             vct1_sub_mul_vct( &temp[shift_init_Mi[id]-1-overlap_low], &Q_bar[i*n_ham+shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol, dot);
           }
           else
           {
             vct1_sub_mul_vct( &temp[shift_init_Mi[id]-1-overlap_low], &Q_bar[i*n_ham+shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low, dot);
           }
         }
     
         //rq = b - temp2;     here b = r_device_ls
         if(id == 0)
         {
           vct1_sub_vct2_asg_vct3( &r[shift_init_Mi[id]-1], &temp[shift_init_Mi[id]-1], &rq[shift_init_Mi[id]-1], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high);
         }
         else if(id != 0 && id != num_procs-1)
         { 
           vct1_sub_vct2_asg_vct3( &r[shift_init_Mi[id]-1-overlap_low], &temp[shift_init_Mi[id]-1-overlap_low], &rq[shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol);
         }
         else
         {
           vct1_sub_vct2_asg_vct3( &r[shift_init_Mi[id]-1-overlap_low], &temp[shift_init_Mi[id]-1-overlap_low], &rq[shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low);
         }
     
         // beta=norm(rq);
         //c_zdotu((shift_end_Mi[id]-shift_init_Mi[id])+1, &rq[shift_init_Mi[id]-1], &rq[shift_init_Mi[id]-1], dot_temp_ls);
         cblas_zdotc_sub((shift_end_Mi[id]-shift_init_Mi[id])+1, &rq[shift_init_Mi[id]-1], 1, &rq[shift_init_Mi[id]-1], 1, dot_temp_ls);
     
         MPI_Barrier(upt_comm);
     
         MPI_Allreduce(dot_temp_ls, dot_temp_2, 1 , MPI_DOUBLE_COMPLEX, MPI_SUM, upt_comm);
     
         *beta = sqrt(creal(dot_temp_2[0]));
     
     
         // v(:,1)=rq/beta; 
         if(id == 0)
         {
           vct_div_slr( &rq[shift_init_Mi[id]-1], &v_ls[shift_init_Mi[id]-1], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high, *beta);
         }
         else if(id != 0 && id != num_procs-1)
         { 
           vct_div_slr( &rq[shift_init_Mi[id]-1-overlap_low], &v_ls[shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol, *beta);
         }
         else
         {
           vct_div_slr( &rq[shift_init_Mi[id]-1-overlap_low], &v_ls[shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low, *beta);
         }
     
         memset(h, 0, sizeof(double complex)*restart*(restart+1));
     
         for(j = 0; j < restart; j++)
         {
           // temp = v(:,j);
// //           if(id ==0)
// //           cudaStat1 = cudaMemcpy(&temp_device[shift_init_Mi[id]-1], &v_ls_device[j*n_ham+(shift_init_Mi[id]-1)],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high)*sizeof(double complex)), cudaMemcpyDeviceToDevice);
// //           else if(id != 0 && id != num_procs-1)
// //           cudaStat1 = cudaMemcpy(&temp_device[shift_init_Mi[id]-1-overlap_low], &v_ls_device[j*n_ham+(shift_init_Mi[id]-1-overlap_low)],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol)*sizeof(double complex)), cudaMemcpyDeviceToDevice);
// //           else
// //           cudaStat1 = cudaMemcpy(&temp_device[shift_init_Mi[id]-1-overlap_low], &v_ls_device[j*n_ham+(shift_init_Mi[id]-1-overlap_low)],(size_t)(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low)*sizeof(double complex)), cudaMemcpyDeviceToDevice);
// //     
// //           for(int i = 0; i < k+1; i++)
// //           {
// //           // temp = temp - (dot(Q_bar(:,i),v(:,j))) * Q_bar(:,i); norm(Q_bar) == 1
// //             cubl_status_ls = cublasZdotc(cubl_handle_ls, (shift_end_Mi[id]-shift_init_Mi[id])+1, &Q_bar_device[i*n_ham+shift_init_Mi[id]-1], 1, &v_ls_device[j*n_ham+shift_init_Mi[id]-1], 1, dot_device); 
// //     
// //             cudaStat1 = cudaMemcpy(dot_temp, dot_device,(size_t)(sizeof(double complex)), cudaMemcpyDeviceToHost);
// //     
// //             MPI_Allreduce(dot_temp, dot, 1 , MPI_DOUBLE_COMPLEX, MPI_SUM, upt_comm);
// //     
// //             cudaStat1 = cudaMemcpy(dot_device, dot,(size_t)(sizeof(double complex)), cudaMemcpyHostToDevice);
// //     
// //             if(id == 0)
// //             {
// //               vct1_sub_mul_vct_kernel<<<numBlock, threadPerBlock>>>( &temp_device[shift_init_Mi[id]-1], &Q_bar_device[i*n_ham+shift_init_Mi[id]-1], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high, dot_device);
// //             }
// //             else if(id != 0 && id != num_procs-1)
// //             { 
// //               vct1_sub_mul_vct_kernel<<<numBlock, threadPerBlock>>>( &temp_device[shift_init_Mi[id]-1-overlap_low], &Q_bar_device[i*n_ham+shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol, dot_device);
// //             }
// //             else
// //             {
// //               vct1_sub_mul_vct_kernel<<<numBlock, threadPerBlock>>>( &temp_device[shift_init_Mi[id]-1-overlap_low], &Q_bar_device[i*n_ham+shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low, dot_device);
// //             }
// //           }
// //    
// //          // temp2=A*temp;
// //          // w(:,j)=temp2;

     
         spmv_csr_hybrid((shift_end_Mi[id]-shift_init_Mi[id])+1, rowptr_real, colptr_real, valptr_real, &v_ls[j*n_ham], temp2, shift_init_Mi[id]-1);
     
         spmv_csr_hybrid((shift_end_Mi[id]-shift_init_Mi[id])+1, rowptr_img, colptr_img, valptr_img, &v_ls[j*n_ham], temp1, shift_init_Mi[id]-1);

         MPI_Barrier(upt_comm);
     
         vct_pls_scl_mul_vct(&temp2[shift_init_Mi[id]-1], &temp1[shift_init_Mi[id]-1], jcmpx, (shift_end_Mi[id]-shift_init_Mi[id])+1);

         MPI_Barrier(upt_comm);

     
         if (num_procs > 1)
         {
	   int rank;
           for(rank =0; rank <=num_procs-2; rank++)
           {
             if(id == rank+1)    MPI_Irecv(&temp2[shift_end_Mi[rank]-overlap_low], overlap_low,MPI_DOUBLE_COMPLEX,rank, 311+rank,upt_comm,&reqs_ls[rank+1]);
             if(id == rank)      MPI_Isend(&temp2[shift_end_Mi[rank]-overlap_low], overlap_low,MPI_DOUBLE_COMPLEX,rank+1, 311+rank,upt_comm, &reqs_ls[rank+1]);
     
             if (id == rank+1) MPI_Wait(&reqs_ls[rank+1], &status_ls);
             if (id == rank)   MPI_Wait(&reqs_ls[rank+1], &status_ls);
          }
        }
 
        MPI_Barrier(upt_comm);
     
         if (num_procs > 1)
         {
           int rank;
           for(rank =0; rank <=num_procs-2; rank++)
           {
             if(id == rank)    MPI_Irecv(&temp2[shift_init_Mi[rank+1]-1],overlap_high, MPI_DOUBLE_COMPLEX,rank+1,350+rank,upt_comm, &reqs_ls[rank+1]);
             if(id == rank+1)  MPI_Isend(&temp2[shift_init_Mi[rank+1]-1],overlap_high, MPI_DOUBLE_COMPLEX,rank,350+rank,upt_comm, &reqs_ls[rank+1]);
     
             if (id == rank+1) MPI_Wait(&reqs_ls[rank+1], &status_ls);
             if (id == rank)   MPI_Wait(&reqs_ls[rank+1], &status_ls);
           }
         }
     
         MPI_Barrier(upt_comm);
     
         if(id ==0)
         memcpy(&w_ls[j*n_ham+shift_init_Mi[id]-1], &temp2[shift_init_Mi[id]-1],(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high)*sizeof(double complex)));
         else if(id != 0 && id != num_procs-1)
         memcpy(&w_ls[j*n_ham+shift_init_Mi[id]-1-overlap_low], &temp2[shift_init_Mi[id]-1-overlap_low],(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol)*sizeof(double complex)));
         else
         memcpy(&w_ls[j*n_ham+shift_init_Mi[id]-1-overlap_low], &temp2[shift_init_Mi[id]-1-overlap_low],(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low)*sizeof(double complex)));
     
         for(i = 0; i < k+1; i++)
           {
           // temp = temp - (dot(Q_bar(:,i),w(:,j))) * Q_bar(:,i); norm(Q_bar) == 1
           //c_zdotu((shift_end_Mi[id]-shift_init_Mi[id])+1, &Q_bar[i*n_ham+shift_init_Mi[id]-1], &temp2[shift_init_Mi[id]-1], dot_temp);
           cblas_zdotc_sub((shift_end_Mi[id]-shift_init_Mi[id])+1, &Q_bar[i*n_ham+shift_init_Mi[id]-1], 1, &temp2[shift_init_Mi[id]-1], 1, dot_temp);
    
           MPI_Barrier(upt_comm);
		        
           MPI_Allreduce(dot_temp, dot, 1 , MPI_DOUBLE_COMPLEX, MPI_SUM, upt_comm);

		       
           if(id == 0)
           {
             vct1_sub_mul_vct( &w_ls[j*n_ham+shift_init_Mi[id]-1], &Q_bar[i*n_ham+shift_init_Mi[id]-1], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high, dot);
           }
           else if(id != 0 && id != num_procs-1)
           { 
             vct1_sub_mul_vct( &w_ls[j*n_ham+shift_init_Mi[id]-1-overlap_low], &Q_bar[i*n_ham+shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol, dot);
           }
           else
           {
             vct1_sub_mul_vct( &w_ls[j*n_ham+shift_init_Mi[id]-1-overlap_low], &Q_bar[i*n_ham+shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low, dot);
           }
           }
     
           // w(:,j)=temp2;
		
           for(i = 0; i <= j ; i++)
           {
             // h(i,j)=dot(w(:,j),v(:,i));
             //c_zdotu((shift_end_Mi[id]-shift_init_Mi[id])+1, &w_ls[j*n_ham+shift_init_Mi[id]-1], &v_ls[i*n_ham+shift_init_Mi[id]-1], dot_temp);
             cblas_zdotc_sub((shift_end_Mi[id]-shift_init_Mi[id])+1, &w_ls[j*n_ham+shift_init_Mi[id]-1], 1, &v_ls[i*n_ham+shift_init_Mi[id]-1], 1, dot_temp);
     
             MPI_Barrier(upt_comm);
     
             MPI_Allreduce(dot_temp, &h[(j*(restart+1))+i], 1 , MPI_DOUBLE_COMPLEX, MPI_SUM, upt_comm);

     
             //  w(:,j)=w(:,j)-h(i,j)*v(:,i);
             if(id == 0)
             {
               vct1_sub_mul_vct( &w_ls[j*n_ham+shift_init_Mi[id]-1], &v_ls[i*n_ham+shift_init_Mi[id]-1], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high, &h[(j*(restart+1))+i]);
             }
             else if(id != 0 && id != num_procs-1)
             { 
               vct1_sub_mul_vct( &w_ls[j*n_ham+shift_init_Mi[id]-1-overlap_low], &v_ls[i*n_ham+shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol, &h[(j*(restart+1))+i]);
             }
             else
             {
               vct1_sub_mul_vct( &w_ls[j*n_ham+shift_init_Mi[id]-1-overlap_low], &v_ls[i*n_ham+shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low, &h[(j*(restart+1))+i]);
             }
           }
     
           // h(j+1,j)=norm(w(:,j));
           //c_zdotu((shift_end_Mi[id]-shift_init_Mi[id])+1, &w_ls[j*n_ham+shift_init_Mi[id]-1], &w_ls[j*n_ham+shift_init_Mi[id]-1], dot_temp_ls);
           cblas_zdotc_sub((shift_end_Mi[id]-shift_init_Mi[id])+1, &w_ls[j*n_ham+shift_init_Mi[id]-1], 1, &w_ls[j*n_ham+shift_init_Mi[id]-1], 1, dot_temp_ls);
     
           MPI_Barrier(upt_comm);
     
           MPI_Allreduce(dot_temp_ls, dot_temp_2, 1 , MPI_DOUBLE_COMPLEX, MPI_SUM, upt_comm);
    
           *norm_temp = sqrt(creal(dot_temp_2[0]));

           h[(j*(restart+1))+j+1] = *norm_temp + 0.0000 * I;
     
           //memcpy(creal(h[(j*(restart+1))+j+1]), norm_temp,(sizeof(double))); 
     
           memcpy(test, &h[(j*(restart+1))+j+1],(sizeof(double complex)));

           //if h(j+1,j)==0
           //if(h_device[((j+1)*(restart+1))+j].x == 0.0000)
           if(creal(test[0]) == 0.0000)
             {
             restart=j;
             }
           else
             {
               //v(:,j+1)=w(:,j)/h(j+1,j);
               if(id == 0)
               {
                 vct_1_div_asg_to_vct_2( &w_ls[j*n_ham+shift_init_Mi[id]-1], &v_ls[(j+1)*n_ham+shift_init_Mi[id]-1], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high, creal(h[(j*(restart+1))+j+1]));
               }
               else if(id != 0 && id != num_procs-1)
               { 
                 vct_1_div_asg_to_vct_2( &w_ls[j*n_ham+shift_init_Mi[id]-1-overlap_low], &v_ls[(j+1)*n_ham+shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol, creal(h[(j*(restart+1))+j+1]));
              }
               else
               {
                 vct_1_div_asg_to_vct_2( &w_ls[j*n_ham+shift_init_Mi[id]-1-overlap_low], &v_ls[(j+1)*n_ham+shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low, creal(h[(j*(restart+1))+j+1]));
               }
             }    
       }
             // g(1:m+1,:)=0;
             memset(g, 0, sizeof(double complex)*(restart+1));
 
             // g(1,:)=beta;
             //cudaStat2 = cudaMemcpy(beta, beta_device,(size_t)(sizeof(double)), cudaMemcpyDeviceToHost);
             *test = *beta + 0.0000 * I;
             memcpy(g, test,(sizeof(double complex)));

     
       for(j = 0; j < restart; j++)
       {
         // eigen_vec is kept only on CPU here
         P = (double complex *) calloc((restart+1)*(restart+1), sizeof(double complex));
     
         for(i = 0; i < restart+1; i++)
         {
           // P = eye(restart+1);
           P[(i*(restart+1))+i]= 1.0000 + 0.0000 * I;
         }

     
         // sqrt(h(j+1,j)^2 + h(j,j)^2);
         // h(j+1,j)^2
         complex_temp1 = ((creal(h[(j*(restart+1))+j+1]) * creal(h[(j*(restart+1))+j+1])) - (cimag(h[(j*(restart+1))+j+1]) * cimag(h[(j*(restart+1))+j+1]))) + (((creal(h[(j*(restart+1))+j+1]) * cimag(h[(j*(restart+1))+j+1])) + (cimag(h[(j*(restart+1))+j+1]) * creal(h[(j*(restart+1))+j+1]))) * I);
         //h(j,j)^2
         complex_temp2 = ((creal(h[j*(restart+1)+j]) * creal(h[j*(restart+1)+j])) - (cimag(h[j*(restart+1)+j]) * cimag(h[j*(restart+1)+j]))) + (((creal(h[j*(restart+1)+j]) * cimag(h[j*(restart+1)+j])) + (cimag(h[j*(restart+1)+j]) * creal(h[j*(restart+1)+j]))) * I);
     
         complex_temp = (creal(complex_temp1) + creal(complex_temp2)) + ((cimag(complex_temp1) + cimag(complex_temp2)) * I);

         // sqrt()
	 if(abs(cimag(complex_temp)) != 0.00000000000000000)
	 {
           sqrt_complex_temp = (sqrt((creal(complex_temp) + sqrt(creal(complex_temp)*creal(complex_temp) + cimag(complex_temp)*cimag(complex_temp)))/2)) + (((cimag(complex_temp)/abs(cimag(complex_temp))) * sqrt(((- creal(complex_temp) + sqrt(creal(complex_temp)*creal(complex_temp) + cimag(complex_temp)*cimag(complex_temp)))/2))) * I);
         }
         else
         {
           sqrt_complex_temp = (sqrt((creal(complex_temp) + sqrt(creal(complex_temp)*creal(complex_temp) + cimag(complex_temp)*cimag(complex_temp)))/2)) + (0.0000000000 * I);
         }
         
       
         // sino=h(j+1,j)/(sqrt(h(j+1,j)^2 + h(j,j)^2));
         sino = (((creal(h[(j*(restart+1))+j+1]) * creal(sqrt_complex_temp)) + (cimag(h[(j*(restart+1))+j+1]) * cimag(sqrt_complex_temp))) / ((creal(sqrt_complex_temp) * creal(sqrt_complex_temp)) + (cimag(sqrt_complex_temp) * cimag(sqrt_complex_temp)))) + ((((cimag(h[(j*(restart+1))+j+1]) * creal(sqrt_complex_temp)) - (creal(h[(j*(restart+1))+j+1]) * cimag(sqrt_complex_temp))) / ((creal(sqrt_complex_temp) * creal(sqrt_complex_temp)) + (cimag(sqrt_complex_temp) * cimag(sqrt_complex_temp)))) * I);

     
         // coso=h(j,j)/(sqrt(h(j+1,j)^2 + h(j,j)^2));
         coso = (((creal(h[j*(restart+1)+j]) * creal(sqrt_complex_temp)) + (cimag(h[j*(restart+1)+j]) * cimag(sqrt_complex_temp))) / ((creal(sqrt_complex_temp) * creal(sqrt_complex_temp)) + (cimag(sqrt_complex_temp) * cimag(sqrt_complex_temp)))) + ((((cimag(h[j*(restart+1)+j]) * creal(sqrt_complex_temp)) - (creal(h[j*(restart+1)+j]) * cimag(sqrt_complex_temp))) / ((creal(sqrt_complex_temp) * creal(sqrt_complex_temp)) + (cimag(sqrt_complex_temp) * cimag(sqrt_complex_temp)))) * I);

         // P(j,j)=conj(coso); // I think taking conj is wrong some case dosent converge 
         P[j*(restart+1)+j] = coso;
     
         // P(j+1,j+1)=coso;
         P[(j+1)*(restart+1)+j+1] = creal(coso) + (cimag(coso) * I);
     
         // P(j,j+1)=conj(sino); // I think taking conj is wrong some case dosent converge 
         P[(j+1)*(restart+1)+j] = sino;
     
         // P(j+1,j)=-sino;
         P[(j*(restart+1))+j+1] = - (creal(sino) + (cimag(sino) * I));

         // h=P*h;
         cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, restart+1, restart, restart+1, scalar1_ls, P, restart+1, h, restart+1, scalar2_ls, h_temp, restart+1);
         memcpy (h, h_temp, sizeof(double complex)*(restart+1)*restart);
     
         // g=P*g;
         cblas_zgemv(CblasColMajor, CblasNoTrans, restart+1, restart+1, scalar1_ls, P, restart+1, g, 1, scalar2_ls, g_temp, 1);
         memcpy (g, g_temp, sizeof(double complex)*(restart+1));

         
         free(P);
 
       }

            
       // minimizer=h(1:m,1:m)\g(1:m,1);
       for(i = 0; i < restart; i++) { memcpy (&h_min[restart*i], &h[(restart*i)+i], sizeof(double complex)*restart); }
       memcpy(g_min, g, sizeof(double complex)*restart);
       
       for(i = 0; i < restart; i++)
       {
         jpvt[i] = 0;
       }
     
       error_code = LAPACKE_zgelsy(LAPACK_COL_MAJOR, restart, restart, 1, (MKL_Complex16 *)h_min, restart, (MKL_Complex16 *)g_min, restart, jpvt, rcond, eff_rank_h_min);
      
       // xm=x0+v(:,1:m)*minimizer;
       memcpy(minimizer, g_min,(sizeof(double complex)*restart));

       cblas_zgemv(CblasColMajor, CblasNoTrans, n_ham, restart, scalar1_ls, v_ls, n_ham, minimizer, 1, scalar2_ls, temp, 1);

       //vct1_add_vct2_asg_vct3_kernel<<<numBlock, threadPerBlock>>>(X0_device, temp_device, temp2_device, n_ham);
       if(id == 0)
       {
         vct1_add_vct2_asg_vct3(&X0[shift_init_Mi[id]-1], &temp[shift_init_Mi[id]-1], &temp2[shift_init_Mi[id]-1], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high);
       }
       else if(id != 0 && id != num_procs-1)
       { 
         vct1_add_vct2_asg_vct3(&X0[shift_init_Mi[id]-1-overlap_low], &temp[shift_init_Mi[id]-1-overlap_low], &temp2[shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol);
       }
       else
       {
         vct1_add_vct2_asg_vct3(&X0[shift_init_Mi[id]-1-overlap_low], &temp[shift_init_Mi[id]-1-overlap_low], &temp2[shift_init_Mi[id]-1-overlap_low], (shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low);
       }
       
    
       // if abs(g(m+1,1))<tol 
       if(sqrt(creal(g[restart])*creal(g[restart])+cimag(g[restart])*cimag(g[restart]))<ls_tol)
       { 
         // x = xm;
         memcpy(t, temp2,(sizeof(double complex)*n_ham));
     
         *ls_counter =  *ls_counter + restart_count * restart;
     
     
         free(X0);
         free(Xmin);
         free(dot);
         free(temp);
         free(temp2);
         free(rq);
         free(beta);
         free(v_ls);
         free(w_ls);
         free(h);
         free(g);
         free(minimizer);
         free(norm_r_ls);
         free(scalar1_ls);
         free(scalar2_ls);
         free(temp1);
         free(dot_temp_ls);
         free(h_min);
         free(g_min);
         free(jpvt);
         free(eff_rank_h_min);
         free(reqs_ls);
         free(dot_temp);
         free(dot_temp_1);
         free(dot_temp_2);
         free(norm_temp);
         free(test);

         //printf("i am in if\n");
     
         return;
       }
       else
       {
         // x0=xm;  
         if(id ==0)
         memcpy(&X0[shift_init_Mi[id]-1], &temp2[shift_init_Mi[id]-1],(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_high)*sizeof(double complex)));
         else if(id != 0 && id != num_procs-1)
         memcpy(&X0[shift_init_Mi[id]-1-overlap_low], &temp2[shift_init_Mi[id]-1-overlap_low],(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_tol)*sizeof(double complex)));
         else
         memcpy(&X0[shift_init_Mi[id]-1-overlap_low], &temp2[shift_init_Mi[id]-1-overlap_low],(((shift_end_Mi[id]-shift_init_Mi[id])+1+overlap_low)*sizeof(double complex)));
     
         // restart=restart+1;
         restart_count = restart_count + 1;
         // if abs(g(m+1,1)) <= normmin
         if(sqrt(creal(g[restart])*creal(g[restart])+cimag(g[restart])*cimag(g[restart]))<= norm_min)
         {
           // xmin = xm;
           memcpy(Xmin, temp2,(sizeof(double complex)*n_ham));
           norm_min = sqrt(creal(g[restart])*creal(g[restart])+cimag(g[restart])*cimag(g[restart]));
         }
       }
     
       } // end of while

       
     // x = xmin;
     memcpy(t, Xmin,(sizeof(double complex)*n_ham));
     
     *ls_counter =  *ls_counter + restart_count * restart;
     
      free(X0);
      free(Xmin);
      free(dot);
      free(temp);
      free(temp2);
      free(rq);
      free(beta);
      free(v_ls);
      free(w_ls);
      free(h);
      free(g);
      free(minimizer);
      free(norm_r_ls);
      free(scalar1_ls);
      free(scalar2_ls);
      free(temp1);
      free(dot_temp_ls);
      free(h_min);
      free(g_min);
      free(jpvt);
      free(eff_rank_h_min);
      free(reqs_ls);
      free(dot_temp);
      free(dot_temp_1);
      free(dot_temp_2);
      free(norm_temp);
      free(test);

     return;
     }     
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
     void c_zdotu(int size, double complex * array_a, double complex * array_b, double complex * output)
        {
          int i;
          double complex temp = 0.000000000 + 0.000000000 * I;
          for( i =0; i < size; i++)
          {
            temp = temp + (((creal(array_a[i]) * creal(array_b[i])) + cimag(array_a[i]) * cimag(array_b[i])) + (((cimag(array_a[i]) * creal(array_b[i])) - creal(array_a[i]) * cimag(array_b[i])) * I));
          }

         *output = temp;
        }
    
    
    
    
    
    
    
     void vct1_add_vct2_asg_vct3(double complex *vct1, double complex * vct2, double complex * vct3, int n_ham)
    {
      int ID;
      #pragma omp parallel for shared(vct1, vct2, vct3, n_ham) private(ID)
      for(ID = 0; ID < n_ham; ID++)
      {
      vct3[ID] = (creal(vct1[ID]) + creal(vct2[ID])) + ((cimag(vct1[ID]) + cimag(vct2[ID])) * I);
      }
    }
    
    
    
    
     void vct1_sub_vct2_asg_vct3(double complex *vct1, double complex * vct2, double complex * vct3, int n_ham)
    {
      int ID;
      #pragma omp parallel for shared(vct1, vct2, vct3, n_ham) private(ID)
      for(ID = 0; ID < n_ham; ID++)
      {
      vct3[ID] = (creal(vct1[ID]) - creal(vct2[ID])) + ((cimag(vct1[ID]) - cimag(vct2[ID])) * I);
      }
    }
    
    
    
     void vct1_neg_asg_vct2(double complex *vct1, double complex * vct2, int n_ham)
    {
      int ID;
      #pragma omp parallel for shared(vct1, vct2, n_ham) private(ID)
      for(ID = 0; ID < n_ham; ID++)
      {
        vct1[ID] =  -creal(vct2[ID]) + (-cimag(vct2[ID]) * I);
      }
    }
    
    
    
     void vct_div_slr_minus_scl_vct(double complex * scr1, double complex * scr2, double complex * des, int n_ham, double Slr1, double Slr2)
    {
    int ID;
    #pragma omp parallel for shared(des, scr1, scr2, Slr1, Slr2, n_ham) private(ID)
    for(ID = 0; ID < n_ham; ID++)
    {
      des[ID] = ((creal(scr1[ID])/Slr1) - (Slr2 * creal(scr2[ID]))) + (((cimag(scr1[ID])/Slr1) - (Slr2 * cimag(scr2[ID]))) * I) ;
    }
    }
    
    
    
     void vct_div_slr(double complex * scr, double complex * des, int n_ham, double Slr)
    {
    int ID;
    #pragma omp parallel for shared(des, scr, Slr, n_ham) private(ID)
    for(ID = 0; ID < n_ham; ID++)
    {
      des[ID] = (creal(scr[ID]) / Slr) + ((cimag(scr[ID]) / Slr) * I);
    }
    }
    



    
     void vct_1_div_asg_to_vct_2(double complex * vct_1, double complex * vct_2, int n_ham, double beta)
    {
      int ID;
      #pragma omp parallel for shared(vct_1, vct_2, beta, n_ham) private(ID)
      for(ID = 0; ID < n_ham; ID++)
      {
        vct_2[ID] =  (creal(vct_1[ID]) / beta) + ((cimag(vct_1[ID]) / beta) * I);
      }
    }
    
    
    
    
    
    
    
    
     void vct1_sub_mul_vct(double complex * vct1, double complex * vct2, int n_ham, double complex * dot)
    {
      int  ID;
      //double complex * temp_array;
      //temp_array = (double complex *) malloc((n_ham)*sizeof(double complex));
      #pragma omp parallel for shared(vct1, vct2, dot, n_ham) private(ID)
      for(ID = 0; ID < n_ham; ID++)
      {
        vct1[ID] = (creal(vct1[ID]) - ((creal(dot[0]) * creal(vct2[ID])) - (cimag(dot[0]) * cimag(vct2[ID])))) + ((cimag(vct1[ID]) - ((creal(dot[0]) * cimag(vct2[ID])) + (cimag(dot[0]) * creal(vct2[ID])))) * I) ;
        //vct1[ID] = temp_array[ID];
      }
      //free(temp_array);
    }
    
    
    
    
    
    
    
     void vct_sub_scl_mul_vct(double complex *vct1, double complex * vct2, double scalar, int n_ham)
    {
      int ID;
      //double complex * temp_array;
      //temp_array = (double complex *) malloc((n_ham)*sizeof(double complex));
      #pragma omp parallel for shared(vct1, vct2, scalar, n_ham) private(ID)
      for(ID = 0; ID < n_ham; ID++)
      {
        vct1[ID] =  (creal(vct1[ID]) - (scalar * creal(vct2[ID]))) + ((cimag(vct1[ID]) - (scalar * cimag(vct2[ID]))) * I);
      }
      //free(temp_array);
    }
    
    
    
    
    
    void shift_A(float * val, int * row, int * col, int n_ham, float shift, int offset)
    {
     int ID, j;
     #pragma omp parallel for shared(row, n_ham, offset, shift, col, val) private(ID, j)
     for(ID = 0; ID < n_ham; ID++)
      {
         for(j=row[ID+offset]-1; j <= row[ID+1+offset]-1; j++)
         {
          if(col[j] == ID+1+offset)
          {
            val[j] = val[j] + shift; 
            break;
          }
         }
      }
    }
    
    
    
    
     void cpy_vct_1_to_vct_2(double complex * vct_scr, double complex * vct_des, int n_ham)
    {
      int ID;
      #pragma omp parallel for shared(vct_des, vct_scr, n_ham) private(ID)
      for(ID = 0; ID < n_ham; ID++)
      {
        vct_des[ID] = (creal(vct_scr[ID])) + ((cimag(vct_scr[ID])) * I);
      }
    }
    
    
//      void mv(int n_ham, int col, double complex * xVal_kr, double complex * y_kr, double complex * Finalans_kr){
//                             //kernel func varaibles
//                             int ID, row_start, row_end, jj;
//                             double dot_x, dot_img;
//                             for(ID = 0; ID < n_ham; n_ham++){
//                                             dot_x=0.0;
//                                             dot_img=0.0; 
//                                             //row_start = csrRowPtr_kr[row]-1;
//                                             //row_end = csrRowPtr_kr[row+1]-1;
//                                             //(x + yi)(u + vi) = (xu ? yv) + (xv + yu)i. 
//                                             for(jj = 0; jj < col; jj++ ){
//                                                       dot_x += ((xVal_kr[(jj*n_ham)+ID].x * y_kr[jj].x) - (xVal_kr[(jj*n_ham)+ID].y * y_kr[jj].y)); 
//                                                       dot_img  += ((xVal_kr[(jj*n_ham)+ID].x * y_kr[jj].y) + (xVal_kr[(jj*n_ham)+ID].y * y_kr[jj].x)); 
//                                                     }
//                     
//                                             Finalans_kr[ID].x = dot_x;  
//                                             Finalans_kr[ID].y = dot_img;  
//                                             
//                                     }
//                     }
    
    
    
    
    
     void spmv_csr_hybrid(int num_rows, const int* rowPtrs, const int* colIdxs, const float* values, 
                                          const double complex* x, double complex* y, const int offset)
    {
      int i;
      double complex rowSum;
      int j;
      int row_start, row_end;
      //clock_t begin, end;
      //begin = clock();
      #pragma omp parallel for shared(x, values, colIdxs, y, rowPtrs, num_rows) private(j, i, row_start, row_end) reduction(+:rowSum)
      for(i = 0; i < num_rows; i++)
      {
        rowSum = 0.0000000 + 0.00000000 *I;
        row_start = rowPtrs[offset+i]-1;
        row_end = rowPtrs[offset+i+1]-1;
        for (j=row_start; j<row_end; j++)
        {
            rowSum += (values[j] * creal(x[colIdxs[j]-1])) + ((values[j] * cimag(x[colIdxs[j]-1])) * I);
        }
        y[offset+i] = rowSum;
      }
      //end = clock();
      //printf("Elapsed: %f seconds\n", (double)(end - begin) / CLOCKS_PER_SEC);
    }
    
    
     void vct_pls_scl_mul_vct(double complex *vct1, double complex * vct2, double complex scalar, int n_ham)
    {
      int ID;
      //double complex * temp_array;
      //temp_array = (double complex *) malloc((n_ham)*sizeof(double complex));
      //clock_t begin, end;
      //begin = clock();

      #pragma omp parallel for shared(vct1, vct2, scalar, n_ham) private(ID)
      for(ID = 0; ID < n_ham; ID++)
      {
        vct1[ID] =  (creal(vct1[ID]) + creal(scalar) * creal(vct2[ID]) - cimag(scalar) * cimag(vct2[ID])) + ((cimag(vct1[ID]) + creal(scalar) * cimag(vct2[ID]) + cimag(scalar) * creal(vct2[ID])) * I);
        //vct1[ID] = temp_array[ID];
      }
      //end = clock();
      //printf("Elapsed: %f seconds\n", (double)(end - begin) / CLOCKS_PER_SEC);
      //free(temp_array);
    }
    
    
    
    //  void Vector_Dot_Product( const double complex *V1 , const double complex *V2 , double complex *V3, int n_ham, int offset  )
    // {
    //  __shared__ double complex chache[BLOCK_SIZE] ;
    //    
    //  double complex temp ;
    // 
    //  unsigned int tid = blockDim.x * blockIdx.x + threadIdx.x ;
    // 
    //  unsigned int chacheindex = threadIdx.x ;
    // 
    //  temp.x = 0.0; temp.y = 0.0;
    // 
    //  while ( tid < n_ham )
    //  {
    //       temp.x += (V1[tid].x * V2[tid].x) - (V1[tid].y * V2[tid].y) ;
    //       temp.y += (V1[tid].x * V2[tid].y) + (V1[tid].y * V2[tid].x) ;
    //       
    //       tid += blockDim.x * gridDim.x ;
    //  }
    // 
    //   chache[chacheindex] = temp ;
    // 
    //  __syncthreads();
    // 
    //  int i  = blockDim.x / 2 ;
    // 
    // while ( i!=0 )
    //  {
    // 
    //   if ( chacheindex < i )
    //   {
    //          chache[chacheindex].x += chache[chacheindex + i].x ;
    //          chache[chacheindex].y += chache[chacheindex + i].y ;
    //   }
    //  
    // __syncthreads();
    // 
    //    i/=2 ;
    //  }
    // 
    //   if ( chacheindex == 0 )
    //          *V3 = chache [0] ;
    // 
    // 
    // }
    
    // __device__ double cuda_atomicAdd(double *address, double val);
    // 
    //  void Vector_Dot_Product(const double complex  *a, const double complex  *b, double complex  *dot_res, int n_ham, int offset)
    // {
    //     __shared__ double complex cache[BLOCK_SIZE]; //thread shared memory
    //     int global_tid=threadIdx.x + blockIdx.x * blockDim.x;
    //     int i=0,cacheIndex=0;
    //     double complex temp = make_double complex(0.0,0.0);
    //     cacheIndex = threadIdx.x;
    //     while(global_tid < n_ham){
    //         temp.x += (a[global_tid].x * b[global_tid].x) - (a[global_tid].y * b[global_tid].y);
    //         temp.y += (a[global_tid].x * b[global_tid].y) + (a[global_tid].y * b[global_tid].x);
    //         global_tid += blockDim.x * gridDim.x;
    //     }
    //     cache[cacheIndex] = temp;
    //     __syncthreads();
    //     for (i=blockDim.x/2; i>0; i>>=1) {
    //         if (threadIdx.x < i) {
    //             cache[threadIdx.x].x += cache[threadIdx.x + i].x;
    //             cache[threadIdx.x].y += cache[threadIdx.x + i].y;
    //         }
    //         __syncthreads();
    //     }
    //     __syncthreads();
    // 
    //     //dot_res[0] = make_double complex(0.0,0.0);
    //     if (cacheIndex==0) {
    //     double result_x = cuda_atomicAdd(&dot_res[0].x,cache[0].x);
    //     double result_y = cuda_atomicAdd(&dot_res[0].y,cache[0].y);
    //     //dot_res[0].x = dot_res[0].x+cache[0].x;
    //     //dot_res[0].y = dot_res[0].y+cache[0].y;
    //     }
    //     __syncthreads();
    // }
    // 
    // 
    // __device__ double cuda_atomicAdd(double *address, double val)
    // {
    //     double assumed,old=*address;
    //     do {
    //         assumed=old;
    //         old= __longlong_as_double(atomicCAS((unsigned long long int*)address,
    //                     __double_as_longlong(assumed),
    //                     __double_as_longlong(val+assumed)));
    //     }while (assumed!=old);
    // 
    //     return old;
    // }
