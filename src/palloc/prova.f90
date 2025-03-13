program prova
        
 use allocation

 type(mem_log) :: mlog
 integer, dimension(:), allocatable :: int_alloc
 character(5), dimension(:), pointer :: int_p
 !integer, dimension(:,:), pointer :: int_p2

 NULLIFY(int_p)
 
 call openMemLog(6,mlog)
 call resetMemLog(mlog)

 call log_allocatep(int_p,10,mlog)

 call writeMemLog(mlog)

 call log_deallocatep(int_p,mlog)

 !call log_allocatep(int_p2,20,20,mlog)

 !call log_deallocatep(int_p2,mlog)


 call writeMemLog(mlog)
 call closeMemLog(mlog)

end program prova
