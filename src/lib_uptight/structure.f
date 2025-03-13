    basis_loop:DO i_a = 1, n_basis

       IF ( verbose.gt.0 .AND. MOD( i_a, 1000 ) .EQ. 0.0 ) THEN 
       END IF
       !IF ( err .NE. 0 ) CALL alloc_error( 'TB_ham', 'sparse_ham', 'a_ref' )
       !IF ( err .NE. 0 ) CALL alloc_error( 'TB_ham', 'sparse_ham', 'onsite' )
       IF ( d_onsite_shift_flag ) THEN
       END IF
       IF ( potential_flag ) THEN
       ENDIF
       IF ( relat ) THEN
       END IF
       DO i_cpl = 1, n_neig
       END DO

       ! -------------------------------------------------------------------------
       nn_loop:DO i_n = 1, SIZE( current_near%ind ) 
          IF ( i_b .NE. i_b_old ) THEN
             ! Increment coupling bloc counter if required
             IF ( ASSOCIATED( cpl( i_cpl )%mat ) ) THEN
                IF ( i_cpl .LT. n_neig ) i_cpl = i_cpl + 1
             END IF
          END IF
          IF ( current_near%order( i_n ) .NE. 0 ) THEN
             IF ( basis%mat( i_a ) .EQ. basis%mat( i_b ) ) THEN  
                IF (i_parent .EQ. 0) STOP 'ERROR: pair not found'
             ELSE
                IF (i_parent .EQ. 0) STOP 'ERROR: pair not found'
             END IF
             IF (err.NE.0) CALL alloc_error('TB_ham','sparse_ham','ref_pair')
          END IF  ! End of non-onsite test

          ! IF NOT DANGLING BOND:
          IF ( .NOT. is_dg_bond( i_n, basis, current_near ) ) THEN
             IF ( .NOT. ASSOCIATED( cpl( i_cpl )%mat ) ) THEN
                IF (err.NE.0) CALL alloc_error('TB_ham','sparse_ham','cpl%mat')
             END IF
             ! ONSITE
             IF ( current_near%order( i_n ) .EQ. 0 ) THEN
                IF (optmat) THEN
                ELSE
                   DO a_st = 1, n_st( i_a )
                   END DO
                   IF ( .not.syst_rotation .and. & 
                      IF (cf(1).ne.cf(2)) then
                      END IF
                   END IF
                END IF
             ! OFF-SITE ELEMENTS
             ELSE
                IF (optmat) THEN
                ELSE
                END IF
                DO a_st = 1, n_st( i_a )
                   DO b_st = 1, n_st( i_b )
                   END DO
                END DO
             END IF
          ! DANGLING BOND:
          ELSE        
             IF (err.NE.0) CALL alloc_error('TB_ham','sparse_ham','hydro%mat' )
             IF (optmat) THEN
             ELSE
             END IF
             DO a_st = 1, n_st( i_a )
                DO b_st = 1, n_st( i_b )
                END DO
             END DO
          END IF    !   if (.not. is_dg_bond())

          IF ( ASSOCIATED( ref_pair ) ) THEN
             IF ( err .NE. 0 ) &
          END IF
       END DO nn_loop    ! IF i_n = 1, ind
       ! END OF NN LOOP
       
       ! ----------------------------------------------------------------------
       ! STORE VALUES ON THE MATRIX:
       spin_loop:DO i_spin_a = 1, n_spin
          IF (relat) THEN
             IF (optmat) THEN            
             ELSE
             END IF
          END IF
          state_loop:DO a_st = 1, n_st( i_a )
             col_loop:DO i_cpl = 1, n_neig

                ! STORE ONSITE
                IF ( i_b .EQ. i_a ) THEN

                   DO i_spin_b = 1, n_spin
                      DO b_st = 1, n_st( i_b )
                         SELECT CASE ( sparse_format )
                         CASE('U')
                            IF ( n_col .LT. n_row ) CYCLE
                         CASE('L')
                            IF ( n_col .GT. n_row ) CYCLE
                         CASE('F')
                         END SELECT
                         IF ( i_spin_a .EQ. i_spin_b ) THEN
                            IF ( relat ) value = value &
                         ELSE
                         END IF
                         IF (i_spin_a.eq.i_spin_b .and. a_st.eq.b_st) THEN
                         ELSE
                         ENDIF
                         IF ( value_test ) THEN
                         ELSE
                         END IF
                      END DO ! loop over b_st
                   END DO  ! loop over i_spin_b

                   ! Align n_col to end of atom block ( include n_st(i_a) ) 
                   DO i_dg_a = 1, n_dg(i_a)
                      SELECT CASE ( sparse_format )
                      CASE('U')
                         IF ( n_col .LT. n_row ) CYCLE
                      CASE('L')
                         IF ( n_col .GT. n_row ) CYCLE
                      CASE('F')
                      END SELECT
                      IF (optmat) then
                      ELSE
                      END IF
                      IF (i_spin_a.eq.i_spin_b) THEN
                      ELSE
                      ENDIF
                      IF ( value_test ) THEN
                      ELSE
                      END IF
                   END DO
                ! STORE OFF-SITE
                ELSE  !( i_a .NE. i_b )
                   DO b_st = 1, n_st( i_b )
                      SELECT CASE ( sparse_format )
                      CASE('U')
                         IF ( n_col .LT. n_row ) CYCLE
                      CASE('L')
                         IF ( n_col .GT. n_row ) CYCLE
                      CASE('F')
                      END SELECT
                      IF (.not. ASSOCIATED(cpl(i_cpl)%mat)) THEN
                      END IF
                      IF ( value_test ) THEN
                      ELSE
                      END IF
                   END DO
                END IF

             END DO col_loop  ! i_cpl = 1, n_neig
          END DO state_loop  ! a_st = 1, n_st(i_a)
       END DO spin_loop  ! i_spin_a
       dg_loop:DO i_spin_a = 1, n_spin
          DO i_dg_a = 1, n_dg(i_a)
             DO a_st = 1, n_st( i_a )
                SELECT CASE ( sparse_format )
                CASE('U')
                  IF ( n_col .LT. n_row ) CYCLE
                CASE('L')
                  IF ( n_col .GT. n_row ) CYCLE
                CASE('F')
                END SELECT
                IF (optmat) then
                ELSE
                END IF
                IF ( value_test ) THEN
                ELSE
                END IF
             END DO
             IF (optmat) THEN
             ELSE
             END IF
          END DO ! end loop on i_dg_a
       END DO dg_loop ! end loop on i_spin_a
       DO i_cpl = 1, n_neig
          IF ( err .NE. 0 ) &
       END DO
       DO i_dg_a = 1, n_dg(i_a)
          IF ( err .NE. 0 ) &
       END DO
       !IF ( err .NE. 0 ) &

    END DO basis_loop  ! end loop on i_a
    ! END BASIS LOOP
    ! +=============================================

    IF ( i_ind .NE. 0 ) THEN
    END IF
    IF ( i_val .NE. 0 ) THEN
    END IF
    SELECT CASE (TRIM(sparse_format))
    END SELECT
    IF ( ioutput_flag .GE. 2 ) THEN
    END IF


  END SUBROUTINE sparse_ham


    IF ( alloy_random .OR. TRIM(type).EQ.'simple' &
       ! --> random alloy : do not average
       IF ( scale ) then             
       ELSE
       END IF 
       DO i_c = 1, n_cpl          
       ENDDO
    ELSE
       SELECT CASE( TRIM(type) )
          IF (scale) THEN
          ELSE
          ENDIF
          DO i_c = 1, n_cpl
          END DO
          ! to do : implement alloying for quaternary pairs
       END SELECT
    END IF
  END SUBROUTINE harrison_scaling
  ! => no_d_flag - logical : true if no d-states are included.
  END SUBROUTINE koster_slater
    if (poldir.eq.1) THEN
    else if (poldir.eq.2) THEN
    elseif (poldir.eq.3) THEN
    end if
  END SUBROUTINE intratomic_optics
    IF ( i_spin_a .EQ. 1 ) THEN
    ELSE IF ( i_spin_a .EQ. 2 ) THEN
    END IF
  END SUBROUTINE spin_sparse
    do k=1,5
       do l=1,5
          IF ( ABS(diag(k,l)).LT.1.0d-14 ) diag(k,l)=(0.d0,0.d0)
       end do
    end do
  END SUBROUTINE CRYSTAL_FIELD
    DO i_state = 1, SIZE( ion(i_ion)%energy )
       SELECT CASE ( ion(i_ion)%state( i_state ) )
       END SELECT
    END DO
  END SUBROUTINE d_onsite_strain
  END SUBROUTINE hydrogen_coupling
!!$    DO i_ind = 1, SIZE( sparse_ind )
!!$       IF ( sparse_ind( i_ind ) .LT. 0 ) THEN
!!$       ELSE
!!$       END IF
!!$       IF ( i_row .EQ. i_col ) THEN
!!$       ELSE
!!$       END IF
!!$    END DO
!!$  END SUBROUTINE TB_sparse_davidson
    IF ( i_ind .EQ. max_ind_D ) THEN 
    END IF
  END SUBROUTINE store_ind_file
  END SUBROUTINE store_ind_mem
    IF ( i_val .EQ. max_ind_D ) THEN 
    END IF
  END SUBROUTINE store_val_file
  END SUBROUTINE store_val_mem
    DO i_state = 1, SIZE( mat_data( i_mat )%ion( i_ion )%energy )
       SELECT CASE ( mat_data( i_mat )%ion( i_ion )%state( i_state ) )
       END SELECT
    END DO
  END FUNCTION d_orbital_present
    SELECT CASE(ham%sparse_fmt)
       do row=1,n_ham
          do p = Mij(row), Mij(row+1) - 1 
             end if
          enddo
       enddo
    END SELECT
  END SUBROUTINE hermitianize
    SELECT CASE(ham%sparse_fmt)
       do row=1,n_ham
          do p = Mij(row), Mij(row+1) - 1 !row+1,n_ham
             end if
          enddo
       enddo
    END SELECT
  END SUBROUTINE check_if_hermitian
    SELECT CASE(ham%sparse_fmt)
       do row=1,n_ham
          do p = Mij(row), Mij(row+1) - 1 !row+1,n_ham
             end if
          enddo
       enddo
    END SELECT
  END SUBROUTINE check_if_antihermitian
    DO i_a = 1, n_basis
       DO i_n = 1, SIZE( current_near%ind ) 
          IF ( is_dg_bond( i_n, basis, current_near ) ) &
       END DO
    END DO
  END SUBROUTINE get_ion_block_size
