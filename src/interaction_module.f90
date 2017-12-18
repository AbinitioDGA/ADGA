module interaction_module
  use parameters_module
  use hdf5_module
  use aux
  implicit none

contains

!========================================================================================================
! read the v(q) from a hdf5 file with specific format
subroutine read_vq(iq, v,er,erstr)
  use aux
  implicit none
  integer, intent(in) :: iq
  complex(kind=8), intent(out) :: v(ndim2,ndim2)
  integer,intent(out) :: er
  !character(len=*),intent(out) :: erstr
  complex(kind=8) :: vq(ndim,ndim,ndim,ndim)
  integer(hid_t) :: vq_file_id, grp_id, iq_id, iq_space_id
  integer :: err, ind, i, j, k, l, i1, i2
  integer :: nmembers, imembers, itype!,er
  character(len=20) :: name_buffer
  character(len=*) :: erstr
  integer(hsize_t), dimension(2) :: iq_dims, iq_maxdims
  integer(hsize_t), dimension(1) :: vq_dims
  double precision :: vq_tmp_r(nqp), vq_tmp_i(nqp)
  complex(kind=8),parameter :: ci = (0d0,1d0)

  er = 0
  erstr = ''

  call h5fopen_f(filename_vq, h5f_acc_rdonly_f, vq_file_id, err)

  call h5dopen_f(vq_file_id, ".axes/Q-points", iq_id, err)
  call h5dget_space_f(iq_id, iq_space_id, err)
  call h5sget_simple_extent_dims_f(iq_space_id, iq_dims, iq_maxdims, err)
  if(iq_dims(2) .ne. nqp) then
     er = 1
     write(erstr,*) 'Inconsistent number of q-points in V^q!', iq_dims(2),'/',nqp
     return
  endif

  vq = 0.d0

  call h5gn_members_f(vq_file_id, "/", nmembers, err)
  do imembers = 1,nmembers - 1
     call h5gget_obj_info_idx_f(vq_file_id, "/", imembers, name_buffer, itype, err)

     read(name_buffer,'(I5.5)') ind
     call index2component_band(ndim,ind,i,j,k,l)
     call h5dopen_f(vq_file_id, name_buffer, grp_id, err)
     call h5dread_f(grp_id, type_r_id, vq_tmp_r, vq_dims, err)
     call h5dread_f(grp_id, type_i_id, vq_tmp_i, vq_dims, err)

     vq(i,j,k,l) = vq_tmp_r(iq)+ci*vq_tmp_i(iq)

     call h5dclose_f(grp_id, err)
  enddo
  call h5fclose_f(vq_file_id, err)

  v = 0.d0
  i2 = 0
  do l=1,ndim
     do j=1,ndim
        i2 = i2+1
        i1 = 0
        do i=1,ndim
           do k=1,ndim
              i1 = i1+1
              v(i1,i2) = vq(i,j,k,l)
           enddo
        enddo
     enddo
  enddo

end subroutine read_vq
!========================================================================================================

!========================================================================================================
! read a umatrix file with all possible leg combinations
subroutine read_u(u, u_tilde)
  implicit none
  real(kind=8) :: u_tmp(ndim,ndim,ndim,ndim), u_tilde_tmp(ndim,ndim,ndim,ndim)
  real(kind=8) :: u_value
  complex(kind=8), intent(out) :: u(ndim**2, ndim**2), u_tilde(ndim**2, ndim**2)
  integer :: n,i,j,k,l,i1,i2


  open(21,file=filename_umatrix,status='old')
  read(21,*)
  do n=1,ndim**4
     read(21,*) i, j, k, l, u_value
     u_tmp(i,j,k,l) = u_value
     u_tilde_tmp(i,j,l,k) = u_value
  enddo
  close(21)

  !go into compound index:
  u = 0.d0
  u_tilde = 0.d0
  i2 = 0
  do l=1,ndim
     do j=1,ndim
        i2 = i2+1
        i1 = 0
        do i=1,ndim
           do k=1,ndim
              i1 = i1+1
              u(i1,i2) = u_tmp(i,j,k,l)
              u_tilde(i1,i2) = u_tilde_tmp(i,j,k,l)
           enddo
        enddo
     enddo
  enddo

end subroutine read_u
!========================================================================================================

!========================================================================================================
! create the umatrix array from the given config parameters
subroutine create_u(u, u_tilde)
  implicit none
  real(kind=8) :: u_tmp(ndim,ndim,ndim,ndim), u_tilde_tmp(ndim,ndim,ndim,ndim)
  real(kind=8) :: u_value
  complex(kind=8), intent(out) :: u(ndim**2, ndim**2), u_tilde(ndim**2, ndim**2)
  integer :: n,i,j,k,l,i1,i2,ineq

  allocate(Umat(ndim,ndim,ndim,ndim))
  Umat=0.d0

  ineq = 0
  do i=1,ndim
  do j=1,ndim
  do k=1,ndim
  do l=1,ndim

    ineq=index2ineq(nineq,ndims,i,j,k,l) 
    if (ineq .eq. 0) cycle ! not on the same impurity
  
    ! DD - VALUES
    if (index2cor(nineq,ndims,i,j,k,l)) then
      if (i .eq. j .and. k .eq. l)  then
        if (i .eq. k) then
          Umat(i,j,k,l) = Udd(ineq) ! onsite density-density -- 1 1 1 1
          cycle
        else if (interaction_mode(ineq) .eq. 1) then
          Umat(i,j,k,l) = Jdd(ineq) ! kanamori double hopping -- 1 1 2 2
          cycle
        endif
      endif

      if (i .eq. k .and. j .eq. l) then
        Umat(i,j,k,l) = Vdd(ineq) ! screened density density -- 1 2 1 2
        cycle
      endif
      if (i .eq. l .and. j .eq. k .and. interaction_mode(ineq) .eq. 1) then
        Umat(i,j,k,l) = Jdd(ineq) ! kanamori spin flip -- 1 2 2 1
        cycle
      endif

    ! PP - VALUES
    else if(index2uncor(nineq,ndims,i,j,k,l)) then
      if (i .eq. j .and. k .eq. l)  then
        if (i .eq. k) then
          Umat(i,j,k,l) = Upp(ineq) ! onsite density-density -- 1 1 1 1
          cycle
        else if (interaction_mode(ineq) .eq. 1) then
          Umat(i,j,k,l) = Jpp(ineq) ! kanamori double hopping -- 1 1 2 2
          cycle
        endif
      endif

      if (i .eq. k .and. j .eq. l) then
        Umat(i,j,k,l) = Vpp(ineq) ! screened density density -- 1 2 1 2
        cycle
      endif
      if (i .eq. l .and. j .eq. k .and. interaction_mode(ineq) .eq. 1) then
        Umat(i,j,k,l) = Jpp(ineq) ! kanamori spin flip -- 1 2 2 1
        cycle
      endif

    ! DP - VALUES
    else
      if (i .eq. j .and. k .eq. l)  then
        if (i .eq. k) then
          Umat(i,j,k,l) = Udp(ineq) ! onsite density-density -- 1 1 1 1
          ! this band combination can not exist if there should be 2 different kind of bands
          cycle
        else if (interaction_mode(ineq) .eq. 1) then
          Umat(i,j,k,l) = Jdp(ineq) ! kanamori double hopping -- 1 1 2 2
          cycle
        endif
      endif

      if (i .eq. k .and. j .eq. l) then
        Umat(i,j,k,l) = Vdp(ineq) ! screened density density -- 1 2 1 2
        cycle
      endif
      if (i .eq. l .and. j .eq. k .and. interaction_mode(ineq) .eq. 1) then
        Umat(i,j,k,l) = Jdp(ineq) ! kanamori spin flip -- 1 2 2 1
        cycle
      endif
    endif

  enddo
  enddo
  enddo
  enddo

  deallocate(interaction_mode)
  deallocate(Udd,Vdd,Jdd,Upp,Vpp,Jpp,Udp,Vdp,Jdp)

  !go into compound index:
  u = 0.d0
  u_tilde = 0.d0
  i2 = 0
  do l=1,ndim
     do j=1,ndim
        i2 = i2+1
        i1 = 0
        do i=1,ndim
           do k=1,ndim
              i1 = i1+1
              u(i1,i2) = Umat(i,j,k,l)
              u_tilde(i1,i2) = Umat(i,j,l,k)
           enddo
        enddo
     enddo
  enddo
end subroutine create_u

end module interaction_module
