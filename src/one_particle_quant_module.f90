module one_particle_quant_module
  use lapack_module
  use parameters_module
  use aux
  use mpi_org, only: mpi_wrank,master
  implicit none

  contains

! calculation of greens functions
!=========================================================================== 
subroutine get_giw()
  ! calculation of dmft impurity greens function with orbital diagonal self energy
  ! we do the inversion and then only use the diagonal elements
  implicit none
  integer :: ik, iw, i
  complex(kind=8) :: g(ndim,ndim), g2(ndim,ndim)
  complex(kind=8),parameter :: ci = (0d0,1d0)

  giw = 0.d0
  do ik=1,nkp
     g(:,:) = -hk(:,:,ik)
     do iw=0,iwmax-1 !use symmetry of giw(-w)=giw^*(w)
        do i=1,ndim
           g(i,i) = ci*iw_data(iw)+mu-hk(i,i,ik)-dc(1,i)
        enddo
        do i=1,ndim
           g(i,i) = g(i,i)-siw(iw,i) !no spin dependence in single particle Greens function
        enddo
        g2 = g(:,:)
        call inverse_matrix(g2)
        do i=1,ndim
           giw(iw,i) = giw(iw,i)+g2(i,i)
        enddo
     enddo
  enddo

  do iw=0,iwmax-1
     do i=1,ndim
        giw(-iw-1,i) = conjg(giw(iw,i)) ! It is diagonal, so no transpose is needed
     enddo
  enddo

  giw = giw/dble(nkp)

  return
end subroutine get_giw


subroutine get_gkiw(ikq, iwf, iwb, gkiw)
  ! calculation of lattice greens function with dmft selfenergy
  ! at a specific k-point ikq and a fermionic frequency iwf-iwb
  implicit none
  integer :: i
  integer :: iwf, iwb, ikq
  complex(kind=8), intent(out) :: gkiw(ndim,ndim)
  complex(kind=8),parameter :: ci = (0d0,1d0)


  gkiw(:,:) = -hk(:,:,ikq)
  do i=1,ndim
     gkiw(i,i) = ci*iw_data(iwf-iwb)+mu-hk(i,i,ikq)-dc(1,i)
  enddo
  do i=1,ndim
  gkiw(i,i) = gkiw(i,i)-siw(iwf-iwb,i)
  enddo
  call inverse_matrix(gkiw)


end subroutine get_gkiw

subroutine get_gkiw_dga(ik, iwf, skiw, gkiw)
  ! calculation of lattice greens function with any selfenergy
  ! skiw(ndim,ndim) at a specific k-point ik and fermionic frequency iwf
  implicit none
  integer :: i
  integer, intent(in) :: iwf, ik
  complex(kind=8), intent(in) :: skiw(ndim,ndim)
  complex(kind=8), intent(out) :: gkiw(ndim,ndim)
  complex(kind=8),parameter :: ci = (0d0,1d0)


  gkiw(:,:) = -hk(:,:,ik)-skiw(:,:)
  do i=1,ndim
     gkiw(i,i) = ci*iw_data(iwf)+mu-hk(i,i,ik)-dc(1,i)-skiw(i,i)
  enddo
  call inverse_matrix(gkiw)

end subroutine get_gkiw_dga


subroutine get_sigma_g_loc(sigma_sum, sigma_loc, gloc)
  ! calculation of the k-summed dga self-energy
  ! calculation of the k-summed greens function with the dga selfenergy
  ! inside the vertex box and the dmft selfenergy outisde of it
  implicit none
  complex(kind=8),intent(in) :: sigma_sum(ndim,ndim,-iwfmax_small:iwfmax_small-1,nkp)
  complex(kind=8), intent(out) :: sigma_loc(ndim, ndim,-iwfmax_small:iwfmax_small-1)
  complex(kind=8), intent (out) :: gloc(-iwmax:iwmax-1,ndim,ndim)

  integer :: ik, iw, i
  complex(kind=8) :: g(ndim,ndim), g2(ndim,ndim)
  complex(kind=8),parameter :: ci = (0d0,1d0)

  ! k-summed dga self-energy
  sigma_loc = 0.d0
  do ik=1,nkp
     sigma_loc(:,:,:) = sigma_loc(:,:,:)+sigma_sum(:,:,:,ik)
  enddo
  sigma_loc = sigma_loc/dble(nkp)

  ! k-summed dga (+dmft) greens-function
  gloc = 0.d0
  do ik=1,nkp
     ! Within the small box, use sigma_sum (not orbital diagonal)
     do iw=0,iwfmax_small-1
        g(:,:) = -hk(:,:,ik) - sigma_sum(:,:,iw,ik)
        do i=1,ndim
           g(i,i) = ci*iw_data(iw)+mu-hk(i,i,ik)-dc(1,i)-sigma_sum(i,i,iw,ik)
        enddo
        g2 = g(:,:)
        call inverse_matrix(g2)
        gloc(iw,:,:) = gloc(iw,:,:)+g2(:,:)
     enddo

     ! Outside the small box, use siw (from QMC)
     do iw=iwfmax_small,iwmax-1
        g(:,:) = -hk(:,:,ik)
        do i=1,ndim
           g(i,i) = ci*iw_data(iw)+mu-hk(i,i,ik)-dc(1,i)-siw(iw,i)
        enddo
        g2 = g(:,:)
        call inverse_matrix(g2)
        gloc(iw,:,:) = gloc(iw,:,:)+g2(:,:)
     enddo

  enddo

  ! we use G^dagger (nu) = G(-nu)
  do iw=0,iwmax-1
     gloc(-iw-1,:,:) = TRANSPOSE(conjg(gloc(iw,:,:)))
  enddo

  gloc = gloc/dble(nkp)
end subroutine get_sigma_g_loc


! calculation of susceptibilities
!=========================================================================== 
subroutine get_chi0_loc(iwf, iwb, chi0_loc)
  implicit none
  integer :: i, j, k, l, i1, i2
  integer :: iwf, iwb
  complex(kind=8), intent(out) :: chi0_loc(ndim*ndim,ndim*ndim)

  chi0_loc = 0.d0
  i1 = 0
  do i=1,ndim
     do j=1,ndim
        i1=i1+1
        chi0_loc(i1,i1) = -beta*giw(iwf,i)*giw(iwf-iwb,j)
     enddo
  enddo

end subroutine get_chi0_loc


subroutine get_chi0_loc_inv(iwf, iwb, chi0_loc)
  implicit none
  integer :: i, j, k, l, i1, i2
  integer :: iwf, iwb
  complex(kind=8), intent(out) :: chi0_loc(ndim*ndim,ndim*ndim)

  chi0_loc = 0.d0
  i1 = 0
  do i=1,ndim
     do j=1,ndim
        i1=i1+1 ! orbital compound index {ji}
        chi0_loc(i1,i1) = -1.d0/(beta*giw(iwf,i)*giw(iwf-iwb,j))
     enddo
  enddo

end subroutine get_chi0_loc_inv


subroutine accumulate_chi0(ik, ikq, iwf, iwb, chi0)
  implicit none
  integer :: i, j, k, l, i1, i2
  integer :: iwf, iwb, ik, ikq
  complex(kind=8) :: g1(ndim,ndim), g2(ndim,ndim)
  complex(kind=8), intent(inout) :: chi0(ndim*ndim,ndim*ndim)
  complex(KIND=8) :: c
  complex(kind=8),parameter :: ci = (0d0,1d0)

  g1(:,:) = -hk(:,:,ik)
  do i=1,ndim
     g1(i,i) = ci*iw_data(iwf)+mu-hk(i,i,ik)-dc(1,i)-siw(iwf,i)
  enddo


  if (ndim .eq. 1) then
    g1(1,1)=1.d0/g1(1,1)
  else if (ndim .eq. 2) then
    call inverse_matrix_2(g1)
  else if (ndim .eq. 3) then
    call inverse_matrix_3(g1)
  else
    call inverse_matrix(g1)
  end if

  g2(:,:) = -hk(:,:,ikq)
  do i=1,ndim
     g2(i,i) = ci*iw_data(iwf-iwb)+mu-hk(i,i,ikq)-dc(1,i)-siw(iwf-iwb,i)
  enddo

  if (ndim .eq. 1) then
    g2(1,1)=1.d0/g2(1,1)
  else if (ndim .eq. 2) then
    call inverse_matrix_2(g2)
  else if (ndim .eq. 3) then
    call inverse_matrix_3(g2)
  else
    call inverse_matrix(g2)
  end if

  ! Accumulate chi0
  i2=0
  do l=1,ndim
     do k=1,ndim
        i2=i2+1
        i1 = 0
        do i=1,ndim
           c = - beta*g1(i,l)
           do j=1,ndim
              i1=i1+1
              ! the lattice green's functions are not orbital diagonal
              chi0(i1,i2) = chi0(i1,i2) + c*g2(k,j)
           enddo
        enddo
     enddo
  enddo

end subroutine accumulate_chi0


! calculation of occupations
!=========================================================================== 
subroutine get_ndmft()
  implicit none
  integer :: iw,i
  complex(kind=8) :: giw_sum(ndim)

  giw_sum = 0.d0
  n_dmft = 0.d0
  do iw=0,iwmax-1
    giw_sum(:) = giw_sum(:)+giw(iw,:)
  enddo
  n_dmft(:) = 2.d0*real(giw_sum(:))/beta+0.5d0

  if ((verbose .and. (index(verbstr,"Dmft") .ne. 0)) .and. mpi_wrank .eq. master ) then
    open(56, file=trim(output_dir)//"n_dmft.dat", status='unknown')
    write(56,*) '# n_dmft(i,i) [i=1,ndim]'
    write(56,'(100F12.6)') (real(n_dmft(i)),i=1,ndim)
    close(56)
  endif
end subroutine get_ndmft


subroutine get_nfock()
  implicit none
  integer :: ik,iw,i
  complex(kind=8)  :: gkiw(ndim,ndim)

  n_fock = 0.d0
  gkiw = 0.d0
  do ik=1,nkp
     do iw=0,iwmax-1
        call get_gkiw(ik, iw, 0, gkiw)
        n_fock(ik,:,:) = n_fock(ik,:,:)+real(gkiw(:,:))
     enddo
     n_fock(ik,:,:) = 2.d0*n_fock(ik,:,:)/beta
     do i=1,ndim
        n_fock(ik,i,i) = n_fock(ik,i,i)+0.5d0
     enddo
  enddo

  if (mpi_wrank .eq. master .and. text_output) then
    open(110, file=trim(output_dir)//"n_fock.dat", status='unknown')
    write(110,*) '# ik, kx, ky, kz, n_fock(k,i,i) [i=1,ndim]'
    do ik=1,nkp
      write(110,'(I8, 100F12.6)') ik, k_data(1,ik),k_data(2,ik),k_data(3,ik), (real(n_fock(ik,i,i)),i=1,ndim)
    enddo
    close(110)
  endif
end subroutine get_nfock

subroutine get_ndga(sigma_sum)
  implicit none
  integer :: ik,iw,i
  complex(kind=8), intent(in) :: sigma_sum(ndim,ndim,-iwfmax_small:iwfmax_small-1, nkp)
  complex(kind=8) :: gkiw(ndim,ndim)
  complex(kind=8) :: skiw(ndim,ndim)

  ! we calculate n_dga_k with the dga selfenergy in the vertex area and
  ! with the dmft selfenergy outside (diagonal)
  n_dga_k = 0.d0
  n_dga = 0.d0
  gkiw = 0.d0
  do ik=1,nkp
     do iw=0,iwfmax_small-1
        call get_gkiw_dga(ik, iw, sigma_sum(:,:,iw,ik), gkiw)
        n_dga_k(ik,:,:) = n_dga_k(ik,:,:)+real(gkiw(:,:))
     enddo
     do iw=iwfmax_small,iwmax-1
        skiw = 0.d0
        do i=1,ndim
          skiw(i,i) = siw(iw,i) ! orbital diagonal
        enddo
        call get_gkiw_dga(ik, iw, skiw(:,:), gkiw)
        n_dga_k(ik,:,:) = n_dga_k(ik,:,:)+real(gkiw(:,:))
     enddo
     n_dga_k(ik,:,:) = 2.d0*n_dga_k(ik,:,:)/beta
     do i=1,ndim
        n_dga_k(ik,i,i) = n_dga_k(ik,i,i)+0.5d0
        n_dga(i) = n_dga(i) + n_dga_k(ik,i,i)
     enddo
  enddo

  n_dga = n_dga/dble(nkp)

  if (mpi_wrank .eq. master .and. text_output) then
    open(110, file=trim(output_dir)//"n_dga_k.dat", status='unknown')
    write(110,*)  '# ik, kx, ky, kz, n_dga(k,i,i) [i=1,ndim]'
    do ik=1,nkp
      write(110,'(I8,100F12.6)') ik, k_data(1,ik),k_data(2,ik),k_data(3,ik), (real(n_dga_k(ik,i,i)),i=1,ndim)
    enddo
    close(110)

    open(111, file=trim(output_dir)//"n_dga.dat", status='unknown')
    write(111,*) '# n_dga(i,i) [i=1,ndim]'
    write(111,'(100F12.6)') (real(n_dga(i)),i=1,ndim)
    close(111)
  endif
end subroutine get_ndga

end module
