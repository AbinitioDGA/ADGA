! This file is part of the Abinitio Dynamical Vertex Approximation (ADGA)
! package. It is an electronic structure code which allows the inclusion of
! non-local correlations beyond DMFT and the calculation of momentum-dependent
! susceptibilities.
!
! The public repository can be found at
! https://github.com/AbinitioDGA/ADGA
!
! The arXiv publication can be found at
! https://arxiv.org/abs/1710.06651
!
! Copyright (C) <2017, 2018>
! <Anna Galler*, Patrick Thunström, Josef Kaufmann, Matthias Pickem, Jan M. Tomczak, Karsten Held>
! * Corresponding author. E-mail address: galler.anna@gmail.com
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

module eom_module

  use parameters_module
  use one_particle_quant_module
  use kq_tools
  implicit none

contains

!================================================================================================
   subroutine calc_eom_static(kq_ind_eom,iq,v,sigma_dmft,sigma_hf,nonlocal)
   implicit none
   integer,intent(in)             :: kq_ind_eom(nkp_eom,nqp)
   integer,intent(in)             :: iq
   complex(kind=8),intent(in)     :: v(ndim2,ndim2)
   complex(kind=8),intent(inout)  :: sigma_dmft(ndim,ndim,-iwfmax_small:iwfmax_small-1)
   complex(kind=8),intent(inout)  :: sigma_hf(ndim,ndim,nkp_eom)
   logical,intent(in)             :: nonlocal
   complex(kind=8)                :: sigma_tmp(ndim,ndim)
   integer :: dum,i,j,iw,iwf,iwf2,l,k,ik,ikq
   integer :: i1,i2,i3,i4
   complex(kind=8) :: alpha, delta

   ! First of all, do the Hartree local part:
   if(iq==1) then
      !compute local Hartree-Fock contribution and add it to sigma_dmft:
      ! u(i,j,k,l) = u( i1 = {ki}, i2 = {jl} ),
      ! so u(i,j,i,j) = u( i1 = {ii}, i2 = {jj} )
      ! and u(i,j,j,i) = u( i1 = {ji}, i2 = {ji} ) (only for equal spins)
      sigma_tmp = 0
      do i=1,ndim
         i1 = i+(i-1)*ndim ! = {ii}
         do j=1,ndim
            i2 = j+(j-1)*ndim ! = {jj}
            i3 = j+(i-1)*ndim ! = {ji}
            sigma_tmp(i,i) = sigma_tmp(i,i) + (2d0*u(i1,i2)-u(i3,i3))*n_dmft(j)
         enddo
      enddo
      ! Add it to sigma_dmft
      do i=1,ndim
         sigma_dmft(i,i,:) = sigma_dmft(i,i,:) + sigma_tmp(i,i)
      enddo
      ! Calculate the non-local hartree contribution to the self-energy
      if (nonlocal .and. do_vq) then
         ! 1. Compute non-local Hartree-contribution 2*v(q=0)*n_DMFT and add it to sigma_hf (only for q = 0)
         ! vq(i,j,k,l) = vq( i1 = {ki}, i2 = {jl} ),
         ! so vq(i,j,i,j) = vq( i1 = {ii}, i2 = {jj} ), (and vq(i,j,j,i) = vq( i1 = {ji}, i2 = {ji} ) )
         ! Nb: n_dmft is NOT diagonal if we have ineq atoms!
         do i=1,ndim
            i1 = i+(i-1)*ndim ! = {ii}
            do j=1,ndim
               i2 = j+(j-1)*ndim ! = {jj}
               sigma_hf(i,i,:) = sigma_hf(i,i,:) + 2.d0*v(i1,i2)*n_dmft(j) ! FIXME: non-diagonal n_dmft
            enddo
         enddo
      endif
   endif


   if (nonlocal .and. do_vq) then
      ! 1. Compute non-local Fock-contribution v(q)*n_fock(k-q) and add it to sigma_hf
      i2=0
      do l=1,ndim
         do j=1,ndim
            i2=i2+1 ! = {jl}
            i1 = 0
            do i=1,ndim
               do k=1,ndim
                  i1 = i1+1 ! = {ki}
                  ! vq(i,j,k,l) = v(i1 = {ki},i2 = {jl})
                  alpha = v(i1,i2)/nqp
                  do ik=1,nkp_eom
                     ikq = kq_ind_eom(ik,iq)
                     sigma_hf(i,l,ik) = sigma_hf(i,l,ik) - alpha*n_fock(ikq,j,k)
                  enddo
               enddo
            enddo
         enddo
      enddo
   endif

   return
end subroutine calc_eom_static
!=============================================================================================

!=============================================================================================
   subroutine calc_eom_dmft(gammawd,gammawm,iwb,sigma_dmft)
   implicit none
   complex(kind=8),intent(in)     :: gammawd(ndim2,maxdim)
   complex(kind=8),intent(in)     :: gammawm(ndim2,maxdim)
   integer,intent(in)             :: iwb
   complex(kind=8),intent(inout)  :: sigma_dmft(ndim,ndim,-iwfmax_small:iwfmax_small-1)
   integer :: dum,i,j,iw,iwf,iwf2,l,k,ik,ikq
   integer :: i1,i2,i3,i4
   complex(kind=8) :: m_tot(ndim2,maxdim)
   complex(kind=8) :: alpha, delta
   complex(kind=8),allocatable :: gamma_dmft(:,:)

   ! First of all, do the purely local part:
   !DMFT part (for the computation of the DMFT self energy through the local EOM, for comparison)
   if (debug .and. (index(dbgstr,"Dmft_phbar") .ne. 0)) then
      ! sig_phbar = -1/beta*U_tilde*(-1/2*gamma^w_d-3/2*gamma^w_m)
      allocate(gamma_dmft(ndim2,maxdim))
      gamma_dmft = 0.5d0*gammawd + 1.5d0*gammawm
      alpha = 1.d0/beta
      delta = 0.d0
      call zgemm('n', 'n', ndim2, maxdim, ndim2, alpha, u_tilde, ndim2, gamma_dmft, ndim2, delta, m_tot, ndim2)
      deallocate(gamma_dmft)
   else
      ! sig_ph = -1/beta*U*gamma^w_d
      alpha = -1d0/beta
      delta = 0.d0
      call zgemm('n', 'n', ndim2, maxdim, ndim2, alpha, u, ndim2, gammawd, ndim2, delta, m_tot, ndim2)
   endif

   !break up the compound index: Nb: j == k since G_loc is diagonal
   i2 = ndim2*iwfmax_small ! offset
   do iwf=0,iwfmax_small-1
      do l=1,ndim
         do k=1,ndim
            i2 = i2+1 ! = {kl,iwf}
            do i=1,ndim
               i1 = k+(i-1)*ndim ! = {ki}
               ! m_tot_array(i,k,k,l,iwf2) = m_tot(i1,i2)
               sigma_dmft(i,l,iwf) = sigma_dmft(i,l,iwf)+m_tot(i1,i2)*giw(iwf-iwb,k)
            enddo
         enddo
      enddo
      sigma_dmft(:,:,-iwf-1) = TRANSPOSE(conjg(sigma_dmft(:,:,iwf)))
   enddo

   return
end subroutine calc_eom_dmft
!=============================================================================================

!=============================================================================================
   subroutine calc_eom_dynamic(etaqd,etaqm,gammawd,gammaqd,kq_ind_eom,iwb,iq,v,sigma_nl)
   implicit none
   complex(kind=8),intent(in)     :: etaqd(ndim2,maxdim)
   complex(kind=8),intent(in)     :: etaqm(ndim2,maxdim)
   complex(kind=8),intent(in)     :: gammawd(ndim2,maxdim)
   complex(kind=8),intent(in)     :: gammaqd(ndim2,maxdim)
   integer,intent(in)             :: kq_ind_eom(nkp_eom,nqp)
   integer,intent(in)             :: iwb
   integer,intent(in)             :: iq
   complex(kind=8),intent(in)     :: v(ndim2,ndim2)
   complex(kind=8),intent(inout)  :: sigma_nl(ndim,ndim,-iwfmax_small:iwfmax_small-1,nkp_eom)
   integer :: dum,i,j,iw,iwf,iwf2,l,k,ik,ikq
   integer :: i1,i2,i3,i4
   complex(kind=8) :: m_tot(ndim2,maxdim)
   complex(kind=8) :: u_work(ndim2,ndim2)
   complex(kind=8) :: gkiw(ndim,ndim)
   complex(kind=8) :: alpha, delta

   !EQUATION OF MOTION:
   m_tot = 0.d0

   !eta term: density part: -1/beta*(U + V^q - U_tilde/2).eta^q_d
   u_work = v + u - 0.5d0*u_tilde
   alpha = -1.d0/beta/dble(nqp)
   delta = 0.d0
   call zgemm('n', 'n', ndim2, maxdim, ndim2, alpha, u_work, ndim2, etaqd, ndim2, delta, m_tot, ndim2)

   !eta term: magnetic part: -1/beta*(-1.5*U_tilde).eta^q_m
   alpha = 1.5d0/beta/dble(nqp)
   delta = 1.d0  ! Add to m_tot
   call zgemm('n', 'n', ndim2, maxdim, ndim2, alpha, u_tilde, ndim2, etaqm, ndim2, delta, m_tot, ndim2)

   !local part: -1/beta*(-U).gamma^q_d
   alpha = 1.d0/beta/dble(nqp)
   delta = 1.d0  ! Add to m_tot
   call zgemm('n', 'n', ndim2, maxdim, ndim2, alpha, u, ndim2, gammaqd, ndim2, delta, m_tot, ndim2)

   !v(q) part -1/beta*V.gamma^w_d
   if(do_vq) then
      alpha = -1.d0/beta/dble(nqp)
      delta = 1.d0
      call zgemm('n', 'n', ndim2, maxdim, ndim2, alpha, v, ndim2, gammawd, ndim2, delta, m_tot, ndim2)
   endif

   !compute k-dependent self energy (convolution with Greens function gkiw):
   gkiw = 0.d0
   do ik=1,nkp_eom
      ikq = kq_ind_eom(ik,iq)
      i2 = ndim2*iwfmax_small ! offset
      do iwf=0,iwfmax_small-1
         call get_gkiw(ikq, iwf, iwb, gkiw)
         do l=1,ndim
            do k=1,ndim
               i2 = i2+1  ! compound index {kl,iwf} (k fastest)
               i1 = 0
               do i=1,ndim
                  do j=1,ndim
                     i1 = i1+1 ! compound index {ji} (j fastest)
                     sigma_nl(i,l,iwf,ik) = sigma_nl(i,l,iwf,ik)+m_tot(i1,i2)*gkiw(j,k)
                  enddo
               enddo
            enddo
         enddo
         sigma_nl(:,:,-iwf-1,ik) = transpose(conjg(sigma_nl(:,:,iwf,ik)))
      enddo !iwf
   enddo !ik

   return
end subroutine calc_eom_dynamic
!=============================================================================================

!=============================================================================================
subroutine add_siw_dmft(sigma_sum)
  implicit none
  complex(kind=8) :: sigma_sum(ndim, ndim, -iwfmax_small:iwfmax_small-1, nkp_eom)
  integer :: ik, iwf, iband

 ! local contribution is replaced by the DMFT self energy for better asymptotics
    do ik=1,nkp_eom
       do iwf=-iwfmax_small,iwfmax_small-1
          do iband=1,ndim
             sigma_sum(iband, iband, iwf, ik) = sigma_sum(iband, iband, iwf, ik) + siw(iwf, iband)
          enddo
       enddo
    enddo

 end subroutine add_siw_dmft
!=============================================================================================

!==============================================================================================
subroutine output_eom(sigma_sum, sigma_sum_dmft, sigma_sum_hf, sigma_loc, gloc, nonlocal)
   implicit none

   complex(kind=8), intent(in) :: sigma_sum(ndim, ndim, -iwfmax_small:iwfmax_small-1, nkp_eom)
   complex(kind=8), intent(in) :: sigma_sum_dmft(ndim, ndim, -iwfmax_small:iwfmax_small-1)
   complex(kind=8), intent(in) :: sigma_sum_hf(ndim,ndim,nkp_eom)
   complex(kind=8), intent(in) :: sigma_loc(ndim, ndim, -iwfmax_small:iwfmax_small-1)
   complex(kind=8) :: sigma_tmp(ndim*(ndim+1)/2)
   complex(kind=8), intent(in) :: gloc(-iwmax:iwmax-1,ndim,ndim)
   logical,intent(in) :: nonlocal ! Print the non-local data
   integer :: ik, iwf, i, j, iband,ii, i1, i2, i3, i4
   character(len=50) :: eom_format
   character(len=200) :: filename_siwk
   integer            :: er
   character(len=200) :: erstr

   if (nonlocal .and. .not. k_path_eom) then
     open(44, file=trim(output_dir)//"siw_loc_diag.dat", status='unknown')
     write(44,'(A)') '# wf, (real(siw_dga_loc(wf,i,i)), imag(siw_dga_loc(wf,i,i))) [i=1,ndim]'
     do iwf=-iwfmax_small,iwfmax_small-1
        write(44,'(100F12.6)') iw_data(iwf), (real(sigma_loc(i,i,iwf)), aimag(sigma_loc(i,i,iwf)), i=1,ndim)
     enddo
     close(44)
   endif

   open(50, file=trim(output_dir)//"siw_dmft_rep_diag.dat", status='unknown')
   write(50,'(A)') '# wf, (real(siw_dmft_rep(wf,i,i)), imag(siw_dmft_rep(wf,i,i))) [i=1,ndim]'
   do iwf=-iwfmax_small,iwfmax_small-1
      write(50,'(100F12.6)') iw_data(iwf), (real(sigma_sum_dmft(i,i,iwf)),aimag(sigma_sum_dmft(i,i,iwf)), i=1,ndim)
   enddo
   close(50)

   if (nonlocal .and. .not. k_path_eom) then
     open(46, file=trim(output_dir)//"g_loc_diag.dat", status='unknown')
     write(46,'(A)') '# wf, (real(giw_dga_loc(iwf,i,i)), imag(giw_dga_loc(iwf,i,i))) [i=1,ndim]'
     do iwf=-iwmax,iwmax-1
        write(46,'(100F12.6)') iw_data(iwf), (real(gloc(iwf,i,i)), aimag(gloc(iwf,i,i)), i=1,ndim)
     enddo
     close(46)
   endif

   if (ndim .ge. 2) then

     if (nonlocal .and. .not. k_path_eom) then
       open(44, file=trim(output_dir)//"siw_loc_full.dat", status='unknown')
       write(44,'(A)') '# wf, (real(siw_dga_loc(wf,i,j)), imag(siw_dga_loc(iwf,i,j))) [j=i,ndim] [i=1,ndim]'
       do iwf=-iwfmax_small,iwfmax_small-1
         write(44,'(100F12.6)') iw_data(iwf), ((real(sigma_loc(i,j,iwf)), aimag(sigma_loc(i,j,iwf)), j=i,ndim), i=1,ndim)
       enddo
       close(44)
     endif

     open(50, file=trim(output_dir)//"siw_dmft_rep_full.dat", status='unknown')
     write(50,'(A)') '# wf, (real(siw_dmft_rep(wf,i,j)), imag(siw_dmft_rep(iwf,i,j))) [j=i,ndim] [i=1,ndim]'
     do iwf=-iwfmax_small,iwfmax_small-1
       write(50,'(100F12.6)') iw_data(iwf), ((real(sigma_sum_dmft(i,j,iwf)),aimag(sigma_sum_dmft(i,j,iwf)), j=i,ndim), i=1,ndim)
     enddo
     close(50)

     if (nonlocal .and. .not. k_path_eom) then
       open(46, file=trim(output_dir)//"g_loc_full.dat", status='unknown')
       write(46,'(A)') '# wf, (real(giw_dga_loc(wf,i,j)), imag(giw_dga_loc(wf,i,j))) [j=i,ndim] [i=1,ndim]'
       do iwf=-iwmax,iwmax-1
         write(46,'(100F12.6)') iw_data(iwf), ((real(gloc(iwf,i,j)), aimag(gloc(iwf,i,j)), j=i,ndim), i=1,ndim)
       enddo
       close(46)
     endif
   endif

   ! If we only run the dfmt part, then it is highly likely that we want to compare sig_dmft_rep to the original data.
   if (.not. nonlocal .or. (verbose .and. (index(verbstr,"Dmft") .ne. 0))) then
      open(44, file=trim(output_dir)//"siw_dmft_orig.dat", status='unknown')
      write(44,'(A)') '# wf, (real(siw_dmft(wf,i)), imag(siw_dmft(wf,i))) [i=1,ndim]'
      do iwf=-iwfmax_small,iwfmax_small-1
         write(44,'(100F12.6)') iw_data(iwf), (real(siw(iwf,i)), aimag(siw(iwf,i)), i=1,ndim)
      enddo
      close(44)
      open(44, file=trim(output_dir)//"g_dmft_orig.dat", status='unknown')
      write(44,'(A)') '# wf, (real(giw_dmft(wf,i)), imag(giw_dmft(wf,i))) [i=1,ndim]'
      do iwf=-iwmax,iwmax-1
         write(44,'(100F12.6)') iw_data(iwf), (real(giw(iwf,i)), aimag(giw(iwf,i)), i=1,ndim)
      enddo
      close(44)
   endif

   if (nonlocal) then
      if (verbose .and. (index(verbstr,"Siwk") .ne. 0)) then
       do ik=1,nkp_eom
         write(filename_siwk,'(A,F5.3,A,F5.3,A,F5.3,A)') 'siwk_',k_data(1,k_data_eom(ik)),'_', &
           k_data(2,k_data_eom(ik)),'_',k_data(3,k_data_eom(ik)),'.dat'
         open(34, file=trim(output_dir)//filename_siwk)
         write(34,'(A)') '# wf, (real(siwk_dga(wf,ik,i,j)), imag(siwk_dga(wf,ik,i,j))) [j=i,ndim] [i=1,ndim]'
         do iwf=-iwfmax_small,iwfmax_small-1
           write(34,'(100F12.6)') iw_data(iwf), ((real(sigma_sum(i,j,iwf,ik)), aimag(sigma_sum(i,j,iwf,ik)), j=i,ndim), i=1,ndim)
         end do
         close(34)
       end do
      endif


      open(45, file=trim(output_dir)//"siw_all_k.dat",status='unknown')
        write(45,'(A,A)') '# ik, kx, ky, kz, wf, ', &
                        '(real(siwk_dga(wf,ik,i,j)), imag(siwk_dga(wf,ik,i,j)), real(siwk_hf(ik,i,j))) [j=i,ndim] [i=1,ndim]'
        do ik=1,nkp_eom
           do iwf=-iwfmax_small,iwfmax_small-1
              write(45,'(I8, 100F12.6)') ik, k_data(1,k_data_eom(ik)), k_data(2,k_data_eom(ik)), &
                k_data(3,k_data_eom(ik)), iw_data(iwf), &
                ((real(sigma_sum(i,j,iwf,ik)), aimag(sigma_sum(i,j,iwf,ik)), real(sigma_sum_hf(i,j,ik)),j=i,ndim), i=1,ndim)
           enddo
        enddo
      close(45)

   endif

   return
end subroutine output_eom
!================================================================================================

end module
