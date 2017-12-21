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

program main

#ifdef MPI
  use mpi
#endif

  use hdf5_module
  use lapack_module
  use parameters_module
  use config_module
  use one_particle_quant_module
  use eom_module
  use susc_module
  use kq_tools
  use interaction_module
  use aux
  use mpi_org

  implicit none
  integer :: iw, ik, iq, ikq, iwf, iwb, iv, dum, dum1, iwf1, iwf2
  integer :: i, j, k, l, n, i1, i2, i3, i4
  complex(kind=8), allocatable :: gloc(:,:,:)
  complex(kind=8), allocatable :: chi0w(:,:,:), chi0w_inv(:,:,:), chi0(:,:), chi0q(:,:,:), chi0nl(:,:,:)
  complex(kind=8), allocatable :: chi0wFd_slice(:,:)
  complex(kind=8), allocatable ::  chi0wFm_slice(:,:), chi0wFd(:,:), chi0wFm(:,:), chi_loc(:,:)
  complex(kind=8), allocatable ::  oneplusgammawm(:,:), oneplusgammawd(:,:), gammawd(:,:), gammawm(:,:)
  complex(kind=8), allocatable :: chi_qw_dens(:,:,:),chi_qw_magn(:,:,:),bubble_nl(:,:,:),chi_qw_full(:,:,:)
  complex(kind=8), allocatable :: chi_loc_dens(:,:,:),chi_loc_magn(:,:,:),bubble_loc(:,:,:),bubble_loc_tmp(:,:)
  complex(kind=8), allocatable :: chi_loc_qmc(:,:,:)
  integer, allocatable :: kq_ind(:,:), qw(:,:)
  complex(kind=8), allocatable :: bigwork_magn(:,:), etaqd(:,:), etaqm(:,:), rectanglework(:,:)
  complex(kind=8), allocatable :: smallwork(:,:), bigwork_dens(:,:)
  real(kind=8 ):: start, finish
  real(kind=8 ):: tstart, tfinish, timings(7) ! timings: local, chi0, gammas, eta inversion, eta construct, eom, chi
  complex(kind=8) :: alpha, delta
  integer :: iqw
  integer :: offset
  logical :: update_chi_loc_flag
  complex(kind=8), allocatable :: gammaqd(:,:), v(:,:)
  complex(kind=8), allocatable :: sigma_nl(:,:,:,:), sigma_hf(:,:,:), sigma_dmft(:,:,:)
  complex(kind=8), allocatable :: sigma_sum(:,:,:,:), sigma_sum_hf(:,:,:), sigma_sum_dmft(:,:,:), sigma_loc(:,:,:)
! variables for date-time string
  character(20) :: date,time,zone,iwb_str
  character(200) :: output_filename,erstr ! Filename and error string
  integer,dimension(8) :: time_date_values
  logical :: verbose_extra
  integer :: er ! error flag
  logical :: nonlocal ! Do the nonlocal quantities
#ifdef MPI
  character(len=MPI_MAX_PROCESSOR_NAME) :: hostname
#endif

  ! mpi initialization
  call mpi_initialize()

  ! read command line argument -> file name of config file
  if (iargc() .ne. 1) then
    call mpi_stop('The program has to be executed with exactly one argument. (Name of config file)')
  end if

  ! read config settings 
  call read_config(er,erstr)
  if (er .ne. 0) call mpi_stop(erstr,er)

  ! create output folder if not yet existing
  if (mpi_wrank .eq. master) call system('mkdir -p ' // trim(adjustl(output_dir)))
  call mpi_barrier(mpi_comm_world,ierr)

  ! Set up the output file
  verbose_extra = (verbose .and. (index(verbstr,"Extra") .ne. 0))
  if (mpi_wrank .eq. master .or. verbose_extra) then
     ounit = 123
  else
     ounit = 0 ! this is fortrans unit for the stderr
  endif
  if (ounit .gt. 0) then
     if (mpi_wrank .eq. master) then
        output_filename = trim(adjustl(output_dir))//"out"
     else
        write(output_filename,*) mpi_wrank
        output_filename = trim(adjustl(output_dir))//"out."//TRIM(ADJUSTL(output_filename))
     endif
     open(unit=ounit,file=TRIM(ADJUSTL(output_filename)),status='unknown')
  endif

  ! check config settings
  call check_config(er,erstr) 
  if (er .ne. 0) call mpi_stop(erstr,er)


  ! program introduction
  if (ounit .ge. 1) then
    call date_and_time(date,time,zone,time_date_values)
    call mpi_get_processor_name(hostname,i,j)
    write(ounit,*)
    write(ounit,*)                        '/------------------------------------------------------------------\'
    write(ounit,*)                        '|  Ab initio dynamical vertex approximation program (abinitiodga)  |'
    write(ounit,'(1x,a,i11,a,3i6,a)') '|  Running on ',mpi_wsize,' core(s) with ',nkpx, nkpy, nkpz,' k-points |'
    write(ounit,*)                        '|     time             date           host                         |'
    write(ounit,'(" | ",a,7x,a,10x,a,26x,"|")') trim(time),trim(date),hostname(1:i)
    write(ounit,*)                        '\------------------------------------------------------------------/'
    write(ounit,*)
    if (verbose) write(ounit,*) "Verbose string: ",trim(ADJUSTL(verbstr))
    if (debug) write(ounit,*)   "Debug string:   ",trim(ADJUSTL(dbgstr))
  end if


  ! creation of hdf5 output file
  if (mpi_wrank .eq. master) then
    ! generate a date-time string for output file name
    output_filename=trim(output_dir)//'adga-'//trim(date)//'-'//trim(time)//'-output.hdf5'
    ! while the name is generated here to get the correct starting time, 
    ! the file is created later, just before the beginning of the parallel loop.
    ! Thus, the parameters can be written immediately.
  end if


!##################  READ W2DYNAMICS HDF5 OUTPUT FILE ##################
! open the hdf5-fortran interface
  call init_h5()

! read bosonic and fermionic Matsubara axes of one- and two-particle data
  call get_freq_range(iwmax,iwfmax,iwbmax,n3iwf,n3iwb,n2iwb)
  call check_freq_range(mpi_wrank,master,er)
  if (er .ne. 0) call mpi_stop('Frequency range error.')

! after frequencies and dimensions are obtained, arrays can be allocated 
  allocate(siw(-iwmax:iwmax-1,ndim))
  allocate(giw(-iwmax:iwmax-1,ndim))
  allocate(dc(2,ndim)) ! indices: spin band

! read Hamiltonian
  if(read_ext_hk) then
    call read_hk_w2w(er,erstr)
    if (er .ne. 0) call mpi_stop(erstr,er)
  else
    call read_hk_w2dyn(er,erstr)
    if (er .ne. 0) call mpi_stop(erstr,er)
  end if

  call read_mu()   ! w2d chemical potential
  call read_beta() ! w2d inverse temperature
  call read_dc()   ! w2d double counting


  if (ounit .ge. 1) then
    write(ounit,*) 'orb_symmetry = ', orb_sym
  end if

  if (ounit .ge. 1) then
    write(ounit,*) 'beta=', beta
    write(ounit,*) 'mu=', mu
    write(ounit,*) 'dc=', dc
  end if

  call config_init(er,erstr) ! this requires beta, so it has to be called after read_beta()
  if (er .ne. 0) call mpi_stop(erstr,er)
  ! COMPUTE or READ the local single-particle Greens function:
  ! The calculation is necessary if one wants to do calculations with p-bands
  ! since read_giw reads the giw array which only ! contains correlated bands
  call read_siw()  ! w2d self energy
  if(exist_p .or. (debug .and. (index(dbgstr,"Makegiw") .ne. 0))) then
    if (ounit .ge. 1) write(ounit,*) 'Constructig giw from siw. (The QMC self-energy + Ham.hk) '
    call get_giw() ! writes giw_calc.dat if we have the verbose keyword "Dmft"
  else
    if (ounit .ge. 1) write(ounit,*) "Reading giw from file. (The QMC green's function) "
    call read_giw()  ! w2d greens function G_dmft
  endif

  call finalize_h5() ! close the hdf5-fortran interface


  !read umatrix from separate file:
  if (read_ext_u) then
    if (ounit .ge. 1) then
      write(ounit,*) 'Reading the U matrix from file.'
      write(ounit,*) 'U matrix in ', filename_umatrix
    endif
    call read_u(u,u_tilde)
  else
    if (ounit .ge. 1) write(ounit,*) 'Creating U matrix from input parameters.'
    call create_u(u,u_tilde)

    ! test umatrix
    if (mpi_wrank .eq. master .and.  (verbose .and. (index(verbstr,"Umatrix") .ne. 0))) then
      open(unit=10,file=trim(output_dir)//"umatrix.dat")
      write(10,*) 'Umatrix File for the abinitiodga code: band,band,band,band,Uvalue'
      do i=1,ndim
      do j=1,ndim
      do k=1,ndim
      do l=1,ndim
        write(10,'(4I10,F15.8)') i,j,k,l,Umat(i,j,k,l)
      enddo
      enddo
      enddo
      enddo
    endif
    deallocate(Umat)
  endif

  allocate(n_dmft(ndim), n_fock(nkp,ndim,ndim), n_dga(ndim), n_dga_k(nkp,ndim,ndim))

!compute DMFT filling n_dmft
  call get_ndmft() ! writes n_dmft.dat if we have the verbose keyword "Dmft"
!compute k-dependent filling for Fock-term (computed in the EOM):
  call get_nfock() ! writes n_fock.dat

  ! small arrays
  allocate(v(ndim2,ndim2))
  allocate(smallwork(ndim2,ndim2))
  allocate(chi0(ndim2,ndim2))
  ! medium arrays
  allocate(chi0q(ndim2,ndim2,-iwfmax_small:iwfmax_small-1))
  allocate(chi0w_inv(ndim2,ndim2,-iwfmax_small:iwfmax_small-1))
  allocate(gammawd(ndim2,maxdim), gammawm(ndim2,maxdim))
  allocate(oneplusgammawd(ndim2,maxdim))
  allocate(oneplusgammawm(ndim2,maxdim))
  allocate(chi0wFd_slice(ndim2,maxdim))
  allocate(chi0wFm_slice(ndim2,maxdim))
  allocate(etaqd(ndim2,maxdim))
  allocate(etaqm(ndim2,maxdim))
  allocate(rectanglework(ndim2,maxdim))
  ! Huge arrays:
  ! One can save memory by treating one channel at the time, at
  ! the (rather large) cost of recalculating chi_0^q.
  ! A more scalable idea is to use scalapack for the inversion.
  ! A third option is to rewrite the equations to save space 
  ! if we have a lot of non-correlated atoms.
  allocate(chi0wFm(maxdim,maxdim))
  allocate(chi0wFd(maxdim,maxdim))

  if (q_path_susc .and. do_chi .and. (.not. q_vol)) then
    if (do_eom) call mpi_stop('Error: it is currently not possible to use both do_eom and q_path_susc',er)
    call qdata_from_file()
  else
    nqp=nqpx*nqpy*nqpz
    allocate(q_data(nqp))
    call generate_q_vol(nqpx,nqpy,nqpz,q_data)
  end if



  if (ounit .ge. 1) then
    write(ounit,'(1x)')
    if (.not. do_vq) then
      write(ounit,*) 'Running the calculation without V(q)'
    else
      write(ounit,*) 'Running the calculation with V(q)'
      write(ounit,*) 'V(q) data in ', filename_vq
    endif
    write(ounit,'(1x)')
    write(ounit,'(1x,"Frequency information:")')
    write(ounit,*) 'iwmax=', iwmax, ' (number of fermionic matsubara frequencies of one-particle quantities)'
    write(ounit,*) 'iwfmax=', iwfmax, 'iwfmax_small=', iwfmax_small,&
               ' (number of fermionic matsubara frequencies of two-particle quantities)'
    write(ounit,*) 'iwbmax=',iwbmax, 'iwbmax_small=', iwbmax_small, &
               ' (number of bosonic matsubara frequencies of two-particle quantities)'
    write(ounit,'(1x,"k-point information:")')
    write(ounit,*) nkp,'k-points in the mesh'
    if (q_vol) then
      write(ounit,*) nqp,'q-points in the mesh'
    else
      write(ounit,*) nqp,'q-points in the q-path'
    end if
    write(ounit,'(1x)')
  end if


  ! calculate the index of all \vec{k} - \vec{q}
  allocate(kq_ind(nkp,nqp))
  call cpu_time(start)
  call index_kq(kq_ind) ! new method
  call cpu_time(finish)
  if (ounit .ge. 1 .and. (verbose .and. (index(verbstr,"Time") .ne. 0))) then
   write(ounit,*)'TIME: finding k-q index:', finish-start
  endif
!##################### parallel code ##################################################

! define qw compound index for mpi:
  call mpi_distribute()

  allocate(qw(2,nqp*(2*iwbmax_small+1)))
  i1=0
  do iwb=-iwbmax_small,iwbmax_small
     do iq=1,nqp
        i1 = i1+1
        qw(1,i1) = iwb
        qw(2,i1) = iq
     enddo
  enddo

!distribute the qw compound index:
if (do_chi) then
  allocate(chi0w(ndim2,ndim2,-iwfmax_small:iwfmax_small-1))
  allocate(chi0nl(ndim2,ndim2,-iwfmax_small:iwfmax_small-1))
  allocate(chi_qw_dens(ndim2,ndim2,qwstart:qwstop),chi_qw_magn(ndim2,ndim2,qwstart:qwstop))
  allocate(bubble_nl(ndim2,ndim2,qwstart:qwstop))
  allocate(bubble_loc(ndim2,ndim2,-iwbmax_small:iwbmax_small))
  allocate(bubble_loc_tmp(ndim2,ndim2))
  allocate(chi_loc_dens(ndim2,ndim2,-iwbmax_small:iwbmax_small),chi_loc_magn(ndim2,ndim2,-iwbmax_small:iwbmax_small))
  chi_qw_dens=0.d0
  chi_qw_magn=0.d0
  bubble_nl=0.d0
  bubble_loc=0.d0
  bubble_loc_tmp=0.d0
  chi_loc_dens = 0
  chi_loc_magn = 0
  if (external_chi_loc) then
     ! We read the local susceptibilities from file.
     if (ounit .ge. 1)  write(ounit,*) "Reading the local (frequency summed) susceptibilities from file."
     call read_chi_loc(chi_loc_dens,'dens')
     call read_chi_loc(chi_loc_magn,'magn')
  else if (ounit .ge. 1) then
     write(ounit,*) "Constructing the local (frequency summed) susceptibilities from chi^{qvv'}_loc."
  endif
end if

if (external_threelegs) then
  if (ounit .ge. 1) write(ounit,*) "Reading the local threeleg vertex gamma^w from file."
end if

if (do_eom) then
  allocate(gammaqd(ndim2,maxdim))
  allocate(sigma_nl(ndim,ndim,-iwfmax_small:iwfmax_small-1,nkp), sigma_hf(ndim,ndim,nkp))
  allocate(sigma_dmft(ndim,ndim,-iwfmax_small:iwfmax_small-1))
  gammaqd = 0d0
  sigma_nl = 0d0
  sigma_hf = 0d0
  sigma_dmft = 0d0
end if


if (mpi_wrank.eq. master .and. (verbose .and. (index(verbstr,"Kpoints") .ne. 0))) then
  open(unit=256,file=trim(output_dir)//'kdata',status='replace')
  write(256,*) '### kx(ik), ky(ik), kz(ik)'
  do ik=1,nkp
    write(256,*) k_data(:,ik)
  end do
  close(256)
  open(unit=266,file=trim(output_dir)//'qdata',status='replace')
  write(266,*) '### kx(iq), ky(iq), kz(iq)'
  do iq=1,nqp
    write(266,*) k_data(:,q_data(iq))
  end do
  close(266)
end if



if (mpi_wrank .eq. master) then
  write(ounit,*) 'Writing hdf5 output to ',output_filename
  call init_h5_output(output_filename)
end if

  nonlocal = .not. (debug .and. (index(dbgstr,"Onlydmft") .ne. 0)) ! Do the non-local quantities!
  if (ounit .gt. 0) then
    write(ounit,'(1x)')
    if (nonlocal) then
       write(ounit,*) "Starting the main loop:"
    else
       write(ounit,*) "Starting the main loop: (Debug: Only local quantities are calculated!)"
    endif
    if (.not. (verbose .and. (index(verbstr,"Noprogress") .ne. 0))) then
      write(ounit,*) "To supress the progress indicator, use the verbose keyword Noprogress."
      write(ounit,'(1x,"-------------------------------------------------------------------------------------------")')
    endif
    call flush(ounit)
  endif

  timings = 0 ! Initialize the timings array
  iwb = iwbmax_small+3
  do iqw=qwstart,qwstop
     call cpu_time(start)
     ! TIME: local 
     tstart = start

     update_chi_loc_flag = qw(1,iqw) .ne. iwb
     iq = qw(2,iqw)
     iwb = qw(1,iqw)


     !read nonlocal interaction v and go into compound index:
     if(do_vq) then
        call read_vq(iq,v,er,erstr)
        if (er .ne. 0) call mpi_stop(erstr,er)
       ! v = v-u  !otherwise, local U would be included twice
     else
        v = 0.d0
     endif

     !update chi_loc only if iwb is different than the previous one:
     if(update_chi_loc_flag) then

        ! compute local bubble chi0w^{-1}(i1,i2)(orbital compound index i1,i2):
        do iwf=-iwfmax_small,iwfmax_small-1
           call get_chi0_loc_inv(iwf, iwb, chi0w_inv(:,:,iwf))
        enddo
        if (do_chi) then
           ! Use the big box
           bubble_loc_tmp = 0 ! Initialize!
           do iwf=iwstart,iwstop
              call get_chi0_loc(  iwf, iwb, chi0)
              bubble_loc_tmp = bubble_loc_tmp + chi0 ! Accumulate chi0 over iwf
              if (-iwfmax_small .le. iwf .and. iwf .le. iwfmax_small-1) chi0w(:,:,iwf) = chi0 ! only in the small box
           enddo
           bubble_loc_tmp = bubble_loc_tmp/(beta**2) ! Nb: bubble_loc_tmp is also used to obtain bubble_nl = bubble^q - bubble^w 
           if (iq .eq. 1) bubble_loc(:,:,iwb) = bubble_loc(:,:,iwb) + bubble_loc_tmp ! Accumulate bubble_loc only when iq = 1.
        end if

        ! Read chi^w_dens and chi^w_magn (temporarily stored in chi0wF)
        call read_vertex(chi0wFd,chi0wFm,iwb)

        !time reversal symmetry (which is simply a transpose in our compound index)
        ! TODO: Move this to the preprocessing step!
        do i1=1,maxdim
           do i2=i1+1,maxdim
              chi0wFm(i1,i2) = 0.5d0*(chi0wFm(i1,i2)+chi0wFm(i2,i1))
              chi0wFm(i2,i1) = chi0wFm(i1,i2)

              chi0wFd(i1,i2) = 0.5d0*(chi0wFd(i1,i2)+chi0wFd(i2,i1))
              chi0wFd(i2,i1) = chi0wFd(i1,i2)
           enddo
        enddo



        ! Calculate chi0^w.F
        do dum = 0,2*iwfmax_small-1
           iwf = dum - iwfmax_small
           !compute 1+chi0^w.F = chi_loc*chi0_loc_inv (Nb: Here we assume that chi0_loc_inv is diagonal in the compound index):
           do i2=1,ndim2
              i = i2+dum*ndim2 ! compound index (i2,iwf)
              chi0wFm(:,i) = chi0wFm(:,i)*chi0w_inv(i2,i2,iwf)
              chi0wFd(:,i) = chi0wFd(:,i)*chi0w_inv(i2,i2,iwf)
              ! Remove the 1 
              chi0wFm(i,i) = chi0wFm(i,i) - 1d0
              chi0wFd(i,i) = chi0wFd(i,i) - 1d0
           enddo
        enddo



        gammawm = 0.d0
        gammawd = 0.d0

        if (external_threelegs) then ! read gamma^w from external file
          call read_threeleg(gammawm,'magn',iwb)
          call read_threeleg(gammawd,'dens',iwb)

        else ! calculate gamma^w=chi0^w.F^w
          do i1=1,ndim2
             do dum=0,2*iwfmax_small-1
                i = i1+dum*ndim2 ! Compound index (i1,iwf)
                gammawm(i1,:) = gammawm(i1,:)+chi0wFm(i,:)
                gammawd(i1,:) = gammawd(i1,:)+chi0wFd(i,:)
             enddo
          enddo
        end if

        if (iq.eq.1 .and. do_chi .and. .not. external_chi_loc) then
           ! Sum up from the right (chi_loc_dens and chi_loc_ magn have already been initialized to zero)
            call calc_chi_qw(chi_loc_dens(:,:,iwb),gammawd,chi0w)
            call calc_chi_qw(chi_loc_magn(:,:,iwb),gammawm,chi0w)
        end if



        ! Add the identity to gamma^w (intermediate quantity used in several equations)
        oneplusgammawm = gammawm
        oneplusgammawd = gammawd
        do i1=1,ndim2
           do dum=0,2*iwfmax_small-1
              i = i1+dum*ndim2 ! Compound index (i1,iwf)
              oneplusgammawm(i1,i) = oneplusgammawm(i1,i) + 1d0
              oneplusgammawd(i1,i) = oneplusgammawd(i1,i) + 1d0
           enddo
        enddo

     endif !update local quantities

     
     call cpu_time(tfinish) ! TIME: local
     timings(1) = timings(1) + tfinish - tstart
     tstart = tfinish ! TIME: eom
   
     ! Calculate static quantities
     if (do_eom .and. iwb .eq. 0 .and. (iq .eq. 1 .or. do_vq)) call calc_eom_static(kq_ind,iq,v,sigma_dmft,sigma_hf)
     ! Calculate local quantities
     if (do_eom .and. iq .eq. 1) call calc_eom_dmft(gammawd,gammawm,iwb,sigma_dmft)
     call cpu_time(tfinish) ! TIME: eom
     timings(6) = timings(6) + tfinish - tstart
     tstart = tfinish ! TIME: chi0 

     if (nonlocal) then
      ! Calculate chi0q in the big box (which has the same size as the small box when do_chi = .false.)
      chi0q=0.d0
      do iwf=iwstart,iwstop
         ! compute k-summed (but still q-dependent) bubble chi0(i1,i2):
         chi0 = 0
         do ik=1,nkp
            ikq = kq_ind(ik,iq) !Index of G(k+q)
            call accumulate_chi0(ik, ikq, iwf, iwb, chi0) ! This subroutine accumulates chi0 (over k)
         enddo
         chi0 = chi0/dble(nkp)
         if (do_chi) bubble_nl(:,:,iqw) = bubble_nl(:,:,iqw) + chi0 ! Accumulate chi0 over iwf
         if (-iwfmax_small .le. iwf .and. iwf .le. iwfmax_small-1) chi0q(:,:,iwf) = chi0 ! Store in chi0q (small iwf box)
      end do
      call cpu_time(tfinish) ! TIME: chi0
      timings(2) = timings(2) + tfinish - tstart
      tstart = tfinish ! TIME: chi
      if (do_chi) then
         ! Get the non-local bubble
         bubble_nl(:,:,iqw) = bubble_nl(:,:,iqw)/(beta**2)  - bubble_loc_tmp
         ! Add gamma^w.chi0^{nl,q} to the q dependent susceptibility
         chi0nl = chi0q - chi0w ! Temporary array with the iwfmax_small dimension to the right
         call calc_chi_qw(chi_qw_dens(:,:,iqw),gammawd,chi0nl)
         call calc_chi_qw(chi_qw_magn(:,:,iqw),gammawm,chi0nl)
      end if
      call cpu_time(tfinish) ! TIME: chi
      timings(7) = timings(7) + tfinish - tstart
      tstart = tfinish ! TIME: gamma

      ! Now we construct eta^q (in the small box..)
      allocate(bigwork_dens(maxdim,maxdim))
      allocate(bigwork_magn(maxdim,maxdim))
      bigwork_dens = 0.d0
      bigwork_magn = 0.d0
      if (do_eom) gammaqd = 0.d0
      ! Loop over the left fermionic matsubara frequency (array index)
      do dum1= 0,2*iwfmax_small-1
         iwf = dum1 - iwfmax_small ! left fermionic matsubara frequency
         offset = dum1*ndim2 ! fermionic matsubara offset
         !get horizontal slice of chi_loc for matmul (interm2):
         chi0wFm_slice = chi0wFm(offset+1:offset+ndim2,:)
         chi0wFd_slice = chi0wFd(offset+1:offset+ndim2,:)

         ! compute intermediate quantity (chi0*chi0_loc_inv - 1) and store it in smallwork (chi0w_inv is diagonal)
         do i2=1,ndim2
            smallwork(:,i2) = chi0q(:,i2,iwf)*chi0w_inv(i2,i2,iwf) ! Nb: chi0w_inv is here assumed to be diagonal
            smallwork(i2,i2) = smallwork(i2,i2) - 1d0
         enddo

         ! compute intermediate quantity chi0^{nl,q}.F^w = (chi0^q-chi0^w).F^w = (chi0*[chi0^w]^{-1}-1)*(chi0^w.F^w):
         alpha = 1.d0
         delta = 0.d0
         call zgemm('n', 'n', ndim2, maxdim, ndim2, alpha, smallwork, ndim2, chi0wFd_slice, ndim2, delta, rectanglework, ndim2)
         ! Store in bigwork_dens (Nb: with a minus sign)
         bigwork_dens(offset+1:offset+ndim2,:) = -rectanglework
         ! Accumulate over iwf in gamma^q_d
         if (do_eom) gammaqd = gammaqd + rectanglework

         !same for the magn channel: chi0^{nl,q}.F^w
         alpha = 1.d0
         delta = 0.d0
         call zgemm('n', 'n', ndim2, maxdim, ndim2, alpha, smallwork, ndim2, chi0wFm_slice, ndim2, delta, rectanglework, ndim2)
         ! Store in bigwork_magn (Nb: with a minus sign)
         bigwork_magn(offset+1:offset+ndim2,:) = -rectanglework

         ! compute part containing nonlocal interaction v (only in density channel)
         if (do_vq) then
            ! 2*beta^(-2)*chi0.v
            smallwork = 2d0/(beta**2)*MATMUL(chi0q(:,:,iwf),v)
            ! beta^(-2)*(chi0*v)*(1+gamma^w_d)
            alpha = 1.d0
            delta = 0.d0
            call zgemm('n','n',ndim2,maxdim,ndim2,alpha,smallwork,ndim2,oneplusgammawd,ndim2,delta,rectanglework,ndim2)
            ! add it to our big work array (Nb: with a minus sign):
            bigwork_dens(offset+1:offset+ndim2,:) = bigwork_dens(offset+1:offset+ndim2,:) - rectanglework
         endif
      enddo ! dum1
      ! TIME: gamma
      call cpu_time(tfinish)
      timings(3) = timings(3) + tfinish - tstart
      tstart = tfinish ! TIME: inv

      ! We need to add the identity before the inversion: 
      ! bigwork_dens = [1 - chi0^{nl,q}.F^w_d - 2*beta^{-2}*chi0^q.v^q.(1 + gamma^w)]
      ! bigwork_magn = [1 - chi0^{nl,q}.F^w_m ]
      do i=1,maxdim
         bigwork_dens(i,i) = bigwork_dens(i,i) + 1d0
         bigwork_magn(i,i) = bigwork_magn(i,i) + 1d0
      enddo

      call inverse_matrix(bigwork_dens)
      call inverse_matrix(bigwork_magn)

      ! We need to subtract the identity before the multiplication from the left with (1 + gamma^w): 
      do i=1,maxdim
         bigwork_dens(i,i) = bigwork_dens(i,i) - 1d0
         bigwork_magn(i,i) = bigwork_magn(i,i) - 1d0
      enddo
      ! TIME: inv
      call cpu_time(tfinish)
      timings(4) = timings(4) + tfinish - tstart
      tstart = tfinish ! TIME: eta

      ! Multiply with (1 + gamma^w) from the left
      etaqd = 0.d0
      etaqm = 0.d0
      alpha = 1.d0
      delta = 0.d0
      call zgemm('n', 'n', ndim2, maxdim, maxdim, alpha, oneplusgammawd &
              , ndim2, bigwork_dens, maxdim, delta, etaqd, ndim2)
      call zgemm('n', 'n', ndim2, maxdim, maxdim, alpha, oneplusgammawm &
                , ndim2, bigwork_magn, maxdim, delta, etaqm, ndim2)

      ! Now we are done with the bigwork arrays
      deallocate(bigwork_dens, bigwork_magn)
      ! TIME: eta
      call cpu_time(tfinish)
      timings(5) = timings(5) + tfinish - tstart
      tstart = tfinish ! TIME: chi or EOM

      if (do_chi) then
         ! Add eta^q.chi0^q to the q dependent susceptibility
         call calc_chi_qw(chi_qw_dens(:,:,iqw),etaqd,chi0q)
         call calc_chi_qw(chi_qw_magn(:,:,iqw),etaqm,chi0q)
         call cpu_time(tfinish)  ! TIME: chi
         timings(7) = timings(7) + tfinish - tstart
         tstart = tfinish ! restart the timer
      end if
      if (do_eom) then
         !equation of motion
         call calc_eom_dynamic(etaqd,etaqm,gammawd,gammaqd,kq_ind,iwb,iq,v,sigma_nl) 
         ! TIME: eom
         call cpu_time(tfinish)
         timings(6) = timings(6) + tfinish - tstart
         tstart = tfinish ! restart the timer
      end if
     endif ! non-local
     call cpu_time(finish)

     !Output the calculation progress
     if (ounit .gt. 0 .and. .not. (verbose .and. (index(verbstr,"Noprogress") .ne. 0))) then
      if (mod(iqw-qwstart,max((qwstop-qwstart)/10,2)) .eq. 0) then
         write(ounit,'(1x,"Core:",I5,"  Completed qw-point: ",I7," (from ",I7," to ",I7,")  Time per point: ",F8.4)') &
               mpi_wrank, iqw, qwstart, qwstop, finish-start
         call flush(ounit)
      endif
     endif
  enddo !iqw
  if (ounit .ge. 1) then
    if (.not. (verbose .and. (index(verbstr,"Noprogress") .ne. 0))) then
      write(ounit,'(1x,"-------------------------------------------------------------------------------------------")')
      write(ounit,'(1x)')
    endif
  endif

  ! Gather the timing
  if (verbose .and. (index(verbstr,"Time") .ne. 0)) then
#ifdef MPI
     if (ounit .ge. 1) then
       write(ounit,'(1x,"TIME: Wall time per qw-point: (Rank ",i6,")")') mpi_wrank
       write(ounit,'(1x,"       Local,       Chi0,     Gammas,  Inversion,   Eta rest,        EOM,        Chi")') 
       write(ounit,'(1x,7f12.5)') timings/(qwstop-qwstart)
     endif
     call MPI_allreduce(MPI_IN_PLACE,timings,size(timings), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
     if (ounit .ge. 1) then
       write(ounit,'(1x,"TIME: Wall time per qw-point:")')
       write(ounit,'(1x,"       Local,       Chi0,     Gammas,  Inversion,   Eta rest,        EOM,        Chi")') 
       write(ounit,'(1x,7f12.5)') timings/(nqp*(2*iwbmax_small+1))
       write(ounit,'(1x)')
       call flush(ounit)
     endif
  endif

  if (do_chi) then
     deallocate(chi0w,chi0nl,bubble_loc_tmp)
  endif
  if (do_eom) then
     deallocate(gammaqd)
  endif
  deallocate(etaqd,etaqm)
  deallocate(u,u_tilde,chi0w_inv,chi0q,smallwork,gammawd,gammawm,stat=er)
  if (er .ne. 0) call mpi_stop('Error in the deallocation:',er=er)
  deallocate(oneplusgammawd,oneplusgammawm)
  deallocate(chi0wFm,chi0wFd,v,rectanglework)
  deallocate(chi0wFd_slice,chi0wFm_slice)
  deallocate(chi0)

  ! MPI reduction and output
  if (do_eom) then
     allocate(sigma_sum(ndim, ndim, -iwfmax_small:iwfmax_small-1, nkp),sigma_sum_hf(ndim,ndim,nkp))
     allocate(sigma_sum_dmft(ndim, ndim, -iwfmax_small:iwfmax_small-1))
#ifdef MPI
     call MPI_reduce(sigma_nl, sigma_sum, ndim*ndim*2*iwfmax_small*nkp, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, ierr)
     call MPI_reduce(sigma_hf,sigma_sum_hf,ndim*ndim*nkp,MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, ierr)
     call MPI_reduce(sigma_dmft,sigma_sum_dmft,ndim*ndim*2*iwfmax_small,MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, ierr)
#else
     sigma_sum = sigma
     sigma_sum_hf = sigma_hf
     sigma_sum_dmft = sigma_dmft
#endif
     if (mpi_wrank .eq. master) then
       allocate(gloc(-iwmax:iwmax-1,ndim,ndim))
       allocate(sigma_loc(ndim, ndim, -iwfmax_small:iwfmax_small-1))
       gloc = 0
       sigma_loc = 0
       if (verbose .and. (index(verbstr,"Test") .ne. 0)) then
         write(ounit,'(1x,"Orbital trace of the self-energy at the first negative and positive matsubara: (Re,Im),(Re,Im)")') 
         write(ounit,'(1x,"Tr[QMC Self-energy]:         ",4f24.11)') sum((/(siw(-1,i),i=1,ndim)/)),sum((/(siw(0,i),i=1,ndim)/))
         write(ounit,'(1x,"Tr[Local Self-energy]:       ",4f24.11)') sum( (/ (sigma_sum_dmft(i,i,-1),i=1,ndim) /) ), &
                                                                    sum( (/ (sigma_sum_dmft(i,i,0),i=1,ndim) /) )
         write(ounit,'(1x,"Tr[Non-local Self-energy]:   ",4f24.11)') sum( (/ (sigma_sum(i,i,-1,1),i=1,ndim) /) ), &
                                                                    sum( (/ (sigma_sum(i,i,0,1),i=1,ndim) /) ) 
       endif
       call add_siw_dmft(sigma_sum)  !add the dmft-selfenergy
       call get_sigma_g_loc(sigma_sum, sigma_loc, gloc) ! calculate the k-summed dga selfenergy and k-summed dga(dmft) greens-function
       if (verbose .and. (index(verbstr,"Test") .ne. 0)) then
         write(ounit,'(1x,"Tr[Total Self-energy]:       ",4f24.11)') sum( (/ (sigma_sum(i,i,-1,1),i=1,ndim) /) ), &
                                                                    sum( (/ (sigma_sum(i,i,0,1),i=1,ndim) /) )
         write(ounit,'(1x,"Tr[Total local Self-energy]: ",4f24.11)') sum( (/ (sigma_loc(i,i,-1),i=1,ndim) /) ), &
                                                                     sum( (/ (sigma_loc(i,i,0),i=1,ndim) /) )
         write(ounit,'(1x,"Tr[QMC Greens function]:     ",4f24.11)') sum((/(giw(-1,i),i=1,ndim)/)),sum((/(giw(0,i),i=1,ndim)/)) 
         write(ounit,'(1x,"Tr[Local Greens function]:   ",4f24.11)') sum((/(gloc(-1,i,i),i=1,ndim)/)),&
                                                                    sum((/(gloc(0,i,i),i=1,ndim)/))
         write(ounit,'(1x)')
         call flush(ounit)
       endif

       if (text_output) then
         call output_eom(sigma_sum, sigma_sum_dmft, sigma_sum_hf, sigma_loc, gloc, nonlocal)
       endif
       if (nonlocal) then
         call output_eom_h5(output_filename,sigma_sum,sigma_sum_hf,sigma_loc,sigma_sum_dmft)
         call get_ndga(sigma_sum) ! calculate the k-dependent and k-summed dga occupation
         call output_occ_h5(output_filename)
       endif
       deallocate(gloc,sigma_loc)
     end if
     deallocate(sigma_nl, sigma_sum, sigma_sum_dmft, sigma_sum_hf)
  end if
  deallocate(giw)

  if (do_chi) then
#ifdef MPI
    if (mpi_wrank.eq.master) then
       call MPI_reduce(MPI_IN_PLACE,bubble_loc,ndim2*ndim2*(2*iwbmax_small+1),&
                       MPI_DOUBLE_COMPLEX,MPI_SUM,master,MPI_COMM_WORLD,ierr)
      if (.not. external_chi_loc) then
       call MPI_reduce(MPI_IN_PLACE,chi_loc_dens,ndim2*ndim2*(2*iwbmax_small+1),&
                       MPI_DOUBLE_COMPLEX,MPI_SUM,master,MPI_COMM_WORLD,ierr)
       call MPI_reduce(MPI_IN_PLACE,chi_loc_magn,ndim2*ndim2*(2*iwbmax_small+1),&
                       MPI_DOUBLE_COMPLEX,MPI_SUM,master,MPI_COMM_WORLD,ierr)
      endif
    else
       call MPI_reduce(bubble_loc,bubble_loc,ndim2*ndim2*(2*iwbmax_small+1),&
                       MPI_DOUBLE_COMPLEX,MPI_SUM,master,MPI_COMM_WORLD,ierr)
      if (.not. external_chi_loc) then
       call MPI_reduce(chi_loc_dens,chi_loc_dens,ndim2*ndim2*(2*iwbmax_small+1),&
                       MPI_DOUBLE_COMPLEX,MPI_SUM,master,MPI_COMM_WORLD,ierr)
       call MPI_reduce(chi_loc_magn,chi_loc_magn,ndim2*ndim2*(2*iwbmax_small+1),&
                       MPI_DOUBLE_COMPLEX,MPI_SUM,master,MPI_COMM_WORLD,ierr)
      endif
    end if
#endif
    ! First the local bubble
    if (mpi_wrank .eq. master) then
      if (.not. external_chi_loc) then
         ! Add the local bubble to chi_loc! (Since it is included in the external chi_loc)
         chi_loc_dens = chi_loc_dens + bubble_loc
         chi_loc_magn = chi_loc_magn + bubble_loc
      endif
      if (verbose .and. (index(verbstr,"Test") .ne. 0)) then
         write(ounit,'(1x,"Orbital sum of the local susceptibilities at w = 0:")') 
         write(ounit,'(1x,"Sum Chi_0^w:            ",2f12.7)') & 
               sum((/((bubble_loc(i+(i-1)*ndim,j+(j-1)*ndim,0),i=1,ndim),j=1,ndim)/))
         write(ounit,'(1x,"Sum Chi_d^w:            ",2f12.7)') & 
               sum((/((chi_loc_dens(i+(i-1)*ndim,j+(j-1)*ndim,0),i=1,ndim),j=1,ndim)/))
         write(ounit,'(1x,"Sum Chi_m^w:            ",2f12.7)') & 
               sum((/((chi_loc_magn(i+(i-1)*ndim,j+(j-1)*ndim,0),i=1,ndim),j=1,ndim)/))
         write(ounit,'(1x)')
         call flush(ounit)
      endif

      ! text file is always reduced output
      if (text_output) then
        call output_chi_loc(bubble_loc,'chi_bubble_loc.dat')
        call output_chi_loc(chi_loc_dens,'chi_dens_loc.dat')
        call output_chi_loc(chi_loc_magn,'chi_magn_loc.dat')
      endif

      if (susc_full_output) then
        call output_chi_loc_full_h5(output_filename,'bubble_loc',bubble_loc)
        call output_chi_loc_full_h5(output_filename,'magn',chi_loc_magn)
        call output_chi_loc_full_h5(output_filename,'dens',chi_loc_dens)
      else
        call output_chi_loc_reduced_h5(output_filename,'bubble_loc',bubble_loc)
        call output_chi_loc_reduced_h5(output_filename,'magn',chi_loc_magn)
        call output_chi_loc_reduced_h5(output_filename,'dens',chi_loc_dens)
      endif
    end if
    deallocate(bubble_loc)

    if (nonlocal) then
      if (mpi_wrank .eq. master) then
        allocate(chi_qw_full(ndim2,ndim2,nqp*(2*iwbmax_small+1)))
      else
        allocate(chi_qw_full(1,1,1))
      end if
      chi_qw_full=0.d0

      ! the purely non-local bubble
#ifdef MPI
      call MPI_gatherv(bubble_nl,(qwstop-qwstart+1)*ndim2**2,MPI_DOUBLE_COMPLEX,&
                     chi_qw_full,rct,disp,MPI_DOUBLE_COMPLEX,master,MPI_COMM_WORLD,ierr)
#else
      chi_qw_full = bubble_nl
#endif
      if (mpi_wrank .eq. master) then
        ! call output_chi_qw(chi_qw_full,qw,'bubble_nl.dat')
        if (verbose .and. (index(verbstr,"Test") .ne. 0)) then
           write(ounit,'(1x,"Orbital sum of non-local Chi0 at w = 0 and first q-point (nl = purely non-local):")') 
           write(ounit,'(1x,"Sum Chi_0^q  - Chi_0^w: ",2f12.7)') &
                 sum((/((chi_qw_full(i+(i-1)*ndim,j+(j-1)*ndim,iwbmax_small*nqp+1),i=1,ndim),j=1,ndim)/))
           call flush(ounit)
        endif

        !call output_chi_qw(chi_qw_full,qw,'chi_qw_dens.dat')
        if (text_output) then
          call output_chi_qw(chi_qw_full,qw,'chi_bubble_nl.dat')
        endif

        if (susc_full_output) then
          if (q_vol) then
            call output_chi_qw_full_h5(output_filename,'bubble_nl',chi_qw_full)
          else if (q_path_susc) then
            call output_chi_qpath_full_h5(output_filename,'bubble_nl',chi_qw_full)
          end if
        else
          if (q_vol) then
            call output_chi_qw_reduced_h5(output_filename,'bubble_nl',chi_qw_full)
          else if (q_path_susc) then
            call output_chi_qpath_reduced_h5(output_filename,'bubble_nl',chi_qw_full)
          end if
        endif

      end if

      !--------------------
      !   Density channel
      !--------------------
      if (verbose .and. ((index(verbstr,"Test") .ne. 0) .or. (index(verbstr,"Extra") .ne. 0))) then
         ! the non-local chi^q_d part: Chi^q_d - Chi^q_0
#ifdef MPI
         call MPI_gatherv(chi_qw_dens,(qwstop-qwstart+1)*ndim2**2,MPI_DOUBLE_COMPLEX,&
                        chi_qw_full,rct,disp,MPI_DOUBLE_COMPLEX,master,MPI_COMM_WORLD,ierr)
#else
         chi_qw_full = chi_qw_dens
#endif
         if (mpi_wrank .eq. master) then
           if (verbose .and. (index(verbstr,"Test") .ne. 0)) then
              write(ounit,'(1x,"Sum Chi_d^nl - Chi_0^q: ",2f12.7)') &
                    sum((/((chi_qw_full(i+(i-1)*ndim,j+(j-1)*ndim,iwbmax_small*nqp+1),i=1,ndim),j=1,ndim)/))
              call flush(ounit)
           endif
           ! Print to file
           if (index(verbstr,"Extra") .ne. 0) then 
             if (susc_full_output) then
               if (q_vol) then
                 call output_chi_qw_full_h5(output_filename,'dens-nl',chi_qw_full)
               else if (q_path_susc) then
                 call output_chi_qpath_full_h5(output_filename,'dens-nl',chi_qw_full)
               endif
             else
               if (q_vol) then
                 call output_chi_qw_reduced_h5(output_filename,'dens-nl',chi_qw_full)
               else if (q_path_susc) then
                 call output_chi_qpath_reduced_h5(output_filename,'dens-nl',chi_qw_full)
               endif
             endif
           endif
         endif
      endif
      ! Add the purely non-local bubble
      chi_qw_dens = chi_qw_dens + bubble_nl
#ifdef MPI
      call MPI_gatherv(chi_qw_dens,(qwstop-qwstart+1)*ndim2**2,MPI_DOUBLE_COMPLEX,&
                     chi_qw_full,rct,disp,MPI_DOUBLE_COMPLEX,master,MPI_COMM_WORLD,ierr)
#else
      chi_qw_full = chi_qw_dens
#endif
      if (mpi_wrank .eq. master) then
        ! Add the local susceptibility
        do iw=1,2*iwbmax_small+1
          do iq=1,nqp
            chi_qw_full(:,:,iq+nqp*(iw-1))=chi_qw_full(:,:,iq+nqp*(iw-1))+chi_loc_dens(:,:,iw-iwbmax_small-1)
          end do
        end do
        if (verbose .and. (index(verbstr,"Test") .ne. 0)) then
           write(ounit,'(1x,"Sum Chi_d^q:            ",2f12.7)') &
                 sum((/((chi_qw_full(i+(i-1)*ndim,j+(j-1)*ndim,iwbmax_small*nqp+1),i=1,ndim),j=1,ndim)/))
           call flush(ounit)
        endif

        ! text file is always reduced output
        if (text_output) then
          call output_chi_qw(chi_qw_full,qw,'chi_dens_nl.dat')
        endif
        if (susc_full_output) then
          if (q_vol) then
            call output_chi_qw_full_h5(output_filename,'dens',chi_qw_full)
          else if (q_path_susc) then
            call output_chi_qpath_full_h5(output_filename,'dens',chi_qw_full)
          end if
        else
          if (q_vol) then
            call output_chi_qw_reduced_h5(output_filename,'dens',chi_qw_full)
          else if (q_path_susc) then
            call output_chi_qpath_reduced_h5(output_filename,'dens',chi_qw_full)
          end if
        endif

      end if
      deallocate(chi_loc_dens)

      !--------------------
      !   Magnetic channel
      !--------------------
      if (verbose .and. ((index(verbstr,"Test") .ne. 0) .or. (index(verbstr,"Extra") .ne. 0))) then
         ! the non-local chi^q_m part: Chi^q_m - Chi^q_0
#ifdef MPI
         call MPI_gatherv(chi_qw_magn,(qwstop-qwstart+1)*ndim2**2,MPI_DOUBLE_COMPLEX,&
                     chi_qw_full,rct,disp,MPI_DOUBLE_COMPLEX,master,MPI_COMM_WORLD,ierr)
#else
         chi_qw_full = chi_qw_magn
#endif
         if (mpi_wrank .eq. master) then
           if (verbose .and. (index(verbstr,"Test") .ne. 0)) then
              write(ounit,'(1x,"Sum Chi_m^nl - Chi_0^q:  ",2f12.7)') &
                    sum((/((chi_qw_full(i+(i-1)*ndim,j+(j-1)*ndim,iwbmax_small*nqp+1),i=1,ndim),j=1,ndim)/))
              call flush(ounit)
           endif
           ! Print to file
           if (index(verbstr,"Extra") .ne. 0) then 
             if (susc_full_output) then
               if (q_vol) then
                 call output_chi_qw_full_h5(output_filename,'magn-nl',chi_qw_full)
               else if (q_path_susc) then
                 call output_chi_qpath_full_h5(output_filename,'magn-nl',chi_qw_full)
               endif
             else
               if (q_vol) then
                 call output_chi_qw_reduced_h5(output_filename,'magn-nl',chi_qw_full)
               else if (q_path_susc) then
                 call output_chi_qpath_reduced_h5(output_filename,'magn-nl',chi_qw_full)
               endif
             endif
           endif
         endif
      endif
      ! Add the purely non-local bubble
      chi_qw_magn = chi_qw_magn + bubble_nl
#ifdef MPI
      call MPI_gatherv(chi_qw_magn,(qwstop-qwstart+1)*ndim2**2,MPI_DOUBLE_COMPLEX,&
                     chi_qw_full,rct,disp,MPI_DOUBLE_COMPLEX,master,MPI_COMM_WORLD,ierr)
#else
      chi_qw_full = chi_qw_magn
#endif
      if (mpi_wrank .eq. master) then
        ! Add the local susceptibility
        do iw=1,2*iwbmax_small+1
          do iq=1,nqp
            chi_qw_full(:,:,iq+nqp*(iw-1))=chi_qw_full(:,:,iq+nqp*(iw-1))+chi_loc_magn(:,:,iw-iwbmax_small-1)
          end do
        end do
        if (verbose .and. (index(verbstr,"Test") .ne. 0)) then
           write(ounit,'(1x,"Sum Chi_m^q:            ",2f12.7)') &
                 sum((/((chi_qw_full(i+(i-1)*ndim,j+(j-1)*ndim,iwbmax_small*nqp+1),i=1,ndim),j=1,ndim)/))
           call flush(ounit)
        endif

        ! text file is always reduced output
        if (text_output) then
          call output_chi_qw(chi_qw_full,qw,'chi_magn_nl.dat')
        endif

        if (susc_full_output) then
          if (q_vol) then
            call output_chi_qw_full_h5(output_filename,'magn',chi_qw_full)
          else if (q_path_susc) then
            call output_chi_qpath_full_h5(output_filename,'magn',chi_qw_full)
          end if
        else
          if (q_vol) then
            call output_chi_qw_reduced_h5(output_filename,'magn',chi_qw_full)
          else if (q_path_susc) then
            call output_chi_qpath_reduced_h5(output_filename,'magn',chi_qw_full)
          end if
        endif
      end if
      deallocate(chi_loc_magn)

      deallocate(chi_qw_full)
    endif
    deallocate(chi_qw_magn)
    deallocate(chi_qw_dens)
    deallocate(bubble_nl)
  end if


! Output
  deallocate(iw_data,iwb_data,siw,k_data,q_data,kq_ind,qw)
  if (ounit .ge. 1) then
      write(ounit,'(1x)')
      write(ounit,'(1x,"End of Program")')
  endif
  close(ounit)
  call mpi_close()
end program main
