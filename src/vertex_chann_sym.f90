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

program symmetrize_vertex
!====================================
  use parameters_module
  use aux
  use hdf5
  use hdf5_module
  implicit none

  integer(hid_t)                  :: file_id, new_file_id
  integer(hid_t)                  :: grp_id, g4iw_id, g4err_id
  integer(hid_t)                  :: dset_dens_id, dset_magn_id, dset_err_id
  integer(hsize_t), dimension(3)  :: g4iw_dims, g4iw_maxdims, g4err_dims, g4iw_dims_t, g4err_dims_t
  integer                         :: nmembers, imembers, itype
  integer, parameter              :: rank = 2, rank_iw = 1
  character(len=80)               :: grpname, name_buffer, name_buffer_value, name_buffer_error
  character(len=80)               :: name_buffer2, name_buffer3
  
  double precision, allocatable   :: g4iw_r(:,:,:), g4iw_i(:,:,:), g4err(:,:,:), diff_r(:,:), diff_i(:,:)
  double precision, allocatable   :: g4iw_r_t(:,:,:), g4iw_i_t(:,:,:), g4err_t(:,:,:)
  integer                         :: ind, b1, s1, b2, s2, b3, s3, b4, s4, ineq
  integer, allocatable            :: Nbands(:)
  integer                         :: iwb, iwf, iwf1, iwf2
  logical                         :: su2_only
  logical, allocatable            :: create_comp(:,:)
  integer                         :: ichannel, ntot, ind_band, icount
  integer, allocatable            :: ind_band_list(:)
  character(len=150), allocatable :: filename_vertex_ineq(:)
  character(len=1)                :: arg_sym

  real(kind=8) :: start, finish
!================================================================

! alternatively - read input
! this is definitely more clearer
  write(*,'(A)',advance='no') 'Number of inequivalent atoms: '
  read(*,*) nineq
  allocate(filename_vertex_ineq(nineq))
  allocate(Nbands(nineq))
  do ineq=1,nineq
    if (ineq .eq. 1) then
      write(*,'(A,I1,A)',advance='no') 'Vertex file : '  ! , ineq, ': '
      read(*,*) filename_vertex_ineq(ineq)
    endif
    filename_vertex_ineq(ineq)=filename_vertex_ineq(1)
    write(*,'(A,I1,A)',advance='no') 'Number of correlated bands for inequivalent atom ', ineq, ': '
    read(*,*) Nbands(ineq)
  enddo
  write(*,'(A)',advance='no') 'Outputfile for symmetrized Vertex: '
  read(*,*) filename_vertex_sym
  write(*,*)
  write(*,'(A)',advance='no') 'SU2 symmetry only (s) or SU2 AND orbital symmetry (o)?: '
  read(*,*) arg_sym
  if (arg_sym .eq. 'o' .or. arg_sym .eq. 'O') then
    su2_only = .false.
  elseif (arg_sym .eq. 's' .or. arg_sym .eq. 'S') then
    su2_only = .true.
  else
    su2_only = .false.
    write(*,*) 'Wrong input - Using only SU2 symmetry.'
  endif


!================================================================
!Define orbital symmetry here:
  ! su2_only = .false. 
  write(*,*) 'Symmetrizing ',(filename_vertex_ineq(ineq),ineq=1,nineq),'>>>>>',filename_vertex_sym
  write(*,*) 'Total number of bands: ',sum(Nbands)
  if(su2_only) then
    write(*,*) 'Using only SU2 symmetry'
  else
    write(*,*) 'Using orbital and SU2 symmetry'
  endif
    
!=================================================================

  call cpu_time(start)
  
  call h5open_f(err)
 
  call create_complex_datatype
 


! loop over number of inequivalent atoms
  do ineq=1,nineq


! open vertex file 'vertex_full.hdf5':
    call h5fopen_f(trim(filename_vertex_ineq(ineq)), h5f_acc_rdonly_f, file_id, err)

    if (ineq .eq. 1) then ! only write the axes once
! create new file for the symmetrised vertex:
      call h5fcreate_f(filename_vertex_sym, h5f_acc_trunc_f, new_file_id, err)
! get fermionic and bosonic Matsubara axes: 
      call read_and_write_axes(file_id,new_file_id)
      !call read_axes(file_id, iwb_array, iwf_array, dspace_iwb_id, dspace_iwf_id, dim_iwb, dim_iwf)
      !call write_axes(new_file_id, iwb_array, iwf_array, dspace_iwb_id, dspace_iwf_id, dim_iwb, dim_iwf)
    endif

! write magn/iwb and dens/iwb in the output file vertex_sym.hdf5: 
    call create_channels(new_file_id, ineq)

!=========================================================================

! allocate quantities that are needed in the loop afterwards: 
    g4iw_dims = (/2*iwfmax,2*iwfmax,2*iwbmax+1/)
    g4err_dims = g4iw_dims

    g4iw_dims_t = (/2*iwbmax+1,2*iwfmax,2*iwfmax/)
    g4err_dims_t = g4iw_dims_t

    allocate(g4iw_r(2*iwfmax,2*iwfmax,2*iwbmax+1)) 
    allocate(g4iw_i(2*iwfmax,2*iwfmax,2*iwbmax+1))
    allocate(g4err(2*iwfmax, 2*iwfmax, 2*iwbmax+1))

    allocate(g4iw_r_t(2*iwbmax+1,2*iwfmax,2*iwfmax)) 
    allocate(g4iw_i_t(2*iwbmax+1,2*iwfmax,2*iwfmax))
    allocate(g4err_t(2*iwbmax+1,2*iwfmax, 2*iwfmax))

    allocate(tmp_r_1(g4iw_dims(1), g4iw_dims(2)), tmp_i_1(g4iw_dims(1), g4iw_dims(2)), tmp_err_1(g4iw_dims(1),g4iw_dims(2)))
    allocate( diff_r(g4iw_dims(1), g4iw_dims(2)), diff_i(g4iw_dims(1), g4iw_dims(2)))

    allocate(create_comp(2,Nbands(ineq)**4))
    create_comp = .true.

    allocate(ind_band_list(Nbands(ineq)**2))

! create dataspace:
    dims = (/g4iw_dims(1), g4iw_dims(2)/)
    call h5screate_simple_f(rank, dims, dspace_id, err)
 
!============================================================================
! iterate over all groups(band-spin combinations) in the old vertex file:
    write(name_buffer2,'("/worm-001/ineq-",I3.3,"/g4iw-worm/")') ineq
    call h5gn_members_f(file_id, trim(name_buffer2), nmembers, err)

    do imembers = 0,nmembers-1
       call h5gget_obj_info_idx_f(file_id, trim(name_buffer2), imembers, name_buffer, itype, err)
       read(name_buffer,'(I5.5)') ind
       write(*,*) imembers, trim(name_buffer)
     
       ! read the current group in the old vertex file:
       ! call h5gopen_f(file_id,trim(nam_buffer2) // trim(name_buffer) // "value" ,grp_id,err)

       ! write(name_buffer_value, '((I5.5),A6)') ind, "/value"
       name_buffer_value = trim(name_buffer2) // trim(name_buffer) // "/value"
       call h5dopen_f(file_id,name_buffer_value,g4iw_id,err)
       ! call h5dopen_f(file_id, name_buffer_value, g4iw_id, err)

       ! write(name_buffer_error, '((I5.5),A6)') ind, "/error"
       name_buffer_error = trim(name_buffer2) // trim(name_buffer) // "/error"
       call h5dopen_f(file_id,name_buffer_error,g4err_id,err)
       ! call h5dopen_f(file_id, name_buffer_error, g4err_id, err)

       call h5dread_f(g4iw_id, type_r_id, g4iw_r_t, g4iw_dims_t, err)
       call h5dread_f(g4iw_id, type_i_id, g4iw_i_t, g4iw_dims_t, err)
       call h5dread_f(g4err_id, h5t_native_double, g4err_t, g4err_dims_t, err)

       call h5dclose_f(g4iw_id, err)
       call h5dclose_f(g4err_id, err)
       ! call h5gclose_f(grp_id, err)

       ! Reason for doing this: old w2d format was bff; new w2d format is ffb
       ! Attention: Fortran is implicitly transposed due to the way fortran / C access memory
       ! transpose back
       do iwb=1,2*iwbmax+1
         do iwf=1,2*iwfmax
           g4iw_r(:,iwf,iwb) = g4iw_r_t(iwb,iwf,:)
           g4iw_i(:,iwf,iwb) = g4iw_i_t(iwb,iwf,:)
           g4err(:,iwf,iwb) = g4err_t(iwb,iwf,:)
         enddo
       enddo

  ! get band and spin indices:
       call index2component(Nbands(ineq), ind, b1, s1, b2, s2, b3, s3, b4, s4)
      
  ! get the channel (magnetic: ichannel=1, density: ichannel=2)
       call get_channel(s1, s2, s3, s4, ichannel)
      
       ! get list of ind_band (components to be written in vertex_sym.hdf5)
       if (su2_only) then

          ! without orbital symmetry (only su2):
          call component2index_band(Nbands(ineq), ind_band, b1, b2, b3, b4)
          ntot = 1 
          ind_band_list(1) = ind_band

       else

          !with orbital symmetry: 
          call get_orb_sym(b1, b2, b3, b4, Nbands(ineq), ntot, ind_band_list)

       endif     

       !divisions needed for the channels:
       g4iw_r = g4iw_r/(2.d0*ntot)
       g4iw_i = g4iw_i/(2.d0*ntot)

       if (ichannel == 1) then  !is this correct??? 
          g4err = g4err/(2.d0*ntot)
       else
          g4err = g4err/(4.d0*ntot)
       endif
      
       ! iterates over all components that need to be written: 
       do icount = 1, ntot

         
          if (create_comp(ichannel, ind_band_list(icount))) then
      
             ! creates new groups and datasets if they are not there yet:
             do iwb = 0, 2*iwbmax
                call create_component(new_file_id, ichannel, iwb, ind_band_list(icount), ineq)
             enddo
             
             create_comp(ichannel, ind_band_list(icount)) = .false.

          endif
          
          !adds the data:
          do iwb = 0, 2*iwbmax

             call add_to_component(new_file_id, ichannel, iwb, ind_band_list(icount), g4iw_r, g4iw_i, g4err, ineq)

          enddo
         
       enddo !icount
   
    enddo ! nmembers
  
 
    call h5sclose_f(dspace_id, err)

    deallocate(g4iw_r, g4iw_i, g4err)
    deallocate(g4iw_r_t, g4iw_i_t, g4err_t)
    deallocate(tmp_r_1, tmp_i_1, tmp_err_1, diff_r, diff_i) 
    deallocate(ind_band_list)
    deallocate(create_comp)

    call h5fclose_f(file_id, err)

  enddo

  call h5fclose_f(new_file_id, err)

  call cpu_time(finish)
  write(*,*)'cpu-time =', finish-start
  

end program symmetrize_vertex
