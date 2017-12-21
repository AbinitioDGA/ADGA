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

module mpi_org

  use parameters_module, only: iwbmax_small,nqp,ndim2
#ifdef MPI
    use mpi
#endif
  implicit none
  public

  integer :: mpi_wrank
  integer :: mpi_wsize
  integer :: master
  integer, allocatable :: rct(:),disp(:)
  integer :: ierr
  integer :: qwstart,qwstop

  contains

  subroutine mpi_initialize()
    implicit none
#ifdef MPI
      call MPI_init(ierr)
      call MPI_comm_rank(MPI_COMM_WORLD,mpi_wrank,ierr)
      call MPI_comm_size(MPI_COMM_WORLD,mpi_wsize,ierr)
      master = 0
      allocate(rct(mpi_wsize),disp(mpi_wsize))
      rct = 0
      disp = 0
#else
      mpi_wrank=0
      mpi_wsize=1
      master=0
#endif
  end subroutine mpi_initialize


  subroutine mpi_distribute()
    use parameters_module, ONLY: ounit, verbose, verbstr
    implicit none
    integer :: i,j
#ifdef MPI
      rct=0  ! Elements per rank
      disp=0 ! Offset, i.e. disp(1) = 0
      qwstop = 0
      do i=1,mpi_wsize -1
        rct(i) = (nqp*(2*iwbmax_small+1) - disp(i))/(mpi_wsize+1-i)
        disp(i+1) = disp(i) + rct(i)
      enddo
      rct(mpi_wsize) = (nqp*(2*iwbmax_small+1) - disp(mpi_wsize))
      ! Set qwstart and qsstop
      qwstart = disp(mpi_wrank+1)+1
      qwstop = disp(mpi_wrank+1)+rct(mpi_wrank+1)

      if (ounit .ge.  1) then
         write(ounit,'(1x,"mpi_distribute: average number of iqw points per rank (floored):",i20)'),rct(1)
      endif

      ! Multiply rct and disp with ndim2**2 for future use in mpi_gatherv
      rct = rct*ndim2**2
      disp = disp*ndim2**2
      if (ounit .ge. 1 .and. (verbose .and. (index(verbstr,"Mpi") .ne. 0))) then
        write(ounit,*) 'mpi_distribute: receive count (iqw*ndim**2)',rct
        write(ounit,*) 'mpi_distribute: corresponding displacements',disp
      end if
#else
      qwstart = 1
      qwstop = nqp*(2*iwbmax_small+1)
#endif
  end subroutine mpi_distribute


  subroutine mpi_close()
    implicit none
#ifdef MPI
      deallocate(rct,disp)
      call MPI_finalize(ierr)
#endif
  end subroutine mpi_close

  subroutine mpi_stop(stopmsg,er,mpi_er)
   use parameters_module, ONLY: ounit
   implicit none
   character(len=*),intent(in)  :: stopmsg      !< message to print
   integer,optional,intent(in)  :: er           !< error flag to print
   integer,optional,intent(in)  :: mpi_er          !< rank of the processor printing this message
#ifdef MPI
   character(len=MPI_MAX_ERROR_STRING) :: str      ! temp string
   integer                             :: strlen      ! length
   integer                             :: mpi_errorflag   ! temp error flag
   character(len=5)                    :: fancynr  ! temp string for the mpi rank
   logical                             :: unitopened ! tells if the file unit is open
#endif
      inquire(unit=ounit,opened=unitopened)
      if (unitopened) call flush(ounit)
#if MPI 
      write(fancynr,'(i5)') mpi_wrank
      if (.not. present(mpi_er)) then
         if (present(er)) then
            if (unitopened) then
               write(ounit,'(1x,"Rank ",a,", STOP in abinitiodga: ",a," er=",i5)') TRIM(ADJUSTL(fancynr)),stopmsg,er
               close(ounit)
            endif
            write(*,'(1x,"Rank ",a,", STOP in abinitiodga: ",a," er=",i5)') TRIM(ADJUSTL(fancynr)),stopmsg,er
         else
            if (unitopened) then
               write(ounit,'(1x,"Rank ",a,", STOP in abinitiodga: ",a)') TRIM(ADJUSTL(fancynr)),stopmsg
               close(ounit)
            endif
            write(*,'(1x,"Rank ",a,", STOP in abinitiodga: ",a)') TRIM(ADJUSTL(fancynr)),stopmsg
         endif
         call sleep(1)
         call MPI_Abort(MPI_COMM_WORLD,-1,mpi_errorflag)
      else
         if (unitopened) then
            if (present(er)) then
               write(ounit,'(1x,"Rank ",a,", STOP in abinitiodga: ",a," er=",i25)') TRIM(ADJUSTL(fancynr)),stopmsg,er
            endif
            write(ounit,'(1x,"Rank ",a,", STOP in abinitiodga: ",a," mpi=",i5)') TRIM(ADJUSTL(fancynr)),stopmsg,mpi_er
            close(ounit)
         endif
         if (present(er)) write(*,'(1x,"Rank ",a,", STOP in abinitiodga: ",a," er=",i25)') TRIM(ADJUSTL(fancynr)),stopmsg,er
         write(*,'(1x,"Rank ",a,", STOP in abinitiodga: ",a," mpi=",i5)') TRIM(ADJUSTL(fancynr)),stopmsg,mpi_er
         call MPI_ERROR_STRING(mpi_er, str, strlen,mpi_errorflag)
         if (mpi_errorflag .ne. 0) then
            write(*,*) "Error in MPI_ERROR_STRING:",mpi_errorflag
         else
            write(*,*) "STOP in abinitiodga:" // str(1:strlen)
         endif
         call sleep(1)
         call MPI_Abort(MPI_COMM_WORLD,mpi_er,mpi_errorflag)
      endif
#else
      if (present(er)) then
         write(*,*) "STOP in abinitiodga: " // stopmsg // " er=",er
      else if (.not. present(mpi_er)) then
         write(*,*) "STOP in abinitiodga: " // stopmsg
      endif
      if (present(mpi_er)) then
         write(*,*) "ERROR: Got an mpi error but abinitiodga is compiled without MPI support!"
         write(*,*) "STOP in abinitiodga: " // stopmsg // " mpi=",mpi_er
      endif
#endif

      stop
   end subroutine mpi_stop

end module mpi_org
