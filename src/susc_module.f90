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

module susc_module
  implicit none

  contains

  subroutine calc_chi_qw(chi_qw,interm3,chi0_sum)
     use parameters_module
     implicit none
     complex(kind=8), intent(inout) :: chi_qw(ndim2,ndim2)
     complex(kind=8), intent(in)    :: interm3(ndim2,maxdim)
     complex(kind=8), intent(in)    :: chi0_sum(ndim2,ndim2,-iwfmax_small:iwfmax_small-1)
     complex(kind=8)                :: chi_tmp(ndim2,ndim2)
     integer :: dum,iwf,i3,i1,i2,i
     ! Initialize
     chi_tmp = 0
     do i2=1,ndim2
        do dum=0,2*iwfmax_small-1
           do i3=1,ndim2
              iwf = dum - iwfmax_small
              i = i3 + dum*ndim2 ! = {i3,iwf} 
              chi_tmp(:,i2)=chi_tmp(:,i2)+interm3(:,i)*chi0_sum(i3,i2,iwf)
           end do
        end do
     end do
     ! Add to chi_qw
     chi_qw=chi_qw + chi_tmp/(beta**2)
     return
  end subroutine calc_chi_qw

! subroutine to output susceptibility
  subroutine output_chi_qw(chi_qw,qw,filename_output)
    use parameters_module
    implicit none
    complex(kind=8),intent(in)  :: chi_qw(ndim2,ndim2,nqp*(2*iwbmax_small+1))
    integer,intent(in)          :: qw(2,nqp*(2*iwbmax+1))
    character(len=*),intent(in) :: filename_output
    integer :: iwb,iq,i,j,i1

    if (ounit .ge. 1 .and. (verbose .and. (index(verbstr,"Output") .ne. 0))) then
      write(ounit,*) 'Output chi_qw'
    endif

    ! open file and write head line
    open(unit=10,file=trim(output_dir)//filename_output)
    ! the file header should contain important information
    write(10,'(A)') '#iwbmax, nqp, ndim, beta, mu'
    write(10,'(I5,2X,I5,2X,I5,2X,E14.7E2,2X,E14.7E2)') iwbmax_small,nqp,ndim,beta,mu
    write(10,'(A)') '#iwb, wb, iq, qx, qy, qz, (real(chi_qw(wb,q,i,i,j,j)), imag(chi_qw(wb,q,i,i,j,j))) [j=i,ndim] [i=1,ndim]'

    ! loop over all entries
    do i1=1,nqp*(2*iwbmax_small+1)
       iq = qw(2,i1)
       iwb = qw(1,i1)
        write(10,'(I8,E14.7E2,I8,9999E18.7E2)') iwb,iwb_data(iwb),iq,k_data(:,q_data(iq)), &
                   ((real(chi_qw((i-1)*ndim+1,(j-1)*ndim+1,i1)), real(chi_qw((i-1)*ndim+1,(j-1)*ndim+1,i1)), j=i,ndim), i=1,ndim)
       ! insert an empty line after each omega block. could be useful for plotting.
       if (mod(i1,nqp).eq.0) then
          write(10,*) ' '
       end if

    end do

    close(10)
  end subroutine output_chi_qw

  subroutine output_chi_loc(chi_w,filename_output)
    use parameters_module
    implicit none
    complex(kind=8),intent(in) :: chi_w(ndim2,ndim2,2*iwbmax_small+1)
    character(len=*),intent(in) :: filename_output
    integer :: iwb,i,j,i1

    if (ounit .ge. 1 .and. (verbose .and. (index(verbstr,"Output") .ne. 0))) then
     write(ounit,*) 'output chi_loc'
    endif

    ! open file and write head line
    open(unit=10,file=trim(output_dir)//filename_output)
    ! the file header should contain important information
    write(10,'(A)') '#iwbmax, nqp, ndim, beta, mu'
    write(10,'(I5,2X,I5,2X,I5,2X,E14.7E2,2X,E14.7E2)') iwbmax_small,0,ndim,beta,mu !nqp=0 can be used for the postprocessing
    write(10,'(A)') '#iwb, wb, (real(chi_loc(wb,i,i,j,j)), imag(chi_loc(wb,i,i,j,j))) [j=i,ndim] [i=1,ndim]'

    ! loop over all entries
    do i1=1,2*iwbmax_small+1
      iwb = i1-iwbmax_small-1
      write(10,'(I8,9999E18.7E2)') iwb,iwb_data(iwb), ((real(chi_w((i-1)*ndim+1,(j-1)*ndim+1,i1)), aimag(chi_w((i-1)*ndim+1,(j-1)*ndim+1,i1)), j=i,ndim), i=1,ndim)
    end do
    close(10)
  end subroutine output_chi_loc

end module
