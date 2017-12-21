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

module kq_tools
  use parameters_module
  implicit none

  interface k_vector
    module procedure k_vector_1, k_vector_3
  end interface k_vector
  interface k_index
    module procedure k_index_1, k_index_3
  end interface k_index

contains

subroutine generate_q_vol(nqpx,nqpy,nqpz,qdata)
  implicit none
  integer :: nqpx,nqpy,nqpz
  integer :: qdata(nqpx*nqpy*nqpz),i,j,k,i1

  i1=0
  do i=0,nqpx-1
    do j=0,nqpy-1
      do k=0,nqpz-1
        i1 = i1+1
        qdata(i1)=k_index(i*nkpx/nqpx,j*nkpy/nqpy,k*nkpz/nqpz)
      enddo
    enddo
  enddo

end subroutine generate_q_vol

subroutine index_kq(ind)
      implicit none
      integer ikp,jkp
      integer :: ind(nkp,nqp)
      ind = 0

      do ikp=1,nkp
        do jkp=1,nqp
          ind(ikp,jkp)=k_minus_q(ikp,q_data(jkp))
        end do
      end do

end subroutine index_kq 

! The following function calculates the index of \vec{k} - \vec{q}.
! It uses only integers
! \vec{k} is associated to (ix,iy,iz)
! \vec{q} is associated to (lx,ly,lz)
! k-space is assumed to have nkpx*nkpy*nkpz points
! q-space is assumed to have nqpx*nqpy*nqpz points,
! where each element of the q-space has to be an element of the k-space.
! subtractions are done in integers, 
! fold-back to BZ is achieved by modulo division.
function k_minus_q(ik,iq)
  implicit none
  integer :: ik,iq,k_minus_q
  integer :: ix,iy,iz,lx,ly,lz

  call k_vector(ik,ix,iy,iz)
  call k_vector(iq,lx,ly,lz)

  k_minus_q=1+mod(nkpz+iz-lz,nkpz) + &
              mod(nkpy+iy-ly,nkpy)*nkpz + &
              mod(nkpx+ix-lx,nkpx)*nkpy*nkpz

end function k_minus_q

function k_index_1(k)
  implicit none
  integer,intent(in) :: k(3)
  integer :: k_index_1

  k_index_1 = 1 + k(3) + k(2)*nkpz + k(1)*nkpy*nkpz

end function k_index_1

function k_index_3(kx,ky,kz)
  implicit none
  integer :: k_index_3,kx,ky,kz

  k_index_3 = 1 + kz + ky*nkpz + kx*nkpy*nkpz

end function k_index_3

subroutine k_vector_1(ik,k)
  implicit none
  integer,intent(in) :: ik
  integer,intent(out) :: k(3)

  k(3)=mod(ik-1,nkpz)
  k(2)=mod((ik-1)/nkpz,nkpy)
  k(1)=(ik-1)/(nkpy*nkpz)
end subroutine k_vector_1

subroutine k_vector_3(ik,kx,ky,kz)
  implicit none
  integer :: ik,kx,ky,kz

  kz=mod(ik-1,nkpz)
  ky=mod((ik-1)/nkpz,nkpy)
  kx=(ik-1)/(nkpy*nkpz)
end subroutine k_vector_3

subroutine qdata_from_file()
  use parameters_module

  implicit none
  integer :: iostatus,iq
  character(100) :: str_tmp
  real(kind=8) :: qx,qy,qz

  iostatus=0
  open(unit=101,file=filename_qdata)
  
  nqp=-1
  do while (iostatus.eq.0)
    read(101,*,iostat=iostatus) str_tmp
    nqp=nqp+1
  end do
  close(101)
  !write(*,*) nqp,' q points'

  allocate(q_data(nqp))
  open(unit=101,file=filename_qdata)
  do iq=1,nqp
    ! We read three real numbers.
    ! If all of them are zero, it is the gamma point.
    ! If one of them is larger or equal to 1, the coordinates are cast to integers
    ! and assumed to be given in integer basis [0,nkpi-1]
    ! If neither of above is true, the coordinates are assumed to lie in the interval [0,1).
    read(101,*) qx,qy,qz
    if (qx .eq. 0 .and. qy .eq. 0 .and. qz .eq. 0) then 
      q_data(iq) = k_index(0,0,0) ! gamma point
    else if (qx .ge. 1 .or. qy .ge. 1 .or. qz .ge. 1) then
      q_data(iq) = k_index(int(qx),int(qy),int(qz)) ! cast to integers
    else
      q_data(iq) = k_index(nint(qx*nkpx),nint(qy*nkpy),nint(qz*nkpz)) ! round to nearest integers
    end if
  end do
  close(101)
  !write(*,*) 'q data',q_data

end subroutine qdata_from_file

end module kq_tools
