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
! Copyright (C) <2017-2019>
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

module lapack_module
  implicit none
  private
  public inverse_matrix

  interface inverse_matrix
    module procedure inverse_matrix_d, inverse_matrix_z
  end interface inverse_matrix

  contains

  subroutine inverse_matrix_z(M, erstr, ierr)
    implicit none
    double complex, intent (inout)  :: M(:,:)
    integer, intent(out)            :: ierr
    character(len=200), intent(out) :: erstr
    integer                         :: ndim,lwork
    integer,        allocatable     :: ipiv(:)
    double complex, allocatable     :: work(:)
    double complex                  :: work_query(1)

    ndim = size(M,1)
    allocate(ipiv(ndim))
    call zgetrf(ndim,ndim,M,ndim,ipiv,ierr)
    if (ierr .ne. 0) then
      erstr = "ERROR in ZGETRF"
      return
    endif
    call zgetri(ndim,M,ndim,ipiv,work_query,-1,ierr) ! query for optimal work space
    if (ierr .ne. 0) then
      erstr = "ERROR in ZGETRI at workspace query"
      return
    endif
    lwork = work_query(1)
    allocate(work(lwork))
    call zgetri(ndim,M,ndim,ipiv,work,lwork,ierr)
    if (ierr .ne. 0) then
      erstr = "ERROR in ZGETRI"
      return
    endif
    deallocate(ipiv,work)
  end subroutine inverse_matrix_z

  subroutine inverse_matrix_d(M, erstr, ierr)
    implicit none
    double precision, intent (inout) :: M(:,:)
    integer, intent(out)             :: ierr
    character(len=200), intent(out)  :: erstr
    integer                          :: ndim,lwork
    integer,          allocatable    :: ipiv(:)
    double precision, allocatable    :: work(:)
    double precision                 :: work_query(1)

    ndim = size(M,1)
    allocate(ipiv(ndim))
    call dgetrf(ndim,ndim,M,ndim,ipiv,ierr)
    if(ierr .ne. 0) then
      erstr = "ERROR in DGETRF"
      return
    end if
    call dgetri(ndim,M,ndim,ipiv,work_query,-1,ierr) ! query for optimal work space
    if (ierr .ne. 0) then
      erstr = "ERROR in DGETRI at workspace query"
      return
    endif
    lwork = work_query(1)
    allocate(work(lwork))
    call dgetri(ndim,M,ndim,ipiv,work,lwork,ierr)
    if(ierr .ne. 0) then
      erstr = "ERROR in DGETRI"
      return
    end if
    deallocate(ipiv,work)
  end subroutine inverse_matrix_d

end module lapack_module
