module lapack_module
  implicit none
  private
  public inverse_matrix
  
  interface inverse_matrix
    module procedure inverse_matrix_d, inverse_matrix_z
  end interface inverse_matrix
  
  contains
  
  subroutine inverse_matrix_z(M)
    implicit none
    double complex, intent (inout) :: M(:,:)
    integer                        :: ndim,ierr,lwork
    integer,        allocatable    :: ipiv(:)
    double complex, allocatable    :: work(:)
    double complex                 :: work_query(1)
  
    ndim = size(M,1)
    allocate(ipiv(ndim))
    call zgetrf(ndim,ndim,M,ndim,ipiv,ierr)
    if (ierr .ne. 0) then
      write(*,*)"ERROR in ZGETRF", ierr ; stop
    end if
    call zgetri(ndim,M,ndim,ipiv,work_query,-1,ierr) ! query for optimal work space
    if (ierr .ne. 0) then
      write(*,*)"ERROR in ZGETRI at workspace query", ierr ; stop
    endif
    lwork = work_query(1)
    allocate(work(lwork))
    call zgetri(ndim,M,ndim,ipiv,work,lwork,ierr)
    if (ierr .ne. 0) then
      write(*,*)"ERROR in ZGETRI", ierr ; stop
    endif
    deallocate(ipiv,work)
  end subroutine inverse_matrix_z
  
  subroutine inverse_matrix_d(M)
    implicit none
    double precision, intent (inout) :: M(:,:)
    integer                          :: ndim,ierr,lwork
    integer,          allocatable    :: ipiv(:)
    double precision, allocatable    :: work(:)
    double precision                 :: work_query(1)
  
    ndim = size(M,1)
    allocate(ipiv(ndim))
    call dgetrf(ndim,ndim,M,ndim,ipiv,ierr)
    if(ierr .ne. 0) then
      write(*,*)"ERROR in DGETRF", ierr ; stop
    end if
    call dgetri(ndim,M,ndim,ipiv,work_query,-1,ierr) ! query for optimal work space
    if (ierr .ne. 0) then
      write(*,*)"ERROR in DGETRI at workspace query", ierr ; stop
    endif
    lwork = work_query(1)
    allocate(work(lwork))
    call dgetri(ndim,M,ndim,ipiv,work,lwork,ierr)
    if(ierr .ne. 0) then
      write(*,*)"ERROR in DGETRI", ierr ; stop
    end if
    deallocate(ipiv,work)
  end subroutine inverse_matrix_d
  
end module lapack_module
