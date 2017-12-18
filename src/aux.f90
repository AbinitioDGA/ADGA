module aux
  implicit none

  contains

!==============================================================================
  subroutine get_channel(s1, s2, s3, s4, ichannel)
    implicit none

    integer, intent(in) :: s1, s2, s3, s4
    integer, intent(out) :: ichannel

    if ((s1==s4) .and. (s2==s3) .and. (s1 .ne. s2)) then
       ichannel = 1
    else
       ichannel = 2
    endif

  end subroutine get_channel
!================================================================================


!================================================================================
  subroutine get_orb_sym(b1, b2, b3, b4, nbands, ntot, ind_band_list)
    implicit none
    
    integer, intent(in) :: b1, b2, b3, b4, nbands
    integer,intent(out) :: ntot
    integer,intent(out) :: ind_band_list(:)
    integer :: i, j, icount, ind_band

    if ((b1==b2) .and. (b2==b3) .and. (b3==b4)) then
       ntot = nbands
       
       do i=1,nbands
          call component2index_band(nbands, ind_band, i, i, i, i)
          ind_band_list(i) = ind_band
       enddo

    else
       ntot = nbands**2-nbands
       icount = 1
       do i=1,nbands
          do j=1,nbands
             if (i .ne. j) then
                ind_band = -1
          
                if((b1==b2) .and. (b3==b4) .and. (b1 .ne. b3)) then
                   call component2index_band(nbands, ind_band, i, i, j, j)
                   
                else if((b1==b3) .and. (b2==b4) .and. (b1 .ne. b2)) then
                   call component2index_band(nbands, ind_band, i, j, i, j)

                else if((b1==b4) .and. (b2==b3) .and. (b1 .ne. b2)) then
                   call component2index_band(nbands, ind_band, i, j, j, i)

                endif

                ind_band_list(icount) = ind_band
                icount = icount+1
             endif
          enddo
       enddo

    endif

  end subroutine get_orb_sym
!================================================================================


!================================================================================
! converting a band-spin pattern into an index
  subroutine component2index(Nbands, ind, b1, s1, b2, s2, b3, s3, b4, s4)
    implicit none

    integer,intent(in) :: Nbands
    integer,intent(in) :: b1, s1, b2, s2, b3, s3, b4, s4
    integer,intent(inout) :: ind
    integer :: g1, g2, g3, g4

    g1=2*(b1-1) + s1
    g2=2*(b2-1) + s2
    g3=2*(b3-1) + s3
    g4=2*(b4-1) + s4

    ind =  8*Nbands**3*(g1-1) + 4*Nbands**2*(g2-1) + 2*Nbands*(g3-1) + g4

  end subroutine component2index
!=================================================================================

!================================================================================
! converting an index into a band-spin pattern
  subroutine index2component(Nbands, ind, b1, s1, b2, s2, b3, s3, b4, s4)

    implicit none
    integer,intent(in) :: Nbands,ind
    integer,intent(inout) :: b1, s1, b2, s2, b3, s3, b4, s4
    integer :: tmp1,tmp2,tmp3,ind_tmp
    integer :: g1,g2,g3,g4

    ! the proposed back conversion assumes the indices are
    ! given form 0 to max-1  
    ind_tmp = ind - 1
    tmp1 = 8*Nbands**3
    tmp2 = 4*Nbands**2
    tmp3 = 2*Nbands

    g1 = ind_tmp/tmp1 + 1
    g2 = (ind_tmp-tmp1*(g1-1))/tmp2 + 1
    g3 = (ind_tmp-tmp1*(g1-1)-tmp2*(g2-1))/tmp3 + 1
    g4 = (ind_tmp-tmp1*(g1-1)-tmp2*(g2-1)-tmp3*(g3-1)) + 1

    s1=mod(g1-1,2)+1
    b1=(g1-s1)/2+1

    s2=mod(g2-1,2)+1
    b2=(g2-s2)/2+1

    s3=mod(g3-1,2)+1
    b3=(g3-s3)/2+1

    s4=mod(g4-1,2)+1
    b4=(g4-s4)/2+1

  end subroutine index2component
!=================================================================================

!================================================================================
  subroutine component2index_band(Nbands, ind, b1, b2, b3, b4)
    implicit none
    integer,intent(in) :: Nbands
    integer,intent(in) :: b1, b2, b3, b4
    integer,intent(out) :: ind

    ind =  Nbands**3*(b1-1) + Nbands**2*(b2-1) + Nbands*(b3-1) + b4
  end subroutine component2index_band



  ! converting an index into a band pattern
  subroutine index2component_band(Nbands, ind, b1, b2, b3, b4)
    implicit none
    integer,intent(in) :: Nbands,ind
    integer,intent(out) :: b1, b2, b3, b4
    integer :: tmp1,tmp2,tmp3,ind_tmp

    ! the proposed back conversion assumes the indices are
    ! given form 0 to max-1
    ind_tmp = ind - 1
    tmp1 = Nbands**3
    tmp2 = Nbands**2
    tmp3 = Nbands

    b1 = ind_tmp/tmp1 + 1
    b2 = (ind_tmp-tmp1*(b1-1))/tmp2 + 1
    b3 = (ind_tmp-tmp1*(b1-1)-tmp2*(b2-1))/tmp3 + 1
    b4 = (ind_tmp-tmp1*(b1-1)-tmp2*(b2-1)-tmp3*(b3-1)) + 1
  end subroutine index2component_band


  ! return true if all 4 legs are in the same correlated subspace
  logical function index2cor(nineq,ndims,m,n,o,p)
    implicit none
    integer, intent(in) :: nineq
    integer, intent(in) :: ndims(nineq,2)
    ! band indices from specific beginning or end point of a 1PG
    integer, intent(in) :: m,n,o,p
    ! inequivalent atom number for specific index
    integer :: a,b,c,d
    integer :: dimstart,dimend,ineq, i

    a=0;b=0;c=0;d=0

    do ineq=1,nineq
      dimstart=1
      do i=2,ineq
        dimstart=dimstart+ndims(i-1,1)+ndims(i-1,2)
      enddo
      dimend=dimstart+ndims(ineq,1)-1 ! only in correlated sub space
      if ( m .ge. dimstart .and. m .le. dimend ) a=ineq
      if ( n .ge. dimstart .and. n .le. dimend ) b=ineq
      if ( o .ge. dimstart .and. o .le. dimend ) c=ineq
      if ( p .ge. dimstart .and. p .le. dimend ) d=ineq
    enddo

    ! checking if everything is on the same atom
    ! AND on correlated bands (non correlated lines would have ineq=0)
    if ( (a .eq. b) .and. (c .eq. d) .and. (a .eq. d) .and. (a .ne. 0)) then
      index2cor = .true.
    else
      index2cor = .false.
    endif
  end function index2cor

  ! return true if all 4 legs are in the same uncorrelated subspace
  logical function index2uncor(nineq,ndims,m,n,o,p)
    implicit none
    integer, intent(in) :: nineq
    integer, intent(in) :: ndims(nineq,2)
    ! band indices from specific beginning or end point of a 1PG
    integer, intent(in) :: m,n,o,p
    ! inequivalent atom number for specific index
    integer :: a,b,c,d
    integer :: dimstart,dimend,ineq, i

    a=0;b=0;c=0;d=0

    do ineq=1,nineq
      dimstart=ndims(1,1)+1
      do i=2,ineq
        dimstart=dimstart+ndims(i,1)+ndims(i-1,2)
      enddo
      dimend=dimstart+ndims(ineq,2)-1
      if ( m .ge. dimstart .and. m .le. dimend ) a=ineq
      if ( n .ge. dimstart .and. n .le. dimend ) b=ineq
      if ( o .ge. dimstart .and. o .le. dimend ) c=ineq
      if ( p .ge. dimstart .and. p .le. dimend ) d=ineq
    enddo

    if ( (a .eq. b) .and. (c .eq. d) .and. (a .eq. d) .and. (a .ne. 0)) then
      index2uncor = .true.
    else
      index2uncor = .false.
    endif
  end function index2uncor

  ! returns number of impurity if all 4 legs on the same
  ! otherwise returns 0
  integer function index2ineq(nineq,ndims,m,n,o,p)
    implicit none
    integer, intent(in) :: nineq
    integer, intent(in) :: ndims(nineq,2)
    ! band indices from specific beginning or end point of a 1PG
    integer, intent(in) :: m,n,o,p
    ! inequivalent atom number for specific index
    integer :: a,b,c,d
    integer :: dimstart,dimend,ineq, i

    a=0;b=0;c=0;d=0

    do ineq=1,nineq
      dimstart=1
      do i=2,ineq
        dimstart=dimstart+ndims(i-1,1)+ndims(i-1,2)
      enddo
      dimend=dimstart+ndims(ineq,1)+ndims(ineq,2)-1
      if ( m .ge. dimstart .and. m .le. dimend ) a=ineq
      if ( n .ge. dimstart .and. n .le. dimend ) b=ineq
      if ( o .ge. dimstart .and. o .le. dimend ) c=ineq
      if ( p .ge. dimstart .and. p .le. dimend ) d=ineq
    enddo

    if ( (a .eq. b) .and. (c .eq. d) .and. (a .eq. d) ) then
      index2ineq = a
    else
      index2ineq = 0
    endif
  end function index2ineq

  subroutine inverse_matrix_3(A)
    implicit none
    complex(kind=8) :: A(3,3),B(3,3),detinv

    ! Calculate the inverse determinant of the matrix
    detinv = 1.d0/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
                 - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
                 + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

    ! Calculate the inverse of the matrix
    B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
    B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
    B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
    B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
    B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
    B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
    B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
    B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
    B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))

    A=B
  end subroutine inverse_matrix_3

  subroutine inverse_matrix_2(A)
    implicit none
    complex(kind=8) :: A(2,2),B(2,2),detinv

   ! Calculate the inverse determinant of the matrix
    detinv = 1.d0/(A(1,1)*A(2,2) - A(1,2)*A(2,1))

    ! Calculate the inverse of the matrix
    B(1,1) = +detinv * A(2,2)
    B(2,1) = -detinv * A(2,1)
    B(1,2) = -detinv * A(1,2)
    B(2,2) = +detinv * A(1,1)

    A=B
  end subroutine inverse_matrix_2

end module
