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

module lookup_module
  use parameters_module
  implicit none
  private
  integer :: i, pst
  character(len=150) :: str_temp, str_split
  public :: string_find, int_find, int3_find, float_find, bool_find, group_find, subgroup_find

  contains

  subroutine string_find(search_string, save_string, search_start, search_end)
    character(*), intent(in)  :: search_string
    character(len=150), intent(inout) :: save_string ! keep default string
    integer, intent(in) :: search_start, search_end

    do i=search_start,search_end
      if (index(trim(file_save(i)),trim(search_string)) .ne. 0) then
        str_temp=file_save(i)
        pst=scan(str_temp,seperator)
        save_string=trim(adjustl(str_temp(pst+1:)))
      endif
    enddo
  end subroutine string_find

  subroutine int_find(search_string, save_int, search_start, search_end)
    character(*), intent(in)  :: search_string
    integer, intent(inout) :: save_int ! keep default values
    integer, intent(in) :: search_start, search_end

    do i=search_start,search_end
      if (index(trim(file_save(i)),trim(search_string)) .ne. 0) then
        str_temp=file_save(i)
        pst=scan(str_temp,seperator)
        str_temp=trim(adjustl(str_temp(pst+1:)))
        read(str_temp,*) save_int
      endif
    enddo
  end subroutine int_find

  subroutine int3_find(search_string, save_int1, save_int2, save_int3, search_start, search_end)
    character(*), intent(in)  :: search_string
    integer, intent(inout) :: save_int1, save_int2, save_int3 ! keep default values
    integer, intent(in) :: search_start, search_end

    do i=search_start,search_end
      if (index(trim(file_save(i)),trim(search_string)) .ne. 0) then
        str_temp=file_save(i)
        pst=scan(str_temp,seperator)
        str_temp=trim(adjustl(str_temp(pst+1:)))
        pst=scan(str_temp,multseperator)
        str_split=trim(adjustl(str_temp(:pst-1)))
        read(str_split,*) save_int1
        str_temp=trim(adjustl(str_temp(pst+1:)))
        pst=scan(str_temp,multseperator)
        str_split=trim(adjustl(str_temp(:pst-1)))
        str_temp=trim(adjustl(str_temp(pst+1:)))
        read(str_split,*) save_int2
        read(str_temp,*) save_int3
      endif
    enddo
  end subroutine int3_find

  subroutine float_find(search_string, save_float, search_start, search_end)
    character(*), intent(in)  :: search_string
    real(8), intent(inout) :: save_float ! keep default values
    integer, intent(in) :: search_start, search_end

    do i=search_start,search_end
      if (index(trim(file_save(i)),trim(search_string)) .ne. 0) then
        str_temp=file_save(i)
        pst=scan(str_temp,seperator)
        str_temp=trim(adjustl(str_temp(pst+1:)))
        read(str_temp,*) save_float
      endif
    enddo
  end subroutine float_find

  subroutine bool_find(search_string, save_bool, search_start, search_end)
    character(*), intent(in)  :: search_string
    logical, intent(inout) :: save_bool
    integer, intent(in) :: search_start, search_end

    do i=search_start,search_end
      if (index(trim(file_save(i)),trim(search_string)) .ne. 0) then
        str_temp=file_save(i)
        pst=scan(str_temp,seperator)
        str_temp=trim(adjustl(str_temp(pst+1:)))
        read(str_temp,*) save_bool
      endif
    enddo
  end subroutine bool_find

  subroutine group_find(search_string, save_start, save_end)
    character(*), intent(in) :: search_string
    integer, intent(out) :: save_start, save_end
    save_start=0
    save_end=0

    do i=1,lines
      if (index(trim(file_save(i)),trim(search_string)) .ne. 0) then
        save_start=i+1
        exit
      endif
    enddo

    if (save_start .ge. 1) then
      do i=save_start, lines
        if (index(trim(file_save(i)),'[') .eq. 1) then
          if (index(trim(file_save(i)),'[[') .eq. 1) then ! skip subgroups
            cycle
          endif
          save_end=i-1 ! one above the next session
          exit
        endif
      enddo

      if(save_end .eq. 0) then
        save_end = lines ! if nothing else is found, until the end of the file
      endif
    endif
    return
  end subroutine group_find

  subroutine subgroup_find(search_string, search_start, search_end, save_start, save_end)
    character(*), intent(in) :: search_string
    integer, intent(in) :: search_start, search_end
    integer, intent(out) :: save_start, save_end
    save_start=0
    save_end=0

    do i=search_start, search_end
      if (index(trim(file_save(i)),trim(search_string)) .ne. 0) then
        save_start=i+1
        exit
      endif
    enddo

    do i=save_start, search_end
      if (index(trim(file_save(i)),'[') .eq. 1) then
        save_end=i-1 ! one above the next session
        exit
      endif
    enddo

    if(save_end .eq. 0) then
      save_end = lines ! if nothing else is found, until the end of the file
    endif
  end subroutine subgroup_find

end module lookup_module
