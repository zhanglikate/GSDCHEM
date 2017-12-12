module chem_io_mod

  use mpp_mod,         only : mpp_pe, mpp_root_pe, mpp_broadcast, mpp_scatter, mpp_gather
  use chem_domain_mod, only : chem_pelist, nvl_chem, nip, its, ite, jts, jte, map, indx
  use chem_config_mod, only : config => chem_config, glvl, yyyymmddhhmm
  use chem_config_mod, only : is_nuopc

  implicit none

  integer, parameter :: unitno = 80
  integer, parameter :: filename_len = 80
  logical, parameter :: check_decomp = .true.
  character(len=*), parameter :: subname = 'chem_data_read'
  character(len=*), parameter :: ArchvTimeUnit = 'ts'
  character(len=filename_len), dimension(:), allocatable :: list_outfiles


  private

  public :: chem_data_read, chem_data_read_global
  public :: chem_data_write

  interface chem_data_read
    module procedure chem_data_read_0d
    module procedure chem_data_read_1d
    module procedure chem_data_read_2d
    module procedure chem_data_mread_1d
    module procedure chem_data_mread_2d
    module procedure chem_data_mread_2d_1
  end interface chem_data_read

  interface chem_data_read_global
    module procedure chem_data_read_global_1d
  end interface chem_data_read_global

contains

  subroutine chem_data_read_global_1d(datafile, farray, rc, perm, label)
    character(len=*), intent(in) :: datafile
    real, dimension(:), intent(inout) :: farray
    integer, intent(out) :: rc
    integer, dimension(:), optional, intent(in) :: perm
    character(len=*), optional, intent(in) :: label

    ! -- local variables
    logical :: is_root_pe
    integer :: i

    ! -- begin
    rc = 0

    is_root_pe = (mpp_pe() == mpp_root_pe())

    if (is_root_pe) then
      open(unit=unitno, file=trim(datafile), form='unformatted', action='read', iostat=rc)
      if (rc == 0) rewind(unit=unitno, iostat=rc) 
    end if
    call mpp_broadcast(rc, mpp_root_pe())
    if (rc /= 0) then
      if (is_root_pe) close(unit=unitno)
      return
    end if

    farray = 0.
    if (is_root_pe) then
      read(unit=unitno, iostat=rc) farray
      close(unit=unitno)
    end if
    call mpp_broadcast(rc, mpp_root_pe())
    if (rc /= 0) return

    if (present(perm)) then
      if (is_root_pe) call do_perm(farray, perm, rc)
      call mpp_broadcast(rc, mpp_root_pe())
      if (rc /= 0) return
      if (is_nuopc) then
        if (is_root_pe) call do_perm(farray, indx, rc)
        call mpp_broadcast(rc, mpp_root_pe())
        if (rc /= 0) return
      end if
    end if

    call mpp_broadcast(farray, size(farray), mpp_root_pe())
    if (check_decomp .and. is_root_pe) then
      if (present(label)) then
        write(6,'(a,":",a," - reading  - min/max = ",2g16.6)') subname, trim(label), minval(farray), maxval(farray)
      else
        write(6,'(a,":",a," - reading  - min/max = ",2g16.6)') subname, '', minval(farray), maxval(farray)
      end if
    end if

  end subroutine chem_data_read_global_1d

  subroutine chem_data_read_0d(datafile, ival, rc, label)
    character(len=*), intent(in) :: datafile
    integer, intent(out) :: ival
    integer, intent(out) :: rc
    character(len=*), optional, intent(in) :: label

    ! -- local variables
    logical :: is_root_pe

    ! -- begin
    rc = 0

    is_root_pe = (mpp_pe() == mpp_root_pe())

    if (is_root_pe) then
      open(unit=unitno, file=trim(datafile), form='unformatted', action='read', iostat=rc)
      if (rc == 0) rewind(unit=unitno, iostat=rc) 
    end if
    call mpp_broadcast(rc, mpp_root_pe())
    if (rc /= 0) then
      if (is_root_pe) close(unit=unitno)
      return
    end if

    if (is_root_pe) then
      read(unit=unitno, iostat=rc) ival
      close(unit=unitno)
    end if
    call mpp_broadcast(rc, mpp_root_pe())
    if (rc /= 0) return

    call mpp_broadcast(ival, mpp_root_pe())

    if (check_decomp .and. is_root_pe) then
      if (present(label)) then
        write(6,'(a,":",a," - reading  - value = ",i0)') subname, trim(label), ival
      else
        write(6,'(a,":",a," - reading  - value = ",i0)') subname, '', ival
      end if
    end if
    
  end subroutine chem_data_read_0d

  subroutine chem_data_read_1d(datafile, farray, rc, perm, label)
    character(len=*), intent(in) :: datafile
    real, dimension(jts:jte), intent(inout) :: farray
    integer, intent(out) :: rc
    integer, dimension(:), optional, intent(in) :: perm
    character(len=*), optional, intent(in) :: label

    ! -- local variables
    logical :: is_root_pe
    real, dimension(nip, 1)  :: glob
    real, dimension(jts:jte, 1 ) :: loc

    ! -- begin
    rc = 0

    is_root_pe = (mpp_pe() == mpp_root_pe())

    if (is_root_pe) then
      open(unit=unitno, file=trim(datafile), form='unformatted', action='read', iostat=rc)
      if (rc == 0) rewind(unit=unitno, iostat=rc) 
    end if
    call mpp_broadcast(rc, mpp_root_pe())
    if (rc /= 0) then
      if (is_root_pe) close(unit=unitno)
      return
    end if

    glob = 0.
    if (is_root_pe) then
      read(unit=unitno, iostat=rc) glob(:,1)
      close(unit=unitno)
    end if
    call mpp_broadcast(rc, mpp_root_pe())
    if (rc /= 0) return

    if (present(perm)) then
      if (is_root_pe) call do_perm(glob(:,1), perm, rc)
      call mpp_broadcast(rc, mpp_root_pe())
      if (rc /= 0) return
      if (is_nuopc) then
        if (is_root_pe) call do_perm(glob(:,1), indx, rc)
        call mpp_broadcast(rc, mpp_root_pe())
        if (rc /= 0) return
      end if
    end if

    if (check_decomp .and. is_root_pe) then
      if (present(label)) then
        write(6,'(a,":",a," - reading  - min/max = ",2g16.6)') subname, trim(label), minval(glob), maxval(glob)
      else
        write(6,'(a,":",a," - reading  - min/max = ",2g16.6)') subname, '', minval(glob), maxval(glob)
      end if
    end if

    loc = 0.
    call mpp_scatter(jts, jte, 1, 1, chem_pelist, loc, glob, is_root_pe)
    farray = loc(:,1)

    if (check_decomp) then
      glob = 0.
      loc(:,1) = farray
      call mpp_gather(jts, jte, 1, 1, chem_pelist, loc, glob, is_root_pe)
      if (is_root_pe) then
        if (present(label)) then
          write(6,'(a,":",a," - checking - min/max = ",2g16.6)') subname, trim(label), minval(glob), maxval(glob)
        else
          write(6,'(a,":",a," - checking - min/max = ",2g16.6)') subname, '', minval(glob), maxval(glob)
        end if
      end if
    end if
    
  end subroutine chem_data_read_1d

  subroutine chem_data_mread_2d(datafile, farray, start_rec, num_rec, rc, perm, label)
    character(len=*), intent(in) :: datafile
    integer, intent(in) :: start_rec, num_rec
    real, dimension(jts:jte, num_rec), intent(inout) :: farray
    integer, intent(out) :: rc
    integer, dimension(:), optional, intent(in) :: perm
    character(len=*), optional, intent(in) :: label

    ! -- local variables
    logical :: is_root_pe
    integer :: irec
    real, dimension(nip, 1)  :: glob

    ! -- begin
    rc = 0

    is_root_pe = (mpp_pe() == mpp_root_pe())

    if (is_root_pe) then
      open(unit=unitno, file=trim(datafile), form='unformatted', action='read', iostat=rc)
      if (rc == 0) rewind(unit=unitno, iostat=rc) 
    end if
    call mpp_broadcast(rc, mpp_root_pe())
    if (rc /= 0) then
      if (is_root_pe) close(unit=unitno)
      return
    end if

    irec = 1
    do while (rc == 0 .and. irec < start_rec)
      if (is_root_pe) read(unit=unitno, iostat=rc)
      call mpp_broadcast(rc, mpp_root_pe())
      if (rc /= 0) exit
      irec = irec + 1
    end do

    if (rc /= 0) then
      if (is_root_pe) close(unit=unitno)
      return
    end if


    do irec = start_rec, num_rec
      glob = 0.
      if (is_root_pe) read(unit=unitno, iostat=rc) glob(:,1)
      call mpp_broadcast(rc, mpp_root_pe())
      if (rc /= 0) exit
      if (present(perm)) then
        if (is_root_pe) call do_perm(glob(:,1), perm, rc)
        call mpp_broadcast(rc, mpp_root_pe())
        if (rc /= 0) return
        if (is_nuopc) then
          if (is_root_pe) call do_perm(glob(:,1), indx, rc)
          call mpp_broadcast(rc, mpp_root_pe())
          if (rc /= 0) return
        end if
      end if

      if (check_decomp .and. is_root_pe) then
        if (present(label)) then
          write(6,'(a,":",a," - reading  - min/max = ",2g16.6)') subname, trim(label), minval(glob), maxval(glob)
        else
          write(6,'(a,":",a," - reading  - min/max = ",2g16.6)') subname, '', minval(glob), maxval(glob)
        end if
      end if

      call mpp_scatter(jts, jte, 1, 1, chem_pelist, farray(:,irec:irec), glob, is_root_pe)

      if (check_decomp) then
        glob = 0.
        call mpp_gather(jts, jte, 1, 1, chem_pelist, farray(:,irec:irec), glob, is_root_pe)
        if (is_root_pe) then
          if (present(label)) then
            write(6,'(a,":",a," - checking - min/max = ",2g16.6)') subname, trim(label), minval(glob), maxval(glob)
          else
            write(6,'(a,":",a," - checking - min/max = ",2g16.6)') subname, '', minval(glob), maxval(glob)
          end if
        end if
      end if
    end do

    if (is_root_pe) close(unit=unitno)
    
  end subroutine chem_data_mread_2d

  subroutine chem_data_mread_1d(datafile, farray, irec, rc, perm, label)
    character(len=*), intent(in) :: datafile
    integer, intent(in) :: irec
    real, dimension(jts:jte), intent(inout) :: farray
    integer, intent(out) :: rc
    integer, dimension(:), optional, intent(in) :: perm
    character(len=*), optional, intent(in) :: label

    ! -- local variables
    logical :: is_root_pe
    real, dimension(:,:), allocatable :: loc

    allocate(loc(jts:jte,1))

    call chem_data_mread_2d(datafile, loc, irec, 1, rc, perm=perm, label=label)

    farray = loc(:,1)

  end subroutine chem_data_mread_1d

  subroutine chem_data_read_2d(datafile, farray, rc, perm, label)
    character(len=*), intent(in) :: datafile
    real, dimension(:,:), intent(inout) :: farray
    integer, intent(out) :: rc
    integer, dimension(:), optional, intent(in) :: perm
    character(len=*), optional, intent(in) :: label

    integer :: nlev

    nlev = size(farray, dim=1)
    call chem_data_read_2d_domain(datafile, 1, nlev, jts, jte, chem_pelist, farray, rc, perm=perm, label=label)
    
  end subroutine chem_data_read_2d

  subroutine chem_data_mread_2d_1(datafile, farray, start_rec, rc, label)
    character(len=*), intent(in) :: datafile
    real, dimension(:,:), intent(inout) :: farray
    integer, intent(in) :: start_rec
    integer, intent(out) :: rc
    character(len=*), optional, intent(in) :: label

    integer :: nlev

    nlev = size(farray, dim=1)
    call chem_data_read_2d_domain(datafile, 1, nlev, jts, jte, chem_pelist, farray, rc, start_rec=start_rec, label=label)
    
  end subroutine chem_data_mread_2d_1

  subroutine chem_data_read_2d_domain(datafile, is, ie, js, je, pelist, farray, rc, start_rec, perm, label)
    character(len=*), intent(in) :: datafile
    integer, intent(in) :: is, ie, js, je
    integer, dimension(:), intent(in) :: pelist
    real, dimension(is:ie,js:je), intent(inout) :: farray
    integer, intent(out) :: rc
    integer, dimension(:), optional, intent(in) :: perm
    integer, optional, intent(in) :: start_rec
    character(len=*), optional, intent(in) :: label

    ! -- local variables
    logical :: is_root_pe
    integer :: i
    real, dimension(is:ie, nip)  :: glob

    ! -- begin
    rc = 0

    is_root_pe = (mpp_pe() == mpp_root_pe())

    if (is_root_pe) then
      open(unit=unitno, file=trim(datafile), form='unformatted', status='old', iostat=rc)
      if (rc == 0) rewind(unit=unitno, iostat=rc) 
    end if

    call mpp_broadcast(rc, mpp_root_pe())
    if (rc /= 0) then
      if (is_root_pe) close(unit=unitno)
      return
    end if

    if (present(start_rec)) then
      i = 1
      do while (rc == 0 .and. i < start_rec)
        if (is_root_pe) read(unit=unitno, iostat=rc)
        call mpp_broadcast(rc, mpp_root_pe())
        if (rc /= 0) exit
        i = i + 1
      end do
      if (rc /= 0) then
        if (is_root_pe) close(unit=unitno)
        return
      end if
    end if

    glob = 0.
    if (is_root_pe) then
      read(unit=unitno, iostat=rc) glob
      close(unit=unitno)
    end if
    call mpp_broadcast(rc, mpp_root_pe())
    if (rc /= 0) return

    if (present(perm)) then
      if (is_root_pe) then
        do i = is, ie
          call do_perm(glob(i,:), perm, rc)
          if (rc /= 0) exit
        end do
      end if
      call mpp_broadcast(rc, mpp_root_pe())
      if (rc /= 0) return
      if (is_nuopc) then
        if (is_root_pe) then
          do i = is, ie
            call do_perm(glob(i,:), indx, rc)
            if (rc /= 0) exit
          end do
        end if
        call mpp_broadcast(rc, mpp_root_pe())
        if (rc /= 0) return
      end if
    end if

    if (check_decomp .and. is_root_pe) then
      if (present(label)) then
        write(6,'(a,":",a," - reading  - min/max = ",2g16.6)') subname, trim(label), minval(glob), maxval(glob)
      else
        write(6,'(a,":",a," - reading  - min/max = ",2g16.6)') subname, '', minval(glob), maxval(glob)
      end if
    end if

    call mpp_scatter(is, ie, js, je, pelist, farray, glob, is_root_pe)

    if (check_decomp) then
      glob = 0.
      call mpp_gather(is, ie, js, je, pelist, farray, glob, is_root_pe)
      if (is_root_pe) then
        if (present(label)) then
          write(6,'(a,":",a," - checking - min/max = ",2g16.6)') subname, trim(label), minval(glob), maxval(glob)
        else
          write(6,'(a,":",a," - checking - min/max = ",2g16.6)') subname, '', minval(glob), maxval(glob)
        end if
      end if
    end if
    
  end subroutine chem_data_read_2d_domain

  subroutine do_perm(array, perm, rc)

    real,    dimension(:), intent(inout) :: array
    integer, dimension(:), intent(in   ) :: perm
    integer,               intent(out  ) :: rc

    integer :: i
    real, dimension(:), allocatable :: tmp

    rc = -1
    if (size(array) /= size(perm)) then
      print *,'do_perm: array size /= perm size'
      return
    end if

    allocate(tmp(size(array)))
    tmp = array
    do i = 1, size(perm)
      array(i) = tmp(perm(i))
    end do
    deallocate(tmp)
    rc = 0

  end subroutine do_perm

  ! --------------------------------------------------
  !  write methods
  ! --------------------------------------------------

  subroutine chem_data_write (its, name, field, nlev, scalefactor, accum_start, twodfile, filebase)
    integer,          intent(in) :: its           ! iteration
    character(len=*), intent(in) :: name          ! 4-character field name
    integer,          intent(in) :: nlev          ! number of vertical levels
    real,             intent(in) :: field(nlev,jts:jte)   ! data to be written !sms$distribute end
    real,             intent(in), optional :: scalefactor ! some GRIB fields need a scale factor applied
    integer,          intent(in), optional :: accum_start ! some GRIB fields need an accumulation count
    logical,          intent(in), optional :: twodfile    ! put on 2d file
    character(len=*), intent(in), optional :: filebase

    ! -- local variables
    logical :: is_root_pe, twodfile_local
    real, dimension(:,:), allocatable :: glob
    character(len=filename_len) :: fn, fn_base

    integer :: i, rc

    ! -- begin
    is_root_pe = (mpp_pe() == mpp_root_pe())
    twodfile_local = .false.

    if (is_root_pe) then
      if (present(filebase)) then
        fn_base = trim(filebase)
      else
        fn_base = trim(config % chem_hist_outname)
      end if
      if (present(twodfile)) twodfile_local = twodfile
      fn = ''
      if (twodfile_local) then
        write(fn,'(a,"2D__",i6.6,a)') trim(fn_base), its, ArchvTimeUnit
      else
        write(fn,'(2a,i6.6,a)') trim(fn_base), trim(name), its, ArchvTimeUnit
      end if
      allocate(glob(nlev, nip))
      glob = 0.
    end if

    call mpp_gather(1, nlev, jts, jte, chem_pelist, field, glob, is_root_pe)
    if (is_root_pe) then
      if (present(scalefactor)) glob = scalefactor * glob
    end if

    if (is_root_pe) then
      if (is_file_append(fn)) then
        open(unit=unitno, file=trim(fn), form='unformatted', status='old', position='append')
      else
        open(unit=unitno, file=trim(fn), form='unformatted', status='unknown', position='rewind')
      end if
      write(unit=unitno) header(name, nlev, its)
      if (is_nuopc) then
        do i = 1, nlev
          call do_perm(glob(i,:), map, rc)
        end do
      end if
      write(unit=unitno) glob
      close(unit=unitno)
      if (check_decomp) write(6,'(a,":",a," - written - min/max = ",2g16.6)') 'chem_data_write', name, minval(glob), maxval(glob)
      deallocate(glob)
    end if
    
  end subroutine chem_data_write

  logical function is_file_append(filename)

    character(len=*), intent(in) :: filename
    character(len=filename_len), dimension(:), allocatable :: tmp

    integer :: i, lsize

    is_file_append = .false.

    if (.not.allocated(list_outfiles)) then
      allocate(list_outfiles(1))
      list_outfiles(1) = trim(filename)
      return
    end if

    i = 1
    lsize = size(list_outfiles)
    do while (.not.is_file_append .and. i <= lsize)
      is_file_append = (trim(list_outfiles(i)) == trim(filename))
      i = i + 1
    end do

    if (.not.is_file_append) then
      allocate(tmp(lsize))
      tmp = list_outfiles
      deallocate(list_outfiles)
      allocate(list_outfiles(lsize+1))
      list_outfiles(1:lsize) = tmp
      list_outfiles(lsize+1) = trim(filename)
      deallocate(tmp)
    end if

  end function is_file_append

  function header(varname,levels,its)
    implicit none
    character*(*),intent(in)::varname
    integer,intent(in)::its,levels

    integer,parameter::header_cols=80
    integer,parameter::header_rows=10
    integer,parameter::header_len=header_cols*header_rows

    character(len=header_cols)::line
    character(len=header_len)::h

    character(len=header_len)::header
    integer :: pos

    pos=1
    h=' '
    write (line,'(a,a,a,a)') 'FIM ',varname,' Forecast initial time YYYYMMDDHHMM: ',yyyymmddhhmm
    call append
    write (line,'(a,i0,a,i0,a,i0,a,i0,a,a)') 'Level ',levels,', GLVL= ',glvl,', Step ',its,', ',its,' ',ArchvTimeUnit
    call append
    write (line,'(a,i0,a,i0)') 'dim1=',levels,', nip=',nip
    call append
    write (line,'(i0)') 4
    call append
    write (line,'(i0)') 5
    call append
    write (line,'(i0)') 6
    call append
    write (line,'(i0)') 7
    call append
    write (line,'(i0)') 8
    call append
    write (line,'(i0)') 9
    call append
    write (line,'(i0)') 10
    call append
    header=h

  contains

    subroutine append
    implicit none
    integer::i,j

    if (pos.ge.header_len) then
      write (*,'(a)') 'ERROR: Attempt to write past end of header.'
      stop
    endif
    do i=1,len(line)
      j=pos+i-1
      h(j:j)=line(i:i)
    enddo
    pos=(int((pos+header_cols)/header_cols)*header_cols)+1
    end subroutine append

  end function header

end module chem_io_mod
