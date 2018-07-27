module settings

    use iso_fortran_env, only: r64 => real64
    use globals

    implicit none

    character(2**8) :: outfn = "disk"
    integer :: ngrid = 2**10

    integer, parameter ::   GRID_LINEAR = 1, &
                        &   GRID_LOG    = 2, &
                        &   GRID_ASINH  = 3, &
                        &   GRID_POW2   = 4
    integer :: tgrid = GRID_LOG
    real(r64) :: htop = 100
    logical :: cfg_auto_htop = .true.

contains

!----------------------------- READ_COMMAND_LINE ------------------------------!
!       reads command line arguments and stores them in global variables       !
!------------------------------------------------------------------------------!

  subroutine rdargvgl

    integer :: errno,i
    character(128) :: arg, nextarg

    iterate_arguments : do i=1,command_argument_count()

      call get_command_argument(i,arg)

      select case (arg)
      case("-opacity-ff", "-free-free", "-ff")
        use_opacity_ff = .TRUE.
      case("-no-opacity-ff", "-no-free-free", "-no-ff")
        use_opacity_ff = .FALSE.

      case("-opacity-bf", "-bound-free", "-bf")
        use_opacity_bf = .TRUE.
      case("-no-opacity-bf", "-no-bound-free", "-no-bf")
        use_opacity_bf = .FALSE.

      case("-opacity-planck", "-planck")
        use_opacity_planck = .TRUE.
      case("-no-opacity-planck", "-no-planck")
        use_opacity_planck = .FALSE.

      case("-conduction", "-cond")
        use_conduction = .TRUE.
      case("-no-conduction", "-no-cond")
        use_conduction = .FALSE.

      case ("-lin","-linear","-grid-linear")
        tgrid = GRID_LINEAR
      case ("-log","-grid-log")
        tgrid = GRID_LOG
      case ("-asinh","-grid-asinh")
        tgrid = GRID_ASINH
      case ("-grid-pow2", "-pow2")
        tgrid = GRID_POW2

      case ("-top","-htop","-z-top")
        call get_command_argument(i+1,nextarg)
        read (nextarg,*,iostat = errno) htop
        if ( errno .ne. 0 ) then
          error stop "top must be followed by an argument " &
          & // "(the upper boundary of computation interval)"
        end if
        ! if the height has been given, do not estimate
        cfg_auto_htop = .false.

      case ("-N","-n","-ngrid")
        call get_command_argument(i+1,nextarg)
        read (nextarg,*,iostat = errno) ngrid
        if ( errno .ne. 0 ) then
          error stop "ngrid must be followed by an argument " &
          & // "(number of bins)"
        end if

      case ("-output","-o")
        call get_command_argument(i+1,outfn)
        if ( len_trim(outfn) .eq. 0 ) then
          error stop "output must be followed by an argument " &
          & // "(output filename, without extension)"
        end if

      case ("-equilibrium","-dyfu")
        cfg_temperature_method = EQUATION_DIFFUSION
      case ("-balance","-corona")
        cfg_temperature_method = EQUATION_BALANCE
      case ("-compton")
        cfg_temperature_method = EQUATION_COMPTON
      end select

    end do iterate_arguments

  end subroutine

!--------------------------------------------------------------------------!

  subroutine rdconfgl(cfg)
    use confort
    type(config) :: cfg
    integer :: errno
    character(128) :: buf

    call mincf_get(cfg, 'mbh', buf, errno)
    if ( iand(errno, mincf_not_found) .ne. 0 ) then
      error stop "Key mbh (black hole mass) is REQUIRED!"
    end if
    read (buf,*) mbh

    call mincf_get(cfg, 'mdot', buf, errno)
    if ( iand(errno, mincf_not_found) .ne. 0 ) then
      error stop "Key mdot (accretion rate) is REQUIRED!"
    end if
    read (buf,*) mdot

    if ( mincf_exists(cfg,'output') ) then
      call mincf_get(cfg, 'output', outfn)
    end if
  end subroutine

!--------------------------------------------------------------------------!

  subroutine wpar_gl(u)
    use fileunits, only: fmpare, fmparec, fmparf, fmparl
    integer, intent(in) :: u
    write (u, fmparf) "X", abuX
    write (u, fmparf) "Z", abuZ
    write (u, fmpare) "kappa_abs_0", kapabs0(abuX,abuZ)
    write (u, fmpare) "kabs0", kapabs0(abuX,abuZ)
    write (u, fmpare) "kabp0", kapabp0(abuX,abuZ)
    write (u, fmpare) "htop", htop
    write (u, fmpare) "kappa_es", kapes0(abuX,abuZ)
    write (u, fmparec) "mbh", mbh, "black hole mass (x Msun)"
    write (u, fmparec) "mdot", mdot, "accretion rate (x Ledd)"
    write (u, fmparec) "rschw", rschw, 'Schwarzschild radius'
    write (u, fmparl) "use_opacity_ff", use_opacity_ff
    write (u, fmparl) "use_opacity_bf", use_opacity_bf
    write (u, fmparl) "use_planck", use_opacity_planck
    write (u, fmparl) "use_opacity_bb", .FALSE.
    write (u, fmparl) "conduction", use_conduction
  end subroutine

end module
