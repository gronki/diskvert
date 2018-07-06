module rk4settings

  use iso_fortran_env, only: r64 => real64
  use globals

  implicit none

  !   czy w rownaniach bilansu ma byc czlon comptonowski
  logical :: cfg_compton_term = .true.

  !   czy stosowac schemat eulera zamiast rk4
  logical :: cfg_euler_integration = .false.

  !   czy przeliczyc raz z danymi wartosciami centralnymi czy iterowac
  !   dla spelneia warunkow brzegowych?
  logical :: cfg_single_run = .false.

  !   jaki jest dopuszczalny blad iteracji
  real(r64) :: max_iteration_error = 1e-7

  !   czy umozliwic wylaczenie MRI?
  logical :: cfg_allow_mri_shutdown = .false.

  character, parameter :: EQUATION_MULTIBIL = 'M', &
      EQUATION_COMPTON_2 = 'X'

contains

  subroutine rdargvrk4
    integer :: i,errno
    character(2**8) :: arg, nextarg

    iterate_cmdline_arguments: do i = 1, command_argument_count()
      call get_command_argument(i, arg)
      select case (arg)
      case ("-dynamo-shutdown","-shutdown","-quench")
        cfg_allow_mri_shutdown = .TRUE.
      case ("-no-dynamo-shutdown","-no-shutdown","-no-quench")
        cfg_allow_mri_shutdown = .FALSE.
      case ("-single")
        cfg_single_run = .TRUE.
      case ("-euler","-no-rk4")
        cfg_euler_integration = .TRUE.
      case ("-no-euler","-rk4")
        cfg_euler_integration = .FALSE.
      case ("-compton")
        cfg_temperature_method = EQUATION_COMPTON
      case ("-compton2")
        cfg_temperature_method = EQUATION_COMPTON_2
      case ("-balance-multi","-multibil")
        cfg_temperature_method = EQUATION_MULTIBIL
      case ("-max-iteration-error","-precision")
        call get_command_argument(i+1,nextarg)
        read (nextarg,*,iostat=errno) max_iteration_error
        if ( errno .ne. 0 ) then
          error stop "-precision must be followed by an argument" &
          & // " (relative error to be reached until " &
          & // " iteration stops)"
        end if

      end select
    end do iterate_cmdline_arguments
  end subroutine

  subroutine tempmethod_c(m) bind(C, name = 'tempmethod')
    use iso_c_binding, only: c_char
    character(c_char), value :: m
    select case (m)
    case ('E')
      cfg_temperature_method = EQUATION_DIFFUSION
    case ('B')
      cfg_temperature_method = EQUATION_BALANCE
    case ('C')
      cfg_temperature_method = EQUATION_COMPTON
    case default
      error stop "method must be: E(quilibrium), B(alance) or C(ompton)"
    end select
  end subroutine

end module
