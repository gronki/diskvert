module fileunits

    use iso_fortran_env,        &
    only: uin => input_unit,    &
          uout => output_unit,  &
          uerr => error_unit
    implicit none

    character(*), parameter :: fmpari   = '(A20,1X,I12)'
    character(*), parameter :: fmparic  = '(A20,1X,I12,1X,"#",1X,A)'

    character(*), parameter :: fmparf   = '(A20,1X,F12.6)'
    character(*), parameter :: fmparfc  = '(A20,1X,F12.6,1X,"#",1X,A)'

    character(*), parameter :: fmpare   = '(A20,1X,Es12.4E3)'
    character(*), parameter :: fmparec  = '(A20,1X,Es12.4E3,1X,"#",1X,A)'

    character(*), parameter :: fmpars   = '(A20,1X,12A)'
    character(*), parameter :: fmparsc  = '(A20,1X,12A,1X,"#",1X,A)'

    character(*), parameter :: fmparl   = '(A20,1X,L12)'
    character(*), parameter :: fmparlc  = '(A20,1X,L12,1X,"#",1X,A)'

    character(*), parameter :: fmparg = '(a20, 1x, g12.5)'
    character(*), parameter :: fmpargc = '(a20, 1x, g12.5, 1x, "#", 1x, a)'

    character(*), parameter :: fmhdr  = '("#",1X,17("-"),1X,A)'

    character(*) , parameter :: fmcol = '(2A16)'


end module fileunits
