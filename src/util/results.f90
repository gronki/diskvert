module results

    use iso_fortran_env, only: r64 => real64
    implicit none

    character (len=*) , parameter :: fmtdescr = '(A18,I4,A12,I10,I6)'
    character (len=*) , parameter :: fmtdescr2 = '(2A16)'

    character(len=*), parameter :: fmt_meta_ic  = '(A18,"  ",I11,"  # ",A)'
    character(len=*), parameter :: fmt_meta_i   = '(A18,"  ",I11)'
    character(len=*), parameter :: fmt_meta_fc  = '(A18,"  ",F11.4,"  # ",A)'
    character(len=*), parameter :: fmt_meta_f   = '(A18,"  ",F11.4)'
    character(len=*), parameter :: fmt_meta_ec  = '(A18,"  ",Es11.4,"  # ",A)'
    character(len=*), parameter :: fmt_meta_e   = '(A18,"  ",Es11.4)'
    character(len=*), parameter :: fmt_meta_sc  = '(A18,"  """,A,"""  # ",A)'
    character(len=*), parameter :: fmt_meta_s   = '(A18,"  """,A,"""")'

contains

    subroutine write_results(x,y,labels,outfn)

        real(r64), intent(in), dimension(:,:) :: y
        real(r64), intent(in), dimension(size(y,2)) :: x
        character(len=*), dimension(size(y,1)), intent(in) :: labels
        character(len=*), intent(in) :: outfn
        character(len=50) :: buf
        integer :: i

        integer :: nx,ny

        nx = size(y,2)
        ny = size(y,1)

        open(unit=34, file=(trim(outfn) // '.dat'), action='write')

        write (34,'(300A14)', advance="no") 'N', 'Z', labels

        do i=1,(2 + ny)
            write (buf, '("(",I0,")")') i
            write (34, '(A14)', advance="no") trim(buf)
        end do
        write (34,'(A)') ''

        do i = 1, nx
            write (34,'(I14,120ES14.5E3)')  i, x(i), y(:,i)
        end do

        close(34)

        open(unit=35, file=(trim(outfn) // '.col'), action='write')

        write(35,fmtdescr2) 'N', 'i4'
        write(35,fmtdescr2) 'Z', 'f4'

        do i=1, ny
            write (35,fmtdescr2) TRIM(labels(i)), 'f4'
        enddo

        close(35)


    end subroutine

    subroutine write_results_3(x,nx,y,dy,y_labels,ny,a,a_labels,na,outfn)

        integer, intent(in) :: nx,ny,na
        real(r64), intent(in), dimension(ny,nx) :: y, dy
        real(r64), intent(in), dimension(na,nx) :: a
        real(r64), intent(in), dimension(nx) :: x
        character(len=*), dimension(ny), intent(in) :: y_labels
        character(len=*), dimension(na), intent(in) :: a_labels
        character(len=*), intent(in) :: outfn
        integer :: i
        character(len=50) :: buf

        open(unit=34, file=(trim(outfn) // '.dat'), action='write')

        write (34,'(2A14)', advance="no") 'N', 'Z'

        do i=1,ny
            write (34, '(A14)', advance="no") trim(y_labels(i))
        end do
        do i=1,ny
            write (34, '(A14)', advance="no") "D*" // trim(y_labels(i))
        end do
        do i=1,na
            write (34, '(A14)', advance="no") trim(a_labels(i))
        end do
        write (34,'(A)') ''
        do i=1,(2 + 2*ny + na)
            write (buf, '("(",I0,")")') i
            write (34, '(A14)', advance="no") trim(buf)
        end do
        write (34,'(A)') ''

        do i = 1, nx
            write (34,'(I14,120ES14.5E3)')  i, x(i), y(:,i), dy(:,i), a(:,i)
        end do

        close(34)

        open(unit=35, file=(trim(outfn) // '.col'), action='write')

        write(35,fmtdescr2) 'N', 'i4'
        write(35,fmtdescr2) 'Z', 'f4'

        do i=1, ny
            write (35,fmtdescr2) TRIM(y_labels(i)), 'f4'
        enddo

        do i=1, ny
            write (35,fmtdescr2) "D*" // TRIM(y_labels(i)), 'f4'
        enddo

        do i=1, na
            write (35,fmtdescr2) TRIM(a_labels(i)), 'f4'
        enddo

        close(35)


    end subroutine

end module
