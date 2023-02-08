module test
implicit none
integer, parameter :: a = 1
integer, parameter :: b = 1

end module test

program main
use test
implicit none

  call test1(a)

    
contains
 subroutine test1(a)
 integer a
 
 print*, b
 end subroutine test1

end program main
