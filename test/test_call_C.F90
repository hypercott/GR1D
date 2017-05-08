
program test_call_C

  implicit none

  external testic

  integer :: a,b,c

  a = 3
  b = 4
  c = 0

  call testic(a,b,c)

  write(6,*) c


end program test_call_C
