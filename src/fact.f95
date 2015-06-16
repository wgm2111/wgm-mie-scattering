
module fact
  
  ! public :: f
  
contains 
  
  recursive function f(nreal) result (factorial_result)
    
    real(kind=8), intent(in) :: nreal
    real(kind=8) :: factorial_result
    ! real, intent(in) :: nreal
    ! real :: factorial_result
    
    if (nreal<=0) then
       factorial_result = 1
    else 
       factorial_result = nreal * f(nreal-1)
    endif
  end function f
  
endmodule fact



program test_fact
  ! test fact
  use fact
  
  real(kind=8) :: a,b,c,d,e,g
  
  a = 1
  b = 2
  c = 3
  d = 4
  e = 5
  g = 10
 
  write(*,*) "factorial(",a,") = ", f(a)
  write(*,*) "factorial(",b,") = ", f(b)
  write(*,*) "factorial(",c,") = ", f(c)
  write(*,*) "factorial(",d,") = ", f(d)
  write(*,*) "factorial(",e,") = ", f(e)
  write(*,*) "factorial(",g,") = ", f(g)

end program test_fact
