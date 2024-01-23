program Horners
    implicit none
    integer,parameter :: n=4
    real,external :: f
    real :: x,x0,y,z
    real :: a(0:n)
    integer :: i

    i=1
    x0=1

    a(4)=2
    a(3)=0
    a(2)=-3
    a(1)=3
    a(0)=-4

    
    y=a(n)
    z=a(n)

    do i=n-1,1,-1
        y=x0*y+a(i)
        z=x0*z+y
        
    end do

    y=x0*y+a(0)


    write(*,*)y
    write(*,*)z

    

end program Horners

real function f(x)
real::x
f=2*x*x*x*x-3*x*x+3*x-4
!f=4*x*x*x+5*x*x+6*x+7
return
end function