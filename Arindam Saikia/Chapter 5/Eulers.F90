program eulers
    implicit none
    
    integer :: i,N
    real :: h,a,b,t,yp,y,alpha,w
    real,external :: f

    print*,"Provide the end points a and b:"
    read(*,*)a,b

    print*,"Provide the value of N:"
    read(*,*)N

    print*,"Provide the initial condition, alpha:"
    read(*,*)alpha

    h=(b-a)/N
    t=a
    w=alpha

    print*,"The Output is:"

    do i=1,N,1
        w=w+h*f(t,w)
        t=a+i*h

    
    write(*,*)t,w
    end do

end program eulers

real function f(t,y)
real,intent(in) :: t,y
f=y-t*t+1.0
end function