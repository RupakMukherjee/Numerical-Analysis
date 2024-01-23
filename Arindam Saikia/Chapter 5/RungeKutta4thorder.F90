program rungekutta4thorder
    implicit none
    
    integer :: i,N
    real :: h,a,b,t,yp,y,alpha,w,k1,k2,k3,k4
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

    do i =1,N,1
        k1=h*f(t,w)
        k2=h*f(t+h/2.0,w+k1/2.0)
        k3=h*f(t+h/2.0,w+k2/2.0)
        k4=h*f(t+h,w+k3)

        w=w+(k1+2*k2+2*k3+k4)/6.0
        t=a+i*h

        write(*,*)t,w

    end do


end program rungekutta4thorder

real function f(t,y)
real,intent(in) :: t,y
f=y-t*t+1.0
end function