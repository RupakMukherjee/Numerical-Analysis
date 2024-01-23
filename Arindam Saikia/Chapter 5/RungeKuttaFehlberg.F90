program rungekuttafehlberg
    implicit none

    integer :: i,N
    real :: a,b,t,y,alpha,TOL,hmin,hmax,FLAG,w,h,k1,k2,k3,k4,k5,k6,R,delta
    real,external :: f
    
    print*,"Provide the end points a and b:"
    read(*,*)a,b

    print*,"Provide the initial condition, alpha:"
    read(*,*)alpha

    print*,"Provide the value of TOL:"
    read(*,*)TOL

    print*,"Provide the value of hmax and hmin:"
    read(*,*)hmax,hmin

    t=a
    w=alpha
    h=hmax
    FLAG=1

    print*,"The Output is :"

    do while (FLAG==1)
        k1 = h * f(t, w)
        k2 = h * f(t + h / 4.0, w + k1 / 4.0)
        k3 = h * f(t + 3.0 * h / 8.0, w + 3.0 * k1 / 32.0 + 9.0 * k2 / 32.0)
        k4 = h * f(t + 12.0 * h / 13.0, w + 1932.0 * k1 / 2197.0 - 7200.0 * k2 / 2197.0 + 7296.0 * k3 / 2197.0)
        k5 = h * f(t + h / 2.0, w + 439.0 * k1 / 216.0 - 8.0 * k2 + 3680.0 * k3 / 513.0 - 845.0 * k4 / 4104.0)
        k6 = h * f(t + h, w - 8.0 * k1 / 27.0 + 2.0 * k2 - 3544.0 * k4 / 2565.0 + 1859.0 * k5 / 4104.0 - 11.0 * k6 / 40.0)
        
        R = 1.0 / h * abs(k1 / 360.0 - 128.0 * k3 / 4275.0 - 2197.0 * k4 / 75240.0 + k5 / 50.0 + 2.0 * k6 / 55.0)

    if (R <= TOL) then 
            t = t + h
            w = w + 25.0 * k1 / 216.0 + 1408.0 * k3 / 2565.0 + 2197.0 * k4 / 4104.0 - k5 / 5.0
            print *, t, w, h
        end if

    delta=0.84*(TOL/R)**0.25
    
    if (delta <= 0.1) then
        h=0.1*h
    
    else if (delta >= 4.0) then
        h=4.0*h
    else
        h=delta*h
    end if
   

    if (h > hmax) then
       h=hmax
    end if

    if (t >= b) then
        FLAG = 0
    else if (t + h > b) then
        h = b - t
    else if (h < hmin) then
        FLAG = 0
        print *, "Minimum h exceeded!"
    end if

    end do

    


end program rungekuttafehlberg

real function f(t,y)
real,intent(in) :: t,y
f=y-t*t+1.0
end function