program trapezoidal
    implicit none

    integer :: i,j,N,M 
    real :: alpha,Tol,a,b,h,FLAG,k1,w0,t,w
    real,external :: f,fy

    !INPUT
    a=0
    b=1
    N=5
    alpha=-1
    Tol=10E-4
    M=10

    h=(b-a)/N 
    t=a
    w=alpha

    write(*,*)t,w

    do i=1,N 
        k1=w+h*f(t,w)/2.0
        w0=k1 
        j=1
        FLAG=0

        do while (FLAG==0)
            w=w0-(w0-(h*f(t+h,w0)/2.0)-k1)/(1-h*fy(t+h,w0))
        
        if (abs(w-w0)<Tol) then
                FLAG=1
            else 
                j=j+1
                w0=w
            if (j>M) then
                print*,"The Maximum number of iterations exceeded"
            stop
            end if
        end if
        
        end do
    
        t=a+i*h
        write(*,*)t,w
    end do

    
end program trapezoidal

real function f(t,w)
real,intent(in) :: t,w 
f=5*exp(5*t)*(w-t)*(w-t)+1
end function

real function fy(t,w)
real,intent(in) :: t,w 
fy=10*exp(5*t)*(w-t)
end function