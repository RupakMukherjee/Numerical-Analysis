program simpsondoubleintegral
    implicit none

    integer :: i,n,m,u
    real :: j1,j2,j3,h,a,b,x,y,L,k1,k2,k3,HX,j,Q
    real,external :: f,d,c

    !endpoints : 

    a=1.4
    b=2.0

    !even positive integers :

    m=2
    n=4

    !...................

    h=(b-a)/n
    j1=0
    j2=0
    j3=0

    do i=0,n,1
        x=a+i*h
        HX=(d(x)-c(x))/m
        k1=f(x,c(x))+f(x,d(x))
        k2=0
        k3=0
        do u=1,m-1,1
            y=c(x)+u*HX
            Q=f(x,y)

            if (MOD(u,2)==0) then
                k2=k2+Q

            else 
                k3=k3+Q
                
            end if
            
        end do

        L=(k1+2*k2+4*k3)*HX/3

        if (i==0.or.i==n) then
            j1=j1+L

        else if (MOD(i,2)==0) then
            j2=j2+L
        
        else 
            j3=j3+L
        end if
    

    end do

J=h*(j1+2*j2+4*j3)/3

print*,J

    
end program simpsondoubleintegral

real function f(x,y)
real,intent(in)::x,y
f=log(x+2*y)
end function

real function d(x)
real,intent(in)::x
d=1.5
end function

real function c(x)
real,intent(in)::x
c=1.0
end function