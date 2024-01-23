program adaptivequadrature
implicit none

integer,parameter::n=100
integer :: i
real :: a0,b0,Tol,APP,FD,FE,Z,S1,S2
real, dimension(1:n) :: a, T, h, FA, FC, FB, S, L
real, dimension(1:8) :: v
real,external::f

a0=0
b0=1
Tol=1e-3
APP=0
i=1
T(i)=10*Tol
a(i)=a0
h(i)=(b0-a0)/2
FA(i)=f(a0)
FC(i)=f(a0+h(i))
FB(i)=f(b0)
S(i)=h(i)*(FA(i)+4*FC(i)+FB(i))/3
L(i)=1

do while(i>0)
    FD=f(a(i)+h(i)/2)
    FE=f(a(i)+3*h(i)/2)
    S1=h(i)*(FA(i)+4*FD+FC(i))/6
    S2=h(i)*(FC(i)+4*FE+FB(i))/6

    v(1)=a(i)
    v(2)=FA(i)
    v(3)=FC(i)
    v(4)=FB(i)
    v(5)=h(i)
    v(6)=T(i)
    v(7)=S(i)
    v(8)=L(i)

    i=i-1
    Z=S1+S2-v(7)
    if(abs(Z)<v(6)) then
        APP=APP+(S(1)+S(2))
    else if(v(8)>=n) then
        print*,"Level Exceeded"
        stop
    else
        i=i+1
        a(i)=v(1)+v(5)
        FA(i)=v(3)
        FC(i)=FE
        FB(i)=v(4)
        h(i)=v(5)/2.0
        T(i)=v(6)/2.0
        S(i)=S(2)
        L(i)=v(8)+1

        i=i+1
        a(i)=v(1)
        FA(i)=v(2)
        FC(i)=FD
        FB(i)=v(3)
        h(i)=h(i-1)
        T(i)=T(i-1)
        S(i)=S(1)
        L(i)=L(i-1)

    end if

end do
!if (APP<Tol) Then
    Print*,"output is",APP
!end if

end program adaptivequadrature

real function f(x)
        real, intent(in) :: x
        !f = (100/x**2)*sin(10/x)
        f=sin(x)
    end function f