program gaussiandoubleintegral
    implicit none

    integer,parameter :: n=3,m=3
    integer :: i,u
    real:: h1,h2,J,JX,a,b,c1,d1,x,k1,k2,Q,y
    real,dimension(1:5,1:5) :: r,c0
    real,external::f,c,d

    !Endpoints:
    a=1.4
    b=2.0

    r(2,1)=0.5773502692
    r(2,2)=-0.5773502692
    r(3,1)=0.7745966692
    r(3,2)=0.0000000000
    r(3,3)=-0.7745966692
    r(4,1)=0.8611363116
    r(4,2)=0.3399810436
    r(4,3)=-0.3399810436
    r(4,4)=-0.8611363116
    r(5,1)=0.9061798459
    r(5,2)=0.5384693101
    r(5,3)=0.0000000000
    r(5,4)=-0.5384693101
    r(5,5)=-0.9061798459


    c0(2,1)=1.0000000000
    c0(2,2)=1.0000000000
    c0(3,1)=0.5555555556
    c0(3,2)=0.8888888886
    c0(3,3)=0.5555555556
    c0(4,1)=0.3478548451
    c0(4,2)=0.6521451549
    c0(4,3)=0.6521451549
    c0(4,4)=0.3478548451
    c0(5,1)=0.2369268850
    c0(5,2)=0.4786286705
    c0(5,3)=0.5688888889
    c0(5,4)=0.4786286705
    c0(5,5)=0.2369268850

h1=(b-a)/2
h2=(b+a)/2
J=0

do i=1,m,1
    JX=0
    x=h1*r(m,i)+h2
    d1=d(x)
    c1=c(x)
    k1=(d1-c1)/2
    k2=(d1+c1)/2

    do u=1,n,1
        y=k1*r(n,u)+k2
        Q=f(x,y)
        JX=JX+c0(n,u)*Q
    end do

    J=J+c0(m,i)*k1*JX
end do

J=h1*J

print*, "The Output is :",J

    
end program gaussiandoubleintegral

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