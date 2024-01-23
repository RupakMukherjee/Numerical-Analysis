program Muller
implicit none
real :: Tol
real :: p,p0,p1,p2,h1,h2,d0,D,q1,q2,b,E,h
real,external :: f
integer :: i=3
integer,parameter :: n=1000


p0=0.5
p1=1.0
p2=1.5
h1=p1-p0
h2=p2-p1
q1=(f(p1)-f(p0))/h1
q2=(f(p2)-f(p1))/h2
d0=(q2-q1)/(h2+h1)

Tol=1E-7

do while(i<=n)

   b=q2+h2*d0
   D=(b*b-4*f(p2)*d0)**(1.0/2.0)
    
if (abs(b-D)<abs(b+D)) then
    E=b+D
else  
    E=b-D
end if

h=(-2.0*f(p2))/E
p=p2+h

if (abs(h)<Tol) then
    write(*,*)p
    stop
end if

i=i+1
p0=p1
p1=p2
p2=p
h1=p1-p0
h2=p2-p1
q1=(f(p1)-f(p0))/h1
q2=(f(p2)-f(p1))/h2
d0=(q2-q1)/(h2+h1)

end do

print*,"Method failed after n iterations,n=",n

end program Muller

real function f(x)
real,intent(in)::x
f=x*x*x*x-3*x*x*x+x*x+x+1
return
end function
