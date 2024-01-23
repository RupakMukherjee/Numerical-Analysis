program linearfinitediff
  implicit none
  real,external::p,q,r

  integer,parameter::n=9
  integer::i
   real::a0,b0,h,alpha,beta,x
  real,dimension(0:n+1)::y,w
  real,dimension(n)::a,b,c,d,z,u,l

  a0=1.0
  b0=2.0
  alpha=1.0
  beta=2.0
  y(int(a0))=alpha
  y(int(b0))=beta


 h=(b0-a0)/(n+1)
 x=a0+h
 a(1)=2.0+(h**2)*q(x)
 b(1)=-1.0+(h/2)*p(x)
 d(1)=-(h**2)*r(x)+(1+(h/2.0)*p(x))*alpha



 do i=2,n-1,1
  x=a0+i*h
  a(i)=2.0+(h**2)*q(x)
 b(i)=-1+(h/2.0)*p(x)
 c(i)=-1.0-(h/2.0)*p(x)
 d(i)=-(h**2)*r(x)
 end do

 x=b0-h
 a(n)=2.0+(h**2)*q(x)
 c(n)=-1.0-(h/2)*p(x)
 d(n)=-(h**2)*r(x) +(1-(h/2)*p(x))*beta


 l(1)=a(1)
 u(1)=b(1)/a(1)
 z(1)=d(1)/l(1)

 do i=2,n-1,1
  l(i)=a(i)-c(i)*u(i-1)
  u(i)=b(i)/l(i)
  z(i)=(d(i)-c(i)*z(i-1))/l(i)
   end do

 l(n)=a(n)-c(n)*u(n-1)
 z(n)=(d(n)-c(n)*z(n-1))/l(n)
 w(0)=alpha
 w(n+1)=beta
 w(n)=z(n)

 do i= n-1,1,-1
  w(i)=z(i)-u(i)*w(i+1)
end do


do i=0,n+1,1
  x=a0+i*h

print*,'x(i):',x,'w(i):',w(i)
end do


stop


end program
real function p(g)
implicit none
real::g
p=-2.0/g
return
end function


real function q(g)
implicit none
real::g
q=2.0/(g**2)
return
end function

real function r(g)
implicit none
real::g
r=(sin(log(g)))/(g**2)
return
end function

