program broyden
  implicit none
  integer,parameter::n=3
  integer::i,j,k,NN
  real::ja(n,n),ja0(n,n),ff(n,1),gg(n,1),ff0(n,1),x(n,1),x0(n,1),AA(n,n),v(n,1),u(n,1),BB(n,n),s(n,1),y(n,1),w(n,1),z(n,1)
  real::TOL,p(1,1),l(n,1),v0(n,1),AA0(n,n)
  real::ss(1,n),uu(1,n)
  real,external::f1,f2,f3

  NN=10
  TOL=0.1
  x0(1,1)=0.1
  x0(2,1)=0.1
  x0(3,1)=-0.1
  ja0(1,1)=3
  ja0(1,2)=x0(3,1)*sin(x0(2,1)*x0(3,1))
  ja0(1,3)=x0(2,1)*sin(x0(2,1)*x0(3,1))
  ja0(2,1)=2*x0(1,1)
  ja0(2,2)=-162*(x0(2,1)+0.1)
  ja0(2,3)=cos(x0(3,1))
  ja0(3,1)=-x0(2,1)*exp(-x0(1,1)*x0(2,1))
  ja0(3,2)=-x0(1,1)*exp(-x0(1,1)*x0(2,1))
  ja0(3,3)=20
  ff0(1, 1) = f1(x0(1, 1), x0(2, 1), x0(3, 1))
  ff0(2, 1) = f2(x0(1, 1), x0(2, 1), x0(3, 1))
  ff0(3, 1) = f3(x0(1, 1), x0(2, 1), x0(3, 1))
  ff(1, 1) = f1(x(1, 1), x(2, 1), x(3, 1))
  ff(2, 1) = f2(x(1, 1), x(2, 1), x(3, 1))
  ff(3, 1) = f3(x(1, 1), x(2, 1), x(3, 1))

 AA0=ja0
 v0=ff0

 call inverse(AA0,BB,n)

 AA=bb

 s=-matmul(AA,v0)

 x=x0+s
ff(1, 1) = f1(x(1, 1), x(2, 1), x(3, 1))
ff(2, 1) = f2(x(1, 1), x(2, 1), x(3, 1))
ff(3, 1) = f3(x(1, 1), x(2, 1), x(3, 1))

k=2
 do while(k.le.NN)
  w=v0
  v=ff
  y=v-w

  z=-matmul(AA,y)
  call transpose(s,ss,n)

  p=-matmul(ss,z)
  call transpose(u,uu,n)

  uu=matmul(ss,AA)

  l=(s+z)/p(1,1)
  AA=AA+ (matmul(l,uu))
  s=-matmul(AA,v)
  x=x+s
  if(norm2(s).lt.TOL) then
      print*,x
      stop
  end if
  k=k+1
  print*,"max iterations reached"
  stop
 end do
 end program

 subroutine inverse(a,cc,n)
  implicit none
  integer::i,j,k,n
  real::a(n,n),cc(n,n),coeff
  real::L(n,n),U(n,n),b(n),d(n),x(n)

  L=0.0
  U=0.0
  b=0.0

  do k=1,n-1
  do i=k+1,n
  coeff=a(i,k)/a(k,k)
  L(i,k) = coeff
  do j=k+1,n
  a(i,j) = a(i,j)-coeff*a(k,j)
  end do
 end do
end do

do i=1,n
  L(i,i)=1.0
end do

do j=1,n
  do i=1,j
      U(i,j)=a(i,j)
  end do
end do

do k=1,n
  b(k)=1.0
  d(1)=b(1)
  do i=2,n
      d(i)=b(i)

   do j=1,i-1
      d(i)=d(i)-L(i,j)*d(j)
          end do
end do

x(n)=d(n)/U(n,n)
do i = n-1,1,-1
  x(i) = d(i)
  do j=n,i+1,-1
    x(i)=x(i)-U(i,j)*x(j)
  end do
  x(i) = x(i)/u(i,i)
end do
do i=1,n
  cc(i,k) = x(i)
end do
b(k)=0.0
end do
end subroutine

subroutine transpose(q,r,n)
  implicit none
  integer::i,n
  real::q(n,1),r(1,n)
  do i=1,n
      r(1,i)=0
  end do
      do i=1,n
          r(1,i)= r(i,1) +q(i,1)
      end do
end subroutine













real function f1(x1,x2,x3)
implicit none
real:: x1,x2,x3
f1=3*x1-cos(x2*x3)-(1/2.0)
end function



real function f2(x1,x2,x3)
implicit none
real:: x1,x2,x3
f2=x1*x1-(81*((x2+0.1)**2))+sin(x3)+1.06
end function


real function f3(x1,x2,x3)
implicit none
real:: x1,x2,x3
f3=exp(-x1*x2)+20*x3+(31.4-3)/3.0
end function
