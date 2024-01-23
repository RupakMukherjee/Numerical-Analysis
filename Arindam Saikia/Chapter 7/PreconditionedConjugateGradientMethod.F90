program preconditionedconjugategradient
implicit none

integer::i,j,k,Q
integer,parameter::n=5
real,dimension(1:n)::r,x,w,v,b,u,ax,av
real,dimension(1:n,1:n)::C,a
real::alpha,sum0,t,s,dum,beta,tol


Q=14
TOL=10E-7

 a(1,1)=0.2
 a(1,2)=0.1
 a(1,3)=1
 a(1,4)=1
 a(1,5)=0

 a(2,1)=0.1
 a(2,2)=4
 a(2,3)=-1
 a(2,4)=1
 a(2,5)=-1

 a(3,1)=1
 a(3,2)=-1
 a(3,3)=60
 a(3,4)=0
 a(3,5)=-2

 a(4,1)=1
 a(4,2)=1
 a(4,3)=0
 a(4,4)=8
 a(4,5)=4

 a(5,1)=0
 a(5,2)=-1
 a(5,3)=-2
 a(5,4)=4
 a(5,5)=700

 b(1)=1
 b(2)=2
 b(3)=3
 b(4)=4
 b(5)=5

 x(1)=0
 x(2)=0
 x(3)=0
 x(4)=0
 x(5)=0

 C(1,1)=1
 C(1,2)=0
 C(1,3)=0
 C(1,4)=0
 C(1,5)=0

 C(2,1)=0
 C(2,2)=1
 C(2,3)=0
 C(2,4)=0
 C(2,5)=0

 C(3,1)=0
 C(3,2)=0
 C(3,3)=1
 C(3,4)=0
 C(3,5)=0

 c(4,1)=0
 c(4,2)=0
 c(4,3)=0
 c(4,4)=1
 c(4,5)=0

 c(5,1)=0
 c(5,2)=0
 c(5,3)=0
 c(5,4)=0
 c(5,5)=1

r=b-matmul(a,x)
w=matmul(C,r)
v=matmul(C,w)
alpha=0
do j=1,n
alpha=alpha +w(j)**2
end do

k=1
do while(k.le.Q)
    if(norm2(v(1:n)).lt.TOL) then
        print*,'solution vector:',x(1:n)
        print*,"residue vector:",r(1:n)
        stop
   end if


u=matmul(a,v)
dum=0
do j=1,n
    dum=dum+v(j)*u(j)
end do
t=alpha/dum
x=x+t*v
r=r-t*u
w=matmul(C,r)

beta=0
do j=1,n
    beta=beta+w(j)**2
    end do


if(abs(beta).lt.TOL) then
    if((norm2(r(1:n))).lt.TOL) then
        print*,'sol vector',x(1:n)
        print*,'residue vector',r(1:n)
        stop
    end if
end if

s=(beta/alpha)
v=matmul(C,w)+s*v
alpha=beta
k=k+1
end do


    !step 8
    if (k > n) then 
        print*, "The Maximum number of iterations was exceeded."
        stop
    end if

end program preconditionedconjugategradient
