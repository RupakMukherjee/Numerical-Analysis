program poissonfinite
    implicit none
    integer,parameter::n=6,m=5
    real::a,b,c,d,tol,x(n),y(m),w(n,m),lambda,mu,z,norm,h,k
    real,external::f,g
    integer::Nt,l,i,j


   tol=10E-10
   nt=10000
    a=0.0
    b=2.0
    c=0.0
    d=1.0



    h=(b-a)/n
    k=(d-c)/m

    do i=1,n-1,1
        x(i)=a+i*h
    end do

    do j=1,m-1,1
        y(j)=c+j*k
    end do

    do i=1,n-1,1
        do j=1,m-1,1
            w(i,j)=0
        end do
    end do

   lambda=(h**2.0)/(k**2.0)
   mu=2*(1+lambda)
   l=1

   do while(l<=Nt)

    z=((-h**2)*f(x(1),y(m-1))+g(a,y(m-1))+lambda*g(x(1),d)+lambda*w(1,m-2)+w(2,m-1))/mu
    norm=abs(z-w(1,m-1))
    w(1,m-1)=z

    do i=2,n-2,1
 z=((-h**2)*f(x(i),y(m-1))+lambda*g(x(i),d)+w(i-1,m-1)+w(i+1,m-1)+lambda*w(i,m-2))/mu

 if(abs(w(i,m-1)-z)>Norm) then
    Norm=abs(w(i,m-1)-z)
    w(i,m-1)=z
 end if

 end do

 z=((-h**2)*f(x(n-1),y(m-1))+g(b,y(m-1))+lambda*g(x(n-1),d)+w(n-2,m-1)+lambda*w(n-1,m-2))/mu
if(abs(w(n-1,m-1)-z)>Norm) then
    norm=abs(w(n-1,m-1)-z)
    w(n-1,m-1)=z
end if

do j=m-2,2,-1

z=((-h**2)*f(x(1),y(j))+g(a,y(j))+lambda*w(1,j+1)+lambda*w(1,j-1)+w(2,j))/mu
    if(abs(w(1,j)-z)>Norm) then
        Norm=abs(w(1,j)-z)
        w(1,j)=z
    end if

 do i=2,n-2,1


 z=((-h**2)*f(x(i),y(j))+w(i-1,j)+lambda*w(i,j+1)+w(i+1,j)+lambda*w(i,j-1))/mu
 if (abs(w(i,j)-z)>Norm) then
    norm=abs(w(i,j)-z)
    w(i,j)=z
 end if
end do

z=((-h**2)*f(x(n-1),y(j))+g(b,y(j))+w(n-2,j)+lambda*w(n-1,j+1)+lambda*w(n-1,j-1))/mu
if(abs(w(n-1,j)-z)>Norm) then
    Norm=abs(w(n-1,j)-z)
    w(n-1,j)=z
end if

end do !end at step 13

z=((-h**2)*f(x(1),y(1))+g(a,y(1))+lambda*g(x(1),c)+lambda*w(1,2)+w(2,1))/mu
 if(abs(w(1,1)-z)>Norm) then
    norm=abs(w(1,1)-z)
    w(1,1)=z
 end if

 do i=2,n-2,1
  z=((-h**2)*f(x(i),y(1))+lambda*g(x(i),c)+w(i-1,1)+lambda*w(i,2)+w(i+1,1))/mu
  if(abs(w(i,1)-z)>norm) then
    Norm=abs(w(i,1)-z)
    w(i,1)=z
  end if
  end do

  z=((-h**2)*f(x(n-1),y(1))+g(b,y(1))+lambda*g(x(n-1),c)+w(n-2,1)+lambda*w(n-1,2))/mu
  if(abs(w(n-1,1)-z)>norm) then
    norm=abs(w(n-1,1)-z)
    w(n-1,1)=z
  end if
  !Step 17
  if(Norm<=tol) then
    do i=1,n-1,1
        do j=1,m-1,1
          !  print*,"x=",x(i)
           ! print*,"y=",y(j)
            print*,"i=",i,"j=",j,"w=",w(i,j),'x=',x(i),"y=",y(j)
           ! stop
        end do
    end do
    stop
  end if
  l=l+1
   end do !ends at step 20
   print*,"maximum number of iterations exceeded"
end program

real function f(x,y)
real::x,y
f=x*exp(y)
return
end function

real function g(x,y)
real::x,y
if(x==a) then
    g=a
else if(x==b) then
    g=2*exp(y)
else if(y==c) then
    g=x
else if(y==d) then
    g=exp(1.0)*x
end if
end function


