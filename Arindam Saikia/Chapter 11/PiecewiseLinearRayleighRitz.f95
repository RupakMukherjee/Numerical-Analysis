program piece
    implicit none
    integer,parameter::n=9
    real x(0:n+1),c(n),dum,q(6,n+1),alpha(n),beta(n-1),z(n),zeta(n),f,h,phi,b(n),a(n),ain,bin,ii,inc
    integer i,j,k

    print*,'endpoint a'
    read(*,*) ain
    print*,'endpoint b'
    read(*,*) bin
    h=(bin-ain)/(1.0*n+1)

    do i=0,n+1,1

        x(i)=h*i
    enddo




    do i=1,n+1
        if(i<=n) then
        if(i<=n-1) then
        call cs(x(i),x(i+1),q(1,i),1)
        endif
        q(1,i)=((1/h)**2)*q(1,i)
        call cs(x(i-1),x(i),q(2,i),2)
        q(2,i)=((1/h)**2)*q(2,i)
        call cs(x(i),x(i+1),q(3,i),3)
        q(3,i)=((1/h)**2)*q(3,i)

        call cs(x(i-1),x(i),q(5,i),5)
        q(5,i)=((1/h))*q(5,i)
        call cs(x(i),x(i+1),q(6,i),6)
        q(6,i)=((1/h))*q(6,i)
        endif
        call cs(x(i-1),x(i),q(4,i),4)
        q(4,i)=((1/h)**2)*q(4,i)
    enddo



    do i=1,n-1
        alpha(i)=q(4,i)+q(4,i+1)+q(2,i)+q(3,i)
        beta(i)=q(1,i)-q(4,i+1)
        b(i)=q(5,i)+q(6,i)
    enddo

    alpha(n)=q(4,n)+q(4,n+1)+q(2,n)+q(3,n)
    b(n)=q(5,n)+q(6,n)



    a(1)=alpha(1)
    zeta(1)=beta(1)/alpha(1)
    z(1)=b(1)/a(1)

    do i=2,n-1,1
        a(i)=alpha(i)-beta(i-1)*zeta(i-1)

        zeta(i)=beta(i)/a(i)

        z(i)=(b(i)-beta(i-1)*z(i-1))/a(i)
    enddo
    a(n)=alpha(n)-beta(n-1)*zeta(n-1)
    z(n)=(b(n)-beta(n-1)*z(n-1))/a(n)

    c(n)=z(n)
    print*,x(n),c(n)

    do i=n-1,1,-1
        c(i)=z(i)-zeta(i)*c(i+1)
        print*,x(i),c(i)
    enddo


endprogram


subroutine cs(a,b,xi,hi)

    implicit none

    real XI0,XI1,XI2,X,h,ff,a,b,xi
    integer i,nnn,hi

    nnn=10

    h=(b-a)/nnn

    XI0=ff(a,b,a,hi)+ff(a,b,b,hi)
    XI1=0
    XI2=0

    do i=1,nnn-1,1
    X=a+i*h

        if(mod(i,2)==0)then
            XI2=XI2+ff(a,b,X,hi)
        else
            XI1=XI1+ff(a,b,X,hi)
        end if
    end do

    XI=h*(XI0+2*XI2+4*XI1)/3.0




end subroutine

real function ff(a,b,x,hi)
implicit none
real x,a,b,sq,p,f
integer hi
if(hi==1) then
ff=(b-x)*(x-a)*sq(x)
endif
if(hi==2) then
ff=((x-a)**2)*sq(x)
endif
if(hi==3) then
ff=((b-x)**2)*sq(x)
endif
if(hi==4) then
ff=p(x)
endif
if(hi==5) then
ff=(x-a)*f(x)
endif
if(hi==6) then
ff=(b-x)*f(x)
endif
end function

real function sq(x)
implicit none
real x
sq=3.14**2
end function
real function f(x)
implicit none
real x
f=2*(3.14**2)*sin(3.14*x)
end function
real function p(x)
implicit none
real x
p=1.0
end function

real function phi(i,xx,x,ain,bin,h,n)
implicit none
real ain,bin,x(0:n),xx,h
integer i,n
if(xx.ge.ain .and. xx.le.x(i-1)) then
    phi=0.0
endif
if(xx>x(i-1) .and. xx.le.x(i)) then
    phi=(xx-x(i-1))/h
endif
if(xx>x(i) .and. xx.le.x(i+1)) then
    phi=(x(i+1)-xx)/h
endif
if(xx>x(i+1) .and. xx.le.bin) then
    phi=0.0
endif
endfunction
