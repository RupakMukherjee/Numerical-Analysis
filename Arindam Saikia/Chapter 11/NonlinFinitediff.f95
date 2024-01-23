program nonlinfindif
implicit none

integer,parameter::n=19
real ain,bin,tol,alpha,beta,w(0:n+1),h,t,f,fy,fyy,x,l(n),c(n),a(n),b(n),z(n),u(n),d(n),v(n)
integer i,j,k,m

print*,'endpoint a'
read(*,*) ain

print*,'endpoint b'
read(*,*) bin

print*,'boundary cond. alpha'
read(*,*) alpha

print*,'boundary cond beta'
read(*,*) beta

print*,'tolerance'
read(*,*) tol

print*,'max iterations'
read(*,*) m


h=(bin-ain)/(n+1)
w(0)=alpha
w(n+1)=beta
    do i=1,n,1
        w(i)=alpha+i*h*((beta-alpha)/(bin-ain))
    enddo
    k=1
    do while(k<=m)
        x=ain+h
        t=(w(2)-alpha)/(2*h)
        a(1)=2+(h**2)*fy(x,w(1),t)
        b(1)=-1+(h/2)*fyy(x,w(1),t)
        d(1)=-(2*w(1)-w(2)-alpha+(h**2)*f(x,w(1),t))

        do i=2,n-1,1
            x=ain+i*h
            t=(w(i+1)-w(i-1))/(2*h)
            a(i)=2+(h**2)*fy(x,w(i),t)
            b(i)=-1+(h/2)*fyy(x,w(i),t)
            c(i)=-1-(h/2)*fyy(x,w(i),t)
            d(i)=-(2*w(i)-w(i+1)-w(i-1)+(h**2)*f(x,w(i),t))
        enddo

        x=bin-h
        t=(beta-w(n-1))/(2*h)
        a(n)=2+(h**2)*fy(x,w(n),t)
        c(n)=-1-(h/2)*fyy(x,w(n),t)
        d(n)=-(2*w(n)-w(n-1)-beta+(h**2)*f(x,w(n),t))

        l(1)=a(1)
        u(1)=b(1)/a(1)
        z(1)=d(1)/l(1)

        do i=2,n-1
            l(i)=a(i)-c(i)*u(i-1)
            u(i)=b(i)/l(i)
            z(i)=(d(i)-c(i)*z(i-1))/l(i)
        enddo

        l(n)=a(n)-c(n)*u(n-1)
        z(n)=(d(n)-c(n)*z(n-1))/l(n)

        v(n)=z(n)
        w(n)=w(n)+v(n)

        do i=n-1,1,-1
            v(i)=z(i)-u(i)*v(i+1)
            w(i)=w(i)+v(i)
        enddo

        if(norm2(v)<=tol) then
            do i=0,n+1
                x=ain+i*h
                print*,'x=',x,'w=',w(i)
            enddo
            stop
        endif
        k=k+1
    enddo

    print*, 'max iterations exceeded'

end program nonlinfindif

real function f(x,y,yy)
implicit none
real x,y,yy
f=(32+2*(x**3)-y*yy)/8.0
endfunction
real function fy(x,y,yy)
implicit none
real x,y,yy
fy=-yy/8.0
endfunction
real function fyy(x,y,yy)
implicit none
real x,y,yy
fyy=-y/8.0
endfunction
