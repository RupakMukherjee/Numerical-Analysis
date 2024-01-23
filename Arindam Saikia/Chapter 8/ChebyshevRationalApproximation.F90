program cheby
    implicit none
    integer,parameter::m=2,n=3,NN=m+n
    real q(0:m),p(0:n),a(0:nn+m),b(0:nn,0:nn+1),xm,dum,r,x
    integer i,k,j

    do k=0,nn+m,1
        call cs(0.0,3.14,k,a(k))
        a(k)=2*a(k)/3.14

    enddo
    print *, a

    q(0)=1

    do i=0,NN,1
        do j=0,i
            if(j<=n) then
                b(i,j)=0
            endif
        enddo
        if(i<=n) then
            b(i,i)=1
        endif
        do j=i+1,n,1
            b(i,j)=0
        enddo
        do j=n+1,nn,1
            if(i.ne.0) then
                b(i,j)=-(0.5)*(a(i+j-n)+a(abs(i-j+n)))
            else
                b(i,j)=-(0.5)*a(j-n)
            endif
        enddo
        if(i.ne.0) then
            b(i,nn+1)=a(i)
        else
            b(i,nn+1)=(0.5)*a(i)
        endif



    enddo
    do i=n+1,nn-1,1
        do j=i,nn,1
            if(abs(b(j,i))==maxval(abs(b(i:nn,i)))) then
            k=j
            exit
            endif
        enddo
        if(b(k,i)==0) then
            write(*,*) ' the system is singular'
            stop
        endif
        if(k.ne. i) then
            do j=i,nn+1
                dum=b(i,j)
                b(i,j)=b(k,j)
                b(k,j)=dum
            enddo
        endif
        do j=i+1,nn
            xm=b(j,i)/b(i,i)
            do k=i+1,nn+1
                b(j,k)=b(j,k)-xm*b(i,k)
            enddo
            b(j,i)=0
        enddo
    enddo
    if(b(nn,nn)==0) then
        write(*,*)' the system is singular'
        stop
    endif
    if(m>0) then
        q(m)=b(nn,nn+1)/b(nn,nn)
    endif
    do i=nn-1,n+1,-1
        dum=0
        do j=i+1,nn,1
            dum=dum+b(i,j)*q(j-n)
        enddo
        q(i-n)=(b(i,nn+1)-dum)/b(i,i)
    enddo
    do i=n,0,-1
        dum=0
        do j=n+1,nn,1
            dum=dum+b(i,j)*q(j-n)
        enddo
        p(i)=(b(i,nn+1)-dum)
    enddo
    write(*,*) 'q',q,'p',p
    do i=1,5,1
        x=0.2*i
    r=(p(0)-p(1)*x+p(2)*x**2-p(3)*x**3)/(q(0)+q(1)*x+q(2)*x**2)
    write(*,*) x,r

    enddo
    print *, 'maybe huge error because of composite simpson integral'

end program cheby
    subroutine cs(a,b,k,xi)

    implicit none

    real XI0,XI1,XI2,X,h,f,a,b,xi
    integer i,nnn,k

    nnn=6 !use 6,12,18...

    h=(b-a)/nnn

    XI0=f(a,k)+f(b,k)
    XI1=0
    XI2=0

    do i=1,nnn-1,1
    X=a+i*h

        if(mod(i,2)==0)then
            XI2=XI2+f(X,k)
        else
            XI1=XI1+f(X,k)
        end if
    end do

    XI=h*(XI0+2*XI2+4*XI1)/3




end subroutine

real function f(x,k)
real x,g
integer k
f=g(cos(x))*cos(k*x)
end function
real function g(x)
real x
g=exp(-x)
end function




