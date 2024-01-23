program inversepm
implicit none


integer,parameter::n=3
real a(n,n),aa(n,n),x(n,1),tol,nn,mu,er,y(n,1),dum(n),q,dum1(1,n),dum2(1,1),dum3(1,1),xt(1,n)
integer i,j,k,p,stopper

tol=10e-8
nn=10
a(1,1)=-4
a(1,2)=14
a(1,3)=0
a(2,1)=-5
a(2,2)=13
a(2,3)=0
a(3,1)=-1
a(3,2)=0
a(3,3)=2
x(1,1)=1
x(2,1)=1
x(3,1)=1

    do i=1,n
        xt(1,i)=x(i,1)
    enddo

    call mm(xt,a,dum1,1,n,n)
    call mm(dum1,x,dum2,1,n,1)
    call mm(xt,x,dum3,1,n,1)

    q=dum2(1,1)/dum3(1,1)

        do i=1,n
            do j=1,n
            if(i==j) then
            aa(i,j)=a(i,j)-q
            else
            aa(i,j)=a(i,j)
            endif
            enddo
        enddo

    k=1
    do i=1,n,1
        if(abs(x(i,1))==maxval(abs(x(1:n,1)))) then
        p=i
        exit
        endif
    enddo
    do i=1,n,1
        x(i,1)=x(i,1)/x(p,1)
    enddo
    do while(k<=nn)


        stopper=0
        call ge(aa,x,y,n,stopper)
        print*, y
        if(stopper==1) then
            write(*,*) 'q is an eigenvalue,q=',q
            stop
        endif
        mu=y(p,1)
        do i=1,n,1
            if(abs(y(i,1))==maxval(abs(y(1:n,1)))) then
                p=i
            exit
            endif
        enddo

        do i=1,n,1
            dum(i)=x(i,1)-(y(i,1)/y(p,1))
        enddo
        er=maxval(abs(dum(1:n)))
        do i=1,n,1
            x(i,1)=y(i,1)/y(p,1)
        enddo
        if (er<tol) then
            mu=(1.0/mu)+q
            write(*,*) 'mu=',mu,'x=',x
            stop
        endif
        k=k+1
    enddo
    write(*,*)'max no. of iterations exceeded'

end program inversepm


subroutine mm(a,b,c,n,m,p)
    implicit none
    integer n,m,p
    real a(n,m),b(m,p),c(n,p),dum
    integer i,j,k

    do i=1,n,1
        do j=1,p,1
            dum=0
            do k=1,m,1
                dum=dum+a(i,k)*b(k,j)
            enddo
            c(i,j)=dum
        enddo
    enddo

endsubroutine

subroutine ge(a,b,x,n,stopper)
    implicit none
    integer n
    real dum(n,n+1),dums

    real,dimension(n)::e
    real,dimension(n,n)::a
    real,dimension(n,1)::b
    real,dimension(n,1)::x
    real,dimension(2:n,n-1)::m
    integer i,j,p,stopper

    p=0
    do i=1,n-1,1
        do j=i,n,1
        if(a(j,i).ne.0) then
            p=j
            exit
        endif

        enddo

        if(p==0) then
            stopper=1
            stop
        endif

        if(p/=i) then
        dum(p,1)=a(p,1)
        dum(p,2)=a(p,2)
        dum(p,3)=a(p,3)
        dum(p,4)=b(p,1)
        a(p,1)=a(i,1)
        a(p,2)=a(i,2)
        a(p,3)=a(i,3)
        b(p,1)=b(i,1)
        a(i,1)=dum(p,1)
        a(i,2)=dum(p,2)
        a(i,3)=dum(p,3)
        b(i,1)=dum(p,4)
        endif
        do j=i+1,n,1
            m(j,i)=a(j,i)/a(i,i)
            a(j,1)=a(j,1)-m(j,i)*a(i,1)
            a(j,2)=a(j,2)-m(j,i)*a(i,2)
            a(j,3)=a(j,3)-m(j,i)*a(i,3)
            b(j,1)=b(j,1)-m(j,i)*b(i,1)
        enddo
    enddo
    if(a(n,n)==0) then
        stopper=1
        stop
    endif

    x(n,1)=b(n,1)/a(n,n)
    do i=n-1,1,-1
        dums=0
        do j=i+1,n,1
            dums=dums+a(i,j)*x(j,1)
        enddo
        x(i,1)=(b(i,1)-dums)/a(i,i)
    enddo


end subroutine





