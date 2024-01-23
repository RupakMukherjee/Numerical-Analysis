program wielandt
implicit none
integer,parameter::n=3
real a(n,n),lamda,v(n,1),x(n-1,1),tol,nn,mu,u(n,1),b(n-1,n-1),wt(n-1),w(n,1),dum
integer i,j,k


    print*,'dimension of matrix a,(n), must be assigned in the program'

    do i=1,n,1
        do j=1,n,1
            write(*,*)i,',',j,'element of a'
            read(*,*) a(i,j)
        enddo
    enddo

    write(*,*) 'approx eigenvalue(lamda)'
    read*, lamda
    do i=1,n,1
    print*, i,1,'component of eigenvector'
    read*, v(i,1)
    enddo
    do i=1,n-1
    print*, i,1,'component of vector'
    read*, x(i,1)
    enddo
    print*,'tolerance'
    read*,tol
    print*,'max iteration'
    read*, nn

    do j=1,n
        if(abs(v(i,1))==maxval(abs(v(1:n,1)))) then
        i=j
        exit
        endif
    enddo
    if(i.ne.1) then
        do k=1,i-1
            do j=1,i-1
                b(k,j)=a(k,j)-a(i,j)*v(k,1)/v(i,1)
            enddo
        enddo
    endif
    if(i.ne.1 .and. i.ne.n) then
        do k=i,n-1
            do j=1,i-1
                b(k,j)=a(k+1,j)-a(i,j)*v(k+1,1)/v(i,1)
                b(j,k)=a(j,k+1)-a(i,k+1)*v(j,1)/v(i,1)
            enddo
        enddo
    endif
    if(i.ne.n) then
        do k=i,n-1
            do j=i,n-1
                b(k,j)=a(k+1,j+1)-a(i,j+1)*v(k+1,1)/v(i,1)
            enddo
        enddo
    endif
    call pm(b,x,n-1,mu,wt)
    if(i.ne.1) then
        do k=1,i-1
            w(k,1)=wt(k)
        enddo
    endif
    w(i,1)=0
    if(i.ne.n) then
        do k=i+1,n
            w(k,1)=wt(k-1)
        enddo
    endif
    do k=1,n
        dum=0
        do j=1,n
            dum=dum+a(i,j)*w(j,1)
        enddo
        u(k,1)=(mu-lamda)*w(k,1)+dum*v(k,1)/v(i,1)
    enddo
    print*,'mu=',mu,'u=',u
    stop
end program wielandt



subroutine pm(a,x0,n,mu,x)
    implicit none
    integer n
    real a(n,n),x0(n,1),x(n),tol,nn,mu,er,y(n,1),dum(n)
    integer i,j,k,p

    tol=10e-5
    nn=100


    k=1
    do i=1,n,1
        if(abs(x0(i,1))==maxval(abs(x0(1:n,1)))) then
        p=i
        exit
        endif
    enddo
    do i=1,n,1
        x0(i,1)=x0(i,1)/x0(p,1)
    enddo
    do while(k<=nn)
        call mm(a,x0,y,n,n,1)
        mu=y(p,1)
        do i=1,n,1
            if(abs(y(i,1))==maxval(abs(y(1:n,1)))) then
                p=i
            exit
            endif
        enddo
        if(y(p,1)==0) then
            write(*,*) 'eigenvector',x0
            write(*,*) 'a has the eigenvalue 0, select a new vector x and restart'
            stop
        endif
        do i=1,n,1
            dum(i)=x0(i,1)-(y(i,1)/y(p,1))
        enddo
        er=maxval(abs(dum(1:n)))
        do i=1,n,1
            x0(i,1)=y(i,1)/y(p,1)
            x(i)=x0(i,1)
        enddo
        if (er<tol) then

            exit
        endif
        k=k+1
    enddo
    if(er>tol) then
    write(*,*)'power method failed'
    stop
    endif

end subroutine
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

end subroutine














