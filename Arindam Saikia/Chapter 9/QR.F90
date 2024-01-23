program qr
    implicit none
    integer :: m,k,n,i,j
    real :: Tol,lamda,lamda1,lamda2,mu1,mu2,sigma,bb,dd,cc,shift
    real,dimension(:),allocatable :: a,b,sig,d,c,z,x,y,q,s,r

    n=3
    Tol=10E-6
    m=100

allocate(a(n))
allocate(b(2:n))
allocate(d(n))
allocate(x(n))
allocate(y(n))
allocate(z(n))
allocate(c(n))
allocate(sig(n))
allocate(s(n))
allocate(q(n))
allocate(r(n))

a=(/3.0,3.0,3.0/)
b=(/1.0,1.0/)
    
k=1
shift=0

do while(k<=m)
    if(abs(b(n))<=Tol) then
        lamda=a(n)+shift
        print*,"lamda ->",lamda
        n=n-1
    end if

    if(abs(b(2))<=Tol) then
        lamda=a(1)+shift
        print*,"lamda ->",lamda
        n=n-1
        a(1)=a(2)
        do j=2,n
            a(j)=a(j+1)
            b(j)=b(j+1)
        end do
    end if

    if (n==0) then
        lamda=a(1)+shift
        write(*,*) "lamda ->",lamda
        stop
    end if

    do j=3,n-1
        if(abs(b(j))<+Tol) then
            print*,"Split Into"
            print*,"a:"

            do i=1,j-1
                print*,a(i)
            enddo
            print*,'b:'
            do i=1,j-1
                print*,b(i)
            enddo
            print*,'and'
            print*,'a:'

            do i=j,n
                print*,a(i)
            enddo
            print*,'b:'
            do i=j,n
                print*,b(i)
            enddo
            stop
        endif
    enddo

    bb=-(a(n-1)+a(n))
        cc=a(n)*a(n-1)-(b(n)**2)
        dd=((bb**2)-4*cc)**(0.5)
        if(bb>0) then
            mu1=(-2*cc)/(bb+dd)
            mu2=-(bb+dd)/2.0
        else
            mu1=(dd-bb)/2.0
            mu2=2*cc/(dd-bb)
        endif
        if(n==2) then
            lamda1=mu1+shift
            lamda2=mu2+shift
            print*, 'lamda1 ->',lamda1,'lamda2 ->',lamda2
            stop
        endif

        if(abs(mu2-a(n))<abs(mu1-a(n))) then
            sigma=mu2
        else
            sigma=mu1
        endif
        shift=shift+sigma
        do j=1,n
            d(j)=a(j)-sigma
        enddo
        x(1)=d(1)
        y(1)=b(2)
        do j=2,n
            z(j-1)=((x(j-1)**2)+(b(j)**2))**0.5
            c(j)=x(j-1)/z(j-1)
            sig(j)=b(j)/z(j-1)
            q(j-1)=c(j)*y(j-1)+sig(j)*d(j)
            x(j)=-sig(j)*y(j-1)+c(j)*d(j)
            if(j.ne.n) then
                r(j-1)=sig(j)*b(j+1)
                y(j)=c(j)*b(j+1)
                endif
            enddo
            z(n)=x(n)
            a(1)=sig(2)*q(1)+c(2)*z(1)
            b(2)=sig(2)*z(2)
            do j=2,n-1
                a(j)=sig(j+1)*q(j)+c(j)*c(j+1)*z(j)
                b(j+1)=sig(j+1)*z(j+1)
            enddo
            a(n)=c(n)*z(n)
            k=k+1

        enddo
        print*,'max iterations exceeded'

end program qr