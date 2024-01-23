program crout
    implicit none

    integer,parameter::n=4
    real,dimension(n,n+1)::a
    real,dimension(n)::x,z
    real,dimension(n,n)::l,u
    integer i

    a(1,1)=2
    a(1,2)=-1
    a(1,3)=0
    a(1,4)=0
    a(1,5)=1
    a(2,1)=-1
    a(2,2)=2
    a(2,3)=-1
    a(2,4)=0
    a(2,5)=0
    a(3,1)=0
    a(3,2)=-1
    a(3,3)=2
    a(3,4)=-1
    a(3,5)=0
    a(4,1)=0
    a(4,2)=0
    a(4,3)=-1
    a(4,4)=2
    a(4,5)=1

    !Step 1:
    l(1,1)=a(1,1)
    u(1,2)=a(1,2)/l(1,1)
    z(1)=a(1,n+1)/l(1,1)

    !Step 2:
    do i=2,n-1,1
        l(i,i-1)=a(i,i-1)
        l(i,i)=a(i,i)-l(i,i-1)*u(i-1,i)
        u(i,i+1)=a(i,i+1)/l(i,i)
        z(i)=(a(i,n+1)-l(i,i-1)*z(i-1))/(l(i,i))
    enddo

    !Step 3:
    l(n,n-1)=a(n,n-1)
    l(n,n)=a(n,n)-l(n,n-1)*u(n-1,n)
    z(n)=(a(n,n+1)-l(n,n-1)*z(n-1))/l(n,n)

    !Step 4:
    x(n)=z(n)

    !Step 5:
    do i=n-1,1,-1
        x(i)=z(i)-u(i,i+1)*x(i+1)
    enddo

    !Step 6:
    write(*,*) x
    
end program crout