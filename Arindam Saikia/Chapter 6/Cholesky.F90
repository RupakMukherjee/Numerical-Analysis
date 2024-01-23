program Chol
    implicit none

    integer, parameter::n=3
    real a(n,n),l(n,n),dum
    integer i,j,k

    a(1,1)=4
    a(1,2)=-1
    a(1,3)=1
    a(2,1)=-1
    a(2,2)=4.25
    a(2,3)=2.75
    a(3,1)=1
    a(3,2)=2.75
    a(3,3)=3.5
    
    l(1,1)=((a(1,1)))**(1/2.0)

    do j=2,n,1
    l(j,1)=a(j,1)/l(1,1)
    enddo

    do i=2,n-1,1
    dum=0
        do k=1,i-1,1
            dum=dum+(l(i,k))**2
        enddo

        l(i,i)=((a(i,i)-dum))**(1.0/2)

        do j=i+1,n,1
            dum=0
            do k=1,i-1,1
                dum=dum+l(j,k)*l(i,k)
            enddo
            l(j,i)=(a(j,i)-dum)/l(i,i)
        enddo
    enddo

    dum=0
    do k=1,n-1,1
        dum=dum+(l(n,k))**2
    enddo

    print*, "L Matrix:"
    l(n,n)=sqrt(a(n,n)-dum)
    do i=1,n,1
        do j=1,n,1
        write(*,*) l(j,i)
    enddo
    print*,'. . . . . . . . .'
    enddo

    print*, "L' Matrix:"
    l(n,n)=sqrt(a(n,n)-dum)
    do i=1,n,1
        do j=1,n,1
        write(*,*) l(i,j)
    enddo
    print*,'. . . . . . . . .'
    enddo
end program Chol