program lufact
    implicit none
    integer,parameter :: n=3
    integer :: i,j,k
    real :: dum,dums
    real,dimension(1:n,1:n) :: a,l,u

    a(1,1)=2
    a(1,2)=-1
    a(1,3)=1
    a(2,1)=3
    a(2,2)=3
    a(2,3)=9
    a(3,1)=3
    a(3,2)=3
    a(3,3)=5
    l(1,1)=1
    l(2,2)=1
    l(3,3)=1


    u(1,1)=a(1,1)/l(1,1)

    if(l(1,1)*u(1,1)==0) then
        write(*,*) 'factorization impossible'
        stop
    endif
    

    do j=2,n,1
        u(1,j)=a(1,j)/l(1,1)
        l(j,1)=a(j,1)/u(1,1)
    end do

    do i=2,n-1,1
        dum=0
        do j=1,i-1,1
            dum=dum+l(i,j)*u(j,i)
        enddo
        u(i,i)=(a(i,i)-dum)/l(i,i)
            if(l(i,i)*u(i,i)==0) then
                write(*,*) 'factorization impossible'
                stop
            endif

            do j=i+1,n,1
                dum=0
                dums=0
                do k=1,i-1,1
                    dum=dum+l(i,k)*u(k,j)
                    dums=dums+l(j,k)*u(k,i)
                enddo
                u(i,j)=(1/l(i,i))*(a(i,j)-dum)
                l(j,i)=(1/u(i,i))*(a(j,i)-dums)
            enddo

    enddo

    dum=0
    do k=1,n-1,1
        dum=dum+l(n,k)*u(k,n)
    enddo
    u(n,n)=(a(n,n)-dum)/l(n,n)
    do i=1,n,1
        do j=1,n,1
            print"(2f8.2)", l(i,j)
        enddo
        write(*,*) '_________'
    enddo
    write(*,*) '==========='
    do i=1,n,1
        do j=1,n,1
            print"(2f8.2)", u(i,j)
        enddo
        write(*,*) '_________'
    enddo

end program lufact