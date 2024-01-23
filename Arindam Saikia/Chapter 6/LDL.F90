program ldl
implicit none
    
integer, parameter::n=3
real a(n,n),l(n,n),d(n),v(n),dum
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

    do i=1,n,1
        do j=1,i-1,1
            v(j)=l(i,j)*d(j)
        enddo
        
        dum=0
        
        do j=1,i-1,1
            dum=dum+l(i,j)*v(j)
        enddo
        
        d(i)=a(i,i)-dum
        
        do j=i+1,n,1
            dum=0
        
            do k=1,i-1,1
                dum=dum+l(j,k)*v(k)
            enddo
        
            l(j,i)=(a(j,i)-dum)/d(i)
        
        enddo
    enddo

    do i=1,n,1
        do j=1,i-1,1
            write(*,*) l(i,j)
        enddo
    enddo
    do i=1,n,1
        write(*,*) d(i)
    enddo

end program ldl