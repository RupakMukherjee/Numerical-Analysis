program gaussianelipartialpiv

    implicit none
    integer,parameter::n=3
    real dum(n),dums,nrow(n)
    real,dimension(n,n+1)::a
    real,dimension(n)::x
    real,dimension(2:n,n-1)::m
    integer i,j,p
    a(1,1)=4
    a(1,2)=-1
    a(1,3)=1
    a(1,4)=8
    a(2,1)=2
    a(2,2)=5
    a(2,3)=2
    a(2,4)=3
    a(3,1)=1
    a(3,2)=2
    a(3,3)=4
    a(3,4)=11

    do i=1,n,1
        nrow(i)=i
    enddo
    do i=1,n-1,1
        dums=0
        do j=i,n,1
        if(abs(a(nrow(j),i))>dums) then
            p=j

        endif

        dums=abs(a(nrow(j),i))

        enddo
        if(a(nrow(p),i)==0) then
            write(*,*)'no unique solution exists'
            stop
        endif
        if(nrow(i)/=nrow(p)) then
            dum(i)=nrow(i)
            nrow(i)=nrow(p)
            nrow(p)=dum(i)
            endif
            do j=i+1,n,1
                m(nrow(j),i)=a(nrow(j),i)/a(nrow(i),i)
                a(nrow(j),1)=a(nrow(j),1)-m(nrow(j),i)*a(nrow(i),1)
                a(nrow(j),2)=a(nrow(j),2)-m(nrow(j),i)*a(nrow(i),2)
                a(nrow(j),3)=a(nrow(j),3)-m(nrow(j),i)*a(nrow(i),3)
                a(nrow(j),4)=a(nrow(j),4)-m(nrow(j),i)*a(nrow(i),4)
            enddo
        enddo
        if(a(nrow(n),n)==0) then
            write(*,*)'no unique solution exists'
            stop
        endif
        x(n)=a(nrow(n),n+1)/a(nrow(n),n)
        do i=n-1,1,-1
            dums=0
            do j=i+1,n,1
                dums=dums+a(nrow(i),j)*x(j)
            enddo
            x(i)=(a(nrow(i),n+1)-dums)/a(nrow(i),i)
        enddo
        write(*,*) x

    
end program gaussianelipartialpiv
