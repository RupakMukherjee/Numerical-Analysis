program paderational
    implicit none

    integer i,j,k
    integer,parameter::m=2,n=3,NN=m+n
    real q(0:m),p(0:n),a(0:NN),fact,b(NN,NN+1),dum,xm
    

    
    do i=0,nn
        a(i)=(-1)**i/fact(i)
    enddo

    q(0)=1
    p(0)=a(0)
    do i=1,NN,1
        do j=1,i-1
            if(j<=n) then
                b(i,j)=0
            endif
        enddo
        if(i<=n) then
            b(i,i)=1
        endif
        do j=i+1,NN,1
            b(i,j)=0
        enddo
        do j=1,i,1
            if(j<=m) then
                b(i,n+j)=-a(i-j)
            endif
        enddo
        do j=n+i+1,NN,1
            b(i,j)=0
        enddo
        b(i,NN+1)=a(i)



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
    do i=n,1,-1
        dum=0
        do j=n+1,nn,1
            dum=dum+b(i,j)*q(j-n)
        enddo
        p(i)=(b(i,nn    +1)-dum)
    enddo
    write(*,*) 'q',q,'p',p

    end program paderational

    real function fact(i)
    implicit none
    integer i,j
    if(i==0) then
        fact=1
    else
        fact=1
        do j=1,i,1
            fact=j*fact
        enddo
    endif
end function
