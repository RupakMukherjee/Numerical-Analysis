program continuationalgo
    implicit none
    
    integer,parameter:: n=3
    real x(n,1),h,f,ff(n,1),b(n,1),k1(n,1),k2(n,1),k3(n,1),k4(n,1),Ja(n,n),dum(n,1)
    integer nn,i,j
    
    nn=4
    x(1,1)=0.0
    x(2,1)=0.0
    x(3,1)=0.0
    
    h=0.25
    do i=1,n
    ff(i,1)=f(i,x(1,1),x(2,1),x(3,1))
    enddo
    
    do i=1,n
    b(i,1)=-h*ff(i,1)
    enddo
    
    
    do i=1,nn
        do j=1,n
            dum(j,1)=x(j,1)
        enddo
    
        ja(1,1)=3.0
        ja(1,2)=dum(3,1)*sin(dum(2,1)*dum(3,1))
        ja(1,3)=dum(2,1)*sin(dum(3,1)*dum(2,1))
        ja(2,1)=2*dum(1,1)
        ja(2,2)=-162.0*(dum(2,1)+0.1)
        ja(2,3)=cos(dum(3,1))
        ja(3,1)=-dum(2,1)*exp(-dum(1,1)*dum(2,1))
        ja(3,2)=-dum(1,1)*exp(-dum(1,1)*dum(2,1))
        ja(3,3)=20.0
    
        call ge(ja,b,k1,n)
        do j=1,n
        dum(j,1)=x(j,1)+0.5*k1(j,1)
        enddo
        ja(1,1)=3.0
        ja(1,2)=dum(3,1)*sin(dum(2,1)*dum(3,1))
        ja(1,3)=dum(2,1)*sin(dum(3,1)*dum(2,1))
        ja(2,1)=2*dum(1,1)
        ja(2,2)=-162*(dum(2,1)+0.1)
        ja(2,3)=cos(dum(3,1))
        ja(3,1)=-dum(2,1)*exp(-dum(1,1)*dum(2,1))
        ja(3,2)=-dum(1,1)*exp(-dum(1,1)*dum(2,1))
        ja(3,3)=20.0
    
        call ge(ja,b,k2,n)
    
        do j=1,n
        dum(j,1)=x(j,1)+0.5*k2(j,1)
        enddo
        ja(1,1)=3.0
        ja(1,2)=dum(3,1)*sin(dum(2,1)*dum(3,1))
        ja(1,3)=dum(2,1)*sin(dum(3,1)*dum(2,1))
        ja(2,1)=2*dum(1,1)
        ja(2,2)=-162*(dum(2,1)+0.1)
        ja(2,3)=cos(dum(3,1))
        ja(3,1)=-dum(2,1)*exp(-dum(1,1)*dum(2,1))
        ja(3,2)=-dum(1,1)*exp(-dum(1,1)*dum(2,1))
        ja(3,3)=20.0
    
        call ge(ja,b,k3,n)
        do j=1,n
        dum(j,1)=x(j,1)+k3(j,1)
        enddo
        ja(1,1)=3.0
        ja(1,2)=dum(3,1)*sin(dum(2,1)*dum(3,1))
        ja(1,3)=dum(2,1)*sin(dum(3,1)*dum(2,1))
        ja(2,1)=2*dum(1,1)
        ja(2,2)=-162*(dum(2,1)+0.1)
        ja(2,3)=cos(dum(3,1))
        ja(3,1)=-dum(2,1)*exp(-dum(1,1)*dum(2,1))
        ja(3,2)=-dum(1,1)*exp(-dum(1,1)*dum(2,1))
        ja(3,3)=20.0
    
        call ge(ja,b,k4,n)
        do j=1,n
        x(j,1)=x(j,1)+(k1(j,1)+2*k2(j,1)+2*k3(j,1)+k4(j,1))/6.0
        enddo
        enddo
    write(*,*) x
    end program
    
    subroutine ge(a,b,x,n)
        implicit none
        integer n
        real dum(n,n+1),dums
    
        real,dimension(n)::e
        real,dimension(n,n)::a
        real,dimension(n,1)::b
        real,dimension(n,1)::x
        real,dimension(2:n,n-1)::m
        integer i,j,p
    
    
    
    
    
        p=0
        do i=1,n-1,1
            do j=i,n,1
            if(a(j,i).ne.0) then
                p=j
                exit
            endif
    
            enddo
    
            if(p==0) then
                print*,'error in solving (ge)'
                exit
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
            print*,'error in solving (ge)'
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
    
    real function f(i,x1,x2,x3)
    implicit none
    real x1,x2,x3
    integer i
    if(i==1) then
    f=3*x1-cos(x2*x3)-0.5
    elseif(i==2) then
    f=x1**2-81*(x2+0.1)**2+sin(x3)+1.06
    elseif(i==3) then
    f=exp(-x1*x2)+20*x3+(10*3.1415926-3)/3.0
    endif
    end function
    