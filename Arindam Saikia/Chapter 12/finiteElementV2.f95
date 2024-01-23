program finiteelement
    implicit none
    
    integer, parameter :: n=2 ,m=2, mm=10,nn=10
    integer :: i, j, k, l, t, n, m
    real :: a(10,3), b(10,3), c(10,3), x(11), y(11), delta(11), ntr(10,3), gamma(11)
    real :: alpha(5,6), beta(10), z(10,3,3), H(10,3), jj(10,3,3), II(10,3)
    real ::  jj(10,3,3), ii(10,3), c(10) ,T , g(m,n) ,Ex ,Ey ,dum(4,2)

       
    !Step 0:
    do i=1,mm
        print*,'vertices of T',i
        read(*,*) x(1,i),y(1,i),x(2,i),y(2,i),x(3,i),y(3,i)
    enddo
    do j=1,m
        print*,'node',j
        read(*,*) xx(j),yy(j)
    enddo


    !Step 1:
    do l = n + 1, m
        gamma(l) = g(x(l), y(l))
    end do
    
    !Step 2:
    do i = 1, n
        beta(i) = 0.0
        do j = 1, n
            alpha(i,j) = 0.0
        end do
    end do
    
    !Step 3:
    do i = 1, mm

        delta(i) = x(i,2) * y(i,3) - x(i,3) * y(i,2) - x(i,1) * (y(i,3) - y(i,2)) + y(i,1) * (x(i,3) - x(i,2))
        
        a(i,1) = (x(i,2) * y(i,3) - y(i,2) * x(i,3)) / delta(i)
        b(i,1) = (y(i,2) - y(i,3)) / delta(i)
        c(i,1) = (x(i,3) - x(i,2)) / delta(i)
        a(i,2) = (x(i,3) * y(i,1) - y(i,3) * x(i,1)) / delta(i)
        b(i,2) = (y(i,3) - y(i,1)) / delta(i)
        c(i,2) = (x(i,1) - x(i,3)) / delta(i)
        a(i,3) = (x(i,1) * y(i,2) - y(i,1) * x(i,2)) / delta(i)
        b(i,3) = (y(i,1) - y(i,2)) / delta(i)
        c(i,3) = (x(i,2) - x(i,1)) / delta(i)
        
        ! Add as external function -
        !do j = 1,3,1
        !    NN(i,j) = a(i,j) + b(i,j) * x + c(i,j) * y
        !end do
    end do

    !Step 4:

    do i = 1, MM
        do j = 1,3,1
            do k = 1, j

                call gdi(1,dum(1,2),x(1,i),x(3,i),y(3,i),y(2,i),i,j,k,a,b,c,mm)
                call gdi(2,dum(2,2),x(1,i),x(3,i),y(3,i),y(2,i),i,j,k,a,b,c,mm)
                call gdi(3,dum(3,2),x(1,i),x(3,i),y(3,i),y(2,i),i,j,k,a,b,c,mm)
                
                z(i, j, k) = b(i, j) * b(i, k) * dum(1,2) + c(i, j) * c(i, k) * dum(2,2) - dum(3,2)
                
                call gdi(4,dum(4,2),x(1,i),x(3,i),y(3,i),y(2,i),i,j,k,a,b,c,mm)
                H(i, j) = -dum(4,2)
            end do
        end do
    end do
    
    !Step 5:
    do i = k + 1, NN
        do j = 1,3,1
            do k = 1, j
                !JJ(i, j, k) = i5
                !II(i, j) = i6

                call gdi(5,jj(j,k,i),Ex(2),Ex(4),Ey(2),Ey(4),j,k,i,a,b,c,mm)
                call gdi(6,ii(j,i),Ex(2),Ex(4),Ey(2),Ey(4),j,k,i,a,b,c,mm)
            end do
        end do
    end do
    
    !Step 6:
    do i = 1, MM
        
        !Step 7:
        do k = 1,3,1
        
            !Step 8:
            l = 1
            do while (Ex(l) .ne. x(i,k) .and. Ey(l) .ne. y(i,k))
            Ex(l) =x(i,l)
            Ey(l) =y(i,l)

                if (Ex(l) == x(i,k) .and. Ey(l) == y(i,k)) goto 10 !then
                !exit
                !else
    
                !end if
                l = l + 1
            end do 
            
      10      !Step 9:
            if (k > 1) then
                do j = 1, k - 1
                    
                    !Step 10:
                    t = 1
                    do while (Ex(t) .ne. x(i,j) .and. Ey(t) .ne. y(i,j))
                    Ex(t) =x(i,t)
                    Ey(t) =y(i,t)

                if (Ex(t) == x(i,j) .and. Ey(t) == y(i,j)) goto 11 !then
                !exit
                !else
                t = t + 1
                !end if
            end do
                

                11    !Step 11:
                    if (l <= n) then
                        if (t <= n) then
                            alpha(l, t) = alpha(l, t) + z(i, k, j)
                            alpha(t, l) = alpha(t, l) + z(i, k, j)
                        else
                            beta(l) = beta(l) - gamma(t) * z(i, k, j)
                        
                        end if
                    else 
                        beta(t)=beta(t)-gamma(l)*z(i, k ,j)
                    end if
                end do
            end if
    
    
    
            !Step 12:
            if (l<=n) then
            alpha(l,l)=alpha(l,l)+z(i,k,k)
            beta(l)=beta(l)+H(i,k)
        end do
    end do

    !Step 13:
        do i=k+1,NN
            
            !Step 14:
            do k = 1,3,1
        
            !Step 15:
            l = 1
            do while (Ex(l) .ne. x(i,k) .and. Ey(l) .ne. y(i,k))
            Ex(l) =x(i,l)
            Ey(l) =y(i,l)

                if (E(l) == x(i,k) .and. Ey(l) == y(i,k)) then
                exit
                else
                l = l + 1
                end if
            end do 
            
            
            !Step 16:
            if (k > 1) then
                do j = 1, k - 1,1
                    
                    !Step 17:
                    t = 1
                    do while (Ex(t) .ne. x(i,k) .and. Ey(t) .ne. y(i,k)))
                    Ex(t) = x(i,t)
                    Ey(t) = y(i,t)

                      if (Ex(t) == x(i,k) .and. Ey(t)== y(i,k)) goto 18 !then
                       !exit
                        !else
                    t = t + 1
                      !end if
                    end do 
                

                    !Step 18:
                18  if (l <= n) then
                        if (t <= n) then
                            alpha(l, t) = alpha(l, t) + jj(i, k, j)
                            alpha(t, l) = alpha(t, l) + jj(i, k, j)
                        else
                            beta(l) = beta(l) - gamma(t) *jj(i, k, j)
                        end if
                    else 
                        beta(t)=beta(t)-gamma(l)*jj(i, k ,j)
                    end if
                end do
            end if
    
    
    
            !Step 19:
            if (l<=n) then
            alpha(l,l)=alpha(l,l)+jj(i,k,k)
            beta(l)=beta(l)+ii(i,k)
            
        end do
    end do
    !Step 20:
    call ge(alpha,beta,gamma,n)

    !Step 21:
    do k=1,m
    write(*,*) gamma(k)
    end do

    !Step 22:
    do i = 1, mm
        do j = 1,3,1
            write(*,*) a(i, j)
            write(*,*) b(i, j)
            write(*,*) c(i, j)
        end do
    end do
    
    !Step 23:
    stop

end program finiteelement

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
            print*,'error in ge'
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
        print*,'error in ge'
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

subroutine gdi(chooser,J,a,b,c1,d1,ii,jj,kk,aa,bb,c,mm)
    implicit none
    real a,b,c1,d1,h1,h2,x,y,fff,JX,Q,k1,k2,J,aa(3,mm),bb(3,mm),c(3,mm)
    integer n,i,o,m,chooser,kk,ii,jj,mm

    real, dimension(1:5,1:5):: cc
    real, dimension(1:5,1:5):: rr

    rr(2,1)=0.5773502692
    rr(2,2)=-0.5773502692
    rr(3,1)=0.7745966692
    rr(3,2)=0.0000000000
    rr(3,3)=-0.7745966692
    rr(4,1)=0.8611363116
    rr(4,2)=0.3399810436
    rr(4,3)=-0.3399810436
    rr(4,4)=-0.8611363116
    rr(5,1)=0.9061798459
    rr(5,2)=0.5384693101
    rr(5,3)=0.0000000000
    rr(5,4)=-0.5384693101
    rr(5,5)=-0.9061798459


    cc(2,1)=1.0000000000
    cc(2,2)=1.0000000000
    cc(3,1)=0.5555555556
    cc(3,2)=0.8888888886
    cc(3,3)=0.5555555556
    cc(4,1)=0.3478548451
    cc(4,2)=0.6521451549
    cc(4,3)=0.6521451549
    cc(4,4)=0.3478548451
    cc(5,1)=0.2369268850
    cc(5,2)=0.4786286705
    cc(5,3)=0.5688888889
    cc(5,4)=0.4786286705
    cc(5,5)=0.2369268850


    n=5
    m=5


    if(n>5 .or. m>5) then
        write(*,*) ' the algorithm only works for n,m <= 5'
        stop
    endif


    h1=(b-a)/2
    h2=(b+a)/2
    J=0.0
    do i=1,m,1
        JX=0.0
        x=h1*rr(m,i)+h2
        k1=(d1-c1)/2
        k2=(d1+c1)/2
        do o=1,n,1
            y=k1*rr(n,o)+k2
            Q=fff(chooser,kk,jj,ii,x,y,aa,bb,c,mm)
            JX=JX+cc(n,o)*Q


        enddo
        J=J+cc(m,i)*k1*JX

    enddo
    J=h1*J


end subroutine

real function g1(x,y)
implicit none
real x,y
g1=0.0
endfunction

real function g2(x,y)
implicit none
real x,y
g2=(sqrt(2.0)/2.0)*(y-x)
endfunction

real function nnn(j,i,x,y,a,b,c,mm)
    implicit none
    integer i,j,mm
    real x,y,a(3,mm),b(3,mm),c(3,mm)


        nnn=a(j,i)+b(j,i)*x+c(j,i)*y

endfunction

real function g(x,y)
    implicit none
    real x,y
    g=4.0*x*y
end function

real function fff(cc,k,j,i,x,y,a,b,c,mm)
implicit none
integer i,j,cc,k,mm
real x,y,a(3,mm),b(3,mm),c(3,mm),p,q,r,nnn,f,g1,g2

if(cc==1) then
fff=p(x,y)
elseif(cc==2) then
fff=q(x,y)
elseif(cc==3) then
fff=r(x,y)*nnn(j,i,x,y,a,b,c,mm)*nnn(k,i,x,y,a,b,c,mm)
elseif(cc==4) then
fff=f(x,y)*nnn(j,i,x,y,a,b,c,mm)
elseif(cc==5) then
fff=g1(x,y)*nnn(j,i,x,y,a,b,c,mm)*nnn(k,i,x,y,a,b,c,mm)
elseif(cc==6) then
fff=g2(x,y)*nnn(j,i,x,y,a,b,c,mm)*nnn(k,i,x,y,a,b,c,mm)
endif
end function

real function p(x,y)
implicit none
real x,y
p=y**2
end function

real function q(x,y)
implicit none
real x,y
q=y**2
end function

real function r(x,y)
implicit none
real x,y
r=-y
end function
real function f(x,y)
implicit none
real x,y
f=-x
end function



  
