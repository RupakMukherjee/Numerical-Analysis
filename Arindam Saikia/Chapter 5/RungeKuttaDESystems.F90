program rungekuttadesystems
    implicit none

    integer, parameter :: m=2
    integer :: i,j,N
    real :: a,b,t,h,k(4,m),w(m),alpha(1:m),f
    
    alpha = [0.0, 0.0]


    !a = 0.0
    !b = 1.0
    !N = 5
    !h = (b-a)/N
    h=0.1
    t = a

    do j = 1, m
        w(j)=alpha(j)
    end do
    write(*,*) t, w

    do i=1,N,1
        do j=1,m,1
            k(1,j)=h*f(t,w(1),w(2),j)
        enddo
        do j=1,m,1
            k(2,j)=h*f(t+h/2,w(1)+k(1,1)/2,w(2)+k(1,2)/2,j)
        enddo
        do j=1,m,1
            k(3,j)=h*f(t+h/2,w(1)+k(2,1)/2,w(2)+k(2,2)/2,j)
        enddo
        do j=1,m,1
            k(4,j)=h*f(t+h,w(1)+k(3,1),w(2)+k(3,2),j)
        enddo
        do j=1,m,1
            w(j)=w(j)+(k(1,j)+2*k(2,j)+2*k(3,j)+k(4,j))/6
        enddo

  
    t = a + i * h

    write(*,*) t, w

end do
    
end program rungekuttadesystems

real function f(t,w1,w2,j)
real t,w
integer j

if(j==1) then
f=-4*w1+3*w2+6
elseif(j==2) then
f=0.6*(-4*w1+3*w2+6)-0.2*w2
elseif(j==3) then
f=0
endif
end function