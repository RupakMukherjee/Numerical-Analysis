program ccs
    implicit none
    
    integer :: i,j
    integer,parameter :: n=4
    real :: x(0:n),a(0:n),b(0:n),c(0:n),d(0:n),h(0:n),alpha(0:n),l(0:n), mu(0:n), z(0:n)
    real :: aj(0:n), bj(0:n), cj(0:n), dj(0:n)
    real :: FPO , FPN
    
    write(*,*) 'Enter the values of x:'
    do i = 0, n
        READ(*,*) x(i)
    end do
    
    write(*,*) 'Enter the values of f(x):'
    do i = 0, n
        READ(*,*) a(i)
    end do

    write(*,*) 'Enter the value of f`(x0):'
        READ(*,*) FPO

    write(*,*) 'Enter the value of f`(xn):'
        READ(*,*) FPN
    
    
    do i=0,n-1
        h(i)=x(i+1)-x(i)
    end do
    
    alpha(0)=(3*(a(1)-a(0)))/(h(0)-3*FPO)
    alpha(n)=(3*FPN-3*(a(n)-a(n-1)))/(h(n-1))

    do i=1,n-1
        alpha(i)=3.0/h(i) * (a(i+1) - a(i)) - 3.0/h(i-1) * (a(i) - a(i-1))
    end do
    
    l(0) = 2*h(0)
    mu(0) = 0.5
    z(0) = alpha(0)/l(0)
    
    do i = 1,n-1,1
        l(i) = 2.0 * (x(i+1) - x(i-1)) - h(i-1) * mu(i-1)
        mu(i) = h(i) / l(i)
        z(i) = (alpha(i) - h(i-1)*z(i-1)) / l(i)
    end do
    
    l(n) = h(n-1)*(2-mu(n-1))
    z(n) = (alpha(n)-h(n-1)*z(n-1))/l(n)
    c(n) = z(n)
    
    do j = n-1,0
        c(j) = z(j) - mu(j) * cj(j+1)
        b(j) = (a(j+1) - a(j))/h(j) - h(j)*(cj(j+1) + 2.0*cj(j))/3.0
        d(j) = (cj(j+1) - cj(j))/(3.0*h(j))
    end do
    
    WRITE(*,*) 'j', 'aj', 'bj', 'cj', 'dj'
    do j=0,n-1
    write(*,*)j,a(j),b(j),c(j),d(j)
    end do
    
    end program
    