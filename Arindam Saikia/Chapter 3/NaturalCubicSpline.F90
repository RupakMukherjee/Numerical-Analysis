program ncs
    implicit none

    integer :: i,j
    integer, parameter :: n = 4
    real :: x(0:n), a(0:n), b(0:n), c(0:n), d(0:n), h(0:n), alpha(0:n), l(0:n), mu(0:n), z(0:n)
    real :: aj(0:n), bj(0:n), cj(0:n), dj(0:n)

    write(*,*) 'Enter the values of x:'
    do i = 0, n
        READ(*,*) x(i)
    end do

    write(*,*) 'Enter the values of f(x):'
    do i = 0, n
        READ(*,*) a(i)
    end do

    do i = 0, n-1
        h(i) = x(i+1) - x(i)
    end do

    do i = 1, n-1
        alpha(i) = 3.0/h(i) * (a(i+1) - a(i)) - 3.0/h(i-1) * (a(i) - a(i-1))
    end do

    l(0) = 1.0
    mu(0) = 0.0
    z(0) = 0.0

    do i = 1, n-1
        l(i) = 2.0 * (x(i+1) - x(i-1)) - h(i-1) * mu(i-1)
        mu(i) = h(i) / l(i)
        z(i) = (alpha(i) - h(i-1)*z(i-1)) / l(i)
    end do

    l(n) = 1.0
    z(n) = 0.0
    cj(n) = 0.0

    do j = n-1, 0, -1
        cj(j) = z(j) - mu(j) * cj(j+1)
        bj(j) = (a(j+1) - a(j))/h(j) - h(j)*(cj(j+1) + 2.0*cj(j))/3.0
        dj(j) = (cj(j+1) - cj(j))/(3.0*h(j))
        aj(j) = a(j)
    end do

    WRITE(*,*) 'j', 'aj', 'bj', 'cj', 'dj'
    do j = 0, n-1
        write(*,*) j, aj(j), bj(j), cj(j), dj(j)
    end do

end program ncs
