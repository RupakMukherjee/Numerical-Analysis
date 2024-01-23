program bezier
    implicit none

    integer :: i, t, n
    real, dimension(:,:), allocatable :: x, xp, xm, y, ym, yp
    real, dimension(:,:), allocatable :: a, b
    real :: u, v

    n = 5

    allocate(x(0:n, 0))
    allocate(xp(0:n-1, 0))
    allocate(xm(1:n, 0))

    allocate(y(0:n, 0))
    allocate(yp(0:n-1, 0))
    allocate(ym(1:n, 0))

    allocate(a(1:n, 0:n-1))
    allocate(b(1:n, 0:n-1))

    x(0, 0) = -1
    x(1, 0) = 0
    x(2, 0) = 1
    x(3, 0) = 0
    x(4, 0) = 1

    y(0, 0) = 0
    y(1, 0) = 1
    y(2, 0) = 0.5
    y(3, 0) = 0
    y(4, 0) = -1

    do i = 1, n
        a(i, 0) = x(i - 1, 0)
        b(i, 0) = y(i - 1, 0)
        a(i, 1) = 3 * (xp(i - 1, 0) - x(i - 1, 0))
        b(i, 1) = 3 * (yp(i - 1, 0) - y(i - 1, 0))
        a(i, 2) = 3 * (x(i, 0) + xm(i, 0) - 2 * xp(i - 1, 0))
        b(i, 2) = 3 * (y(i, 0) + xm(i + 1, 0) - 2 * yp(i - 1, 0))
        a(i, 3) = x(i, 0) - x(i - 1, 0) + 3 * xp(i - 1, 0) - 3 * xm(i, 0)
        b(i, 3) = y(i, 0) - y(i - 1, 0) + 3 * yp(i - 1, 0) - 3 * ym(i + 1, 0)
    end do

    t = 0.5

    do i = 1, n
        x(i - 1, t) = a(i, 0) + a(i, 1) * t + a(i, 2) * t**2 + a(i, 3) * t**3
        y(i - 1, t) = b(i, 0) + b(i, 1) * t + b(i, 2) * t**2 + b(i, 3) * t**3
        write(*, *) x(i - 1, t), y(i - 1, t)
    end do

    u = (((64 * t - 352 / 3.0) * t + 60) * t - 14 / 3.0) * t - 1
    v=((((-64/3.0)*t+48)*t-116/3.0)*t+11)*t

    print*, "coordinates:"
    write(*,*) u, v
    
    end program