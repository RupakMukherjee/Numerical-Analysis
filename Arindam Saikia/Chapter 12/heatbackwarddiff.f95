program heat_equation
    implicit none
    integer, parameter :: N = 100, m = 10
    real :: ll, T, alpha, lambda, tt, k, h
    real :: w(m), l(m), u(m), z(m), x
    integer :: i, j
    real, external :: f

    ll = 1.0
    T = 1.0
    alpha = 2.0

    h = ll / m
    k = T / N
    lambda = alpha**2 * k / h**2

    do i = 1, m-1
        w(i) = f(i * h)
    end do

    l(1) = 1 + 2 * lambda
    u(1) = -lambda / l(1)

    do i = 2, m-2
        l(i) = 1 + 2 * lambda + lambda * u(i-1)
        u(i) = -lambda / l(i)
    end do

    l(m-1) = 1 + 2 * lambda + lambda * u(m-2)

    do j = 1, N
        tt = j * k
        z(1) = w(1) / l(1)

        do i = 2, m-1
            z(i) = (w(i) + lambda * z(i-1)) / l(i)
        end do

        w(m-1) = z(m-1)

        do i = m-2, 1, -1
            w(i) = z(i) - u(i) * w(i+1)
        end do

        do i = 1, m-1
            x = i * h
            print *, 't =', tt, 'x =', x, 'w =', w(i)
        end do
    end do
end program

real function f(x)
    real :: x
    f = sin(3.141 * x)
end function
