program cn
    implicit none
    integer, parameter :: m = 10, n = 50
    real :: ll, tt, alpha, w(m), x, t, lamda, k, h, f
    integer :: i, j
    real :: l(m-1), u(m-1), z(m-1)

    ll = 1.0
    tt = 0.5
    alpha = 1.0

    h = ll / m
    k = tt / n
    lamda = alpha**2 * k / (2.0 * h**2)

    w(m) = 0.0

    do i = 1, m-1
        w(i) = f(i * h)
    end do

    l(1) = 1.0 + lamda
    u(1) = -lamda / (2.0 * l(1))

    do i = 2, m-2
        l(i) = 1.0 + lamda + lamda * u(i-1) / 2.0
        u(i) = -lamda / (2.0 * l(i))
    end do

    l(m-1) = 1.0 + lamda + lamda * u(m-2) / 2.0

    do j = 1, n
        t = j * k
        z(1) = ((1.0 - lamda) * w(1) + lamda * w(2) / 2.0) / l(1)

        do i = 2, m-1
            z(i) = ((1.0 - lamda) * w(i) + lamda * (w(i+1) + w(i-1) + z(i-1)) / 2.0) / l(i)
        end do

        w(m-1) = z(m-1)

        do i = m-2, 1, -1
            w(i) = z(i) - u(i) * w(i+1)
        end do

        print *, 't =', t

        do i = 1, m-1
            x = i * h
            print *, 'x =', x, 'w =', w(i)
        end do
    end do

end program

real function f(x)
    implicit none
    real :: x

    f = sin(3.14156 * x)

end function
