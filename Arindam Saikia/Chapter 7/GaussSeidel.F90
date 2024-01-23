program gaussseidel
    implicit none

    integer, parameter :: n = 4
    integer :: i, j, k, iter, iter_max
    real :: tol, error, sum
    real, dimension(:,:), allocatable :: A
    real, dimension(:), allocatable :: b, x, xo

    allocate(A(n, n))
    allocate(b(n))
    allocate(x(n))
    allocate(xo(n))

    A(1, 1) = 10
    A(1, 2) = -1
    A(1, 3) = 2
    A(1, 4) = 0
    A(2, 1) = -1
    A(2, 2) = 11
    A(2, 3) = -1
    A(2, 4) = 3
    A(3, 1) = 2
    A(3, 2) = -1
    A(3, 3) = 10
    A(3, 4) = -1
    A(4, 1) = 0
    A(4, 2) = 3
    A(4, 3) = -1
    A(4, 4) = 8

    b(1) = 6
    b(2) = 25
    b(3) = -11
    b(4) = 15

    xo(1) = 0
    xo(2) = 0
    xo(3) = 0
    xo(4) = 0

    tol = 10E-5
    iter_max = 1000

    k = 1
    do while (k <= iter_max)
        do i = 1, n
            sum = 0
            do j = 1, i-1
                sum = sum + A(i, j)*x(j)
            end do
            do j = i+1, n
                sum = sum + A(i, j)*xo(j)
            end do
            x(i) = (b(i) - sum) / A(i, i)
        end do

        error = maxval(abs(x - xo))
        if (error < tol) then
            print *, "Converged in ", k, " iterations"
            print *, "Solution:"
            print *, x
            stop
        end if

        k = k + 1
        xo = x
    end do

    print *, "Maximum number of iterations exceeded"
    print *, "Solution:"
    print *, x

    deallocate(A)
    deallocate(b)
    deallocate(x)
    deallocate(xo)
end program gaussseidel
