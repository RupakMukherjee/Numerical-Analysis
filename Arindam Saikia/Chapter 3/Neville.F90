program Neville
    implicit none
    
    real :: f
    real :: xo = 8.4
    integer :: i, j, n
    real, dimension(:), allocatable :: x
    real, dimension(:,:), allocatable :: Q

    print *, "Give the value of n:"
    read *, n

    allocate(x(0:n))
    allocate(Q(0:n,0:n))

    do i = 0, n
        print *, "Give the value of x(", i, "):"
        read *, x(i)
        Q(i, 0) = f(x(i))
    end do

    do j = 1, n
        do i = j, n
            Q(i, j) = ((xo - x(i-j)) * Q(i, j-1) - (xo - x(i)) * Q(i-1, j-1)) / (x(i) - x(i-j))
        end do
    end do

    print *, "The Q matrix is:"
    do i = 0, n
        do j = 0, i
            write(*, *) Q(i, j)
        end do
    end do

    deallocate(x)
    deallocate(Q)

end program Neville


real function f(x)
    real, intent(in) :: x
    !f = x*x*x - x
    !f = exp(x)
    f = log(x)
end function f
