program power_method
    implicit none
    integer :: n, p, k, max_iter, i, j
    real :: tol, xp, mu, yp, err
    real, dimension(:,:), allocatable :: A
    real, dimension(:), allocatable :: x, y

    ! INPUTS
    n = 2
    max_iter = 1000
    tol = 10E-4
    allocate(A(n,n))
    !A = (/-2,-3,6,7/)
    A(1,1)=-2
    A(1,2)=-3
    A(2,1)=6
    A(2,2)=7
    
    allocate(x(n))
    !x = (/1,1/)
    x(1)=1
    x(2)=1

    ! Step 1
    k = 1

    do while (k <= max_iter)
        ! Step 2
        p = 1
        xp = abs(x(p))
        do i = 2, n
            if (abs(x(i)) > xp) then
                p = i
                xp = abs(x(p))
            end if
        end do

        ! Step 3
        x = x / xp

        ! Step 4
        y = matmul(A, x)

        ! Step 5
        yp = y(p)

        ! Step 6
        mu = yp

        ! Step 7
        p = 1
        yp = abs(y(p))
        do i = 2, n
            if (abs(y(i)) > yp) then
                p = i
                yp = abs(y(p))
            end if
        end do

        ! Step 8
        if (yp == 0.0) then
            write(*,*) "A has the eigenvalue 0, select a new vector x and restart"
            write(*,*) "Eigenvector: ", x
            stop
        end if

        ! Step 9
        ERR=(norm2(x-(y/yp)))
        x=y/yp

        ! Step 10
        if (ERR<Tol) then
            print*,mu
            print*,x
        end if

        ! Step 11
        k=k+1
    end do

    ! Step 12
    print*,"The maximum number of iteratins exceeded"
    stop

end program power_method