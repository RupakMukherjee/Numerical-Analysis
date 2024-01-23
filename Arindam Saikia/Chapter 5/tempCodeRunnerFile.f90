program rungekuttadesystems
    implicit none

    integer :: i,j,N,m
    real :: a,b,t,h
    real, dimension(:), allocatable :: alpha, w, k1, k2, k3, k4

    ! Allocate arrays
    m = 2
    allocate(alpha(1:m))
    alpha = [1.0, 2.0]
    allocate(w(1:m))
    allocate(k1(1:m))
    allocate(k2(1:m))
    allocate(k3(1:m))
    allocate(k4(1:m))

    ! Set initial conditions
    a = 0.0
    b = 2.0
    N = 10
    h = (b-a)/N
    t = a
    w = alpha

    ! Output initial conditions
    do j = 1, m
        write(*,*) t, w(j)
    end do

    ! Time loop
    do i = 1, N
        ! Compute k values
        do j = 1, m
            k1(j) = h * f(j, t, w)
            k2(j) = h * f(j, t + h / 2.0, w + k1 / 2.0)
            k3(j) = h * f(j, t + h / 2.0, w + k2 / 2.0)
            k4(j) = h * f(j, t + h, w + k3)
        end do

        ! Update w values
        do j = 1, m
            w(j) = w(j) + (k1(j) + 2.0 * k2(j) + 2.0 * k3(j) + k4(j)) / 6.0
        end do

        ! Increment time
        t = a + i * h

        ! Output w values
        do j = 1, m
            write(*,*) t, w(j)
        end do
    end do

contains

    ! Function to compute f values
    function f(j, t, w) result(fj)
        implicit none
        integer, intent(in) :: j
        real, intent(in) :: t, w(:)
        real :: fj

        if (j == 1) then
            fj = 4.0*t + 3.0*w(1) + 6.0
        else
            fj = 2.4*t + 1.6*w(2) + 3.6
        end if

    end function f

end program rungekuttadesystems
