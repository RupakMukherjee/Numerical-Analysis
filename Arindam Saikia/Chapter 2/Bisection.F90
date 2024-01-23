program bisection
    implicit none
    integer :: i = 1
    real :: a, b, c, result

    print *, "Enter the value of a and b: "
    read *, a, b

    do
        c = (a + b) / 2.0
        result = FUNC(c)

        print *, "i=", i, "a=", a, "b=", b, "F(c)=", result

        if (FUNC(a) * result < 0) then
            b = c
        else 
            a = c
        end if

        i = i + 1

        if (abs(result) <= 0.001) exit
    end do

    print *, "Approximate root = ", c

contains

    real function FUNC(x)
        real, intent(in) :: x
        FUNC = x*x*x + 4*x*x - 10
    end function FUNC

end program bisection
