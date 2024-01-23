program fixedpoint
    implicit none
    integer :: i=1,n
    real :: p,p0,Tol

    print*,"Enter the initial approximation value:"
    read*,p0

    print*, "Give the tolerance amount:"
    read*, Tol

    print*,"Enter the maximum number of iterations:"
    read*, n

    do while(i<=n)
        p=g(p0)
        if(abs(p-p0)<Tol) then
            print*,"Approximate Solution is :",p
            stop
        end if

    i=i+1
    p0=p

    write(*,*)p
    end do
    
!-----------------
contains

real function  g(x)
 real,intent(in):: x
 g=x-(x**3+4*x**2-10)/(3*x**2+8*x)
end function

end program fixedpoint