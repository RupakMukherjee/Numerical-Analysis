program Newtons
implicit none
real,external :: f,g
real :: p,p0,Tol
integer :: i=1,n

print*,"Enter the initial approximation value:"
read*,p0

print*,"Enter the tolerance amount:"
read*, Tol

print*,"Enter the Maximum number of iterations:"
read*, n

do while(i<=n)
    p=p0-(f(p0)/g(p0))

    if(abs(p-p0)<Tol) then
        write(*,*)"The approximate solution is :",p
        stop
    end if
i=i+1
p0=p
write(*,*)p
end do

end program Newtons

real function f(x)
real :: x
f=cos(x)-x
end function

real function g(x)
real :: x
g=-sin(x)-1
end function