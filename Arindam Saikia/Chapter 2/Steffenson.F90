program Steffenson
    implicit none
    real,external :: f,g
    real :: p,p0,Tol,p2,p1
    integer :: i,n

    i=1
    Tol=1E-9
    p0=3
    n=1000

    do while(i<=n)
        p1=f(p0)
        p2=f(p1)
        p=p0-(p1-p0)**2/(p2-2*p1+p0)
    if(abs(p-p0)<Tol) then
     write(*,*)"The Approximate solution is",p
        Stop
    end if

    i=i+1
    p0=p
    write(*,*)p

    end do
    
end program Steffenson

real function f(x)
real,intent(in)::x
f=x**3+4*x**2-10
return
end function