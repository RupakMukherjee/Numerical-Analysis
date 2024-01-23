program FalsePosition
    implicit none
    real,external :: f
    real :: p1,p0,Tol,q,q0,q1,p
    integer :: i,n

    Tol=1E-9
    p1=4
    p0=3
    n=1000

    i=2
    q0=f(p0)
    q1=f(p1)

    do while(i<=n)
        p=p1-q1*(p1-p0)/(q1-q0)
    if(abs(p-p1)<Tol) then
     write(*,*)"The Approximate solution is",p
        Stop
    end if

    i=i+1
    q=f(p)
    if(q*q1<0) then
        p0=p1
        q0=q1
    end if
    p1=p
    q1=q
    write(*,*)p
    end do

    write(*,*)"The method failed after n iterations, n=",n



end program FalsePosition

real function f(x)
real::x
f=cos(x)-x
return
end function