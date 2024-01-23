program Secant
implicit none
real,external :: f
real :: p1,p0,Tol,q0,q1,p
integer :: i=1,n

Tol=1E-9
p1=0
p0=2
n=1000

i=2
q0=f(p0)
q1=f(p1)

do while(i<=n)
    p=p1-q1*(p1-p0)/(q1-q0)
    if(abs(p-p1)<Tol) then
     write(*,*)"The approximate solution is:",p
        Stop
    end if
i=i+1
p0=p1
q0=q1
p1=p
q1=f(p)
write(*,*)p
end do

write(*,*)"The method failed after n iterations, n=",n

end program Secant

real function f(x)
real::x
f=cos(x)-x
return
end function