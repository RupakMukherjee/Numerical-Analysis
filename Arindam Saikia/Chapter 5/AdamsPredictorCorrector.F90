program adamspredictorcorrector
    implicit none

    integer :: i,j,N
    real :: a,b,alpha,h,k1,k2,k3,k4,wf,tf
    real, dimension(:), allocatable :: t,w
    real,external :: f
    
    print*,"Provide the end points a and b:"
    read(*,*)a,b
    
    print*,"Provide the value of N:"
    read(*,*)N

    print*,"Provide the initial condition, alpha:"
    read(*,*)alpha

h=(b-a)/N

allocate(t(0:N))
allocate(w(0:N))

t(0)=a
w(0)=alpha

print*,"The Output is:"
do i=1,3,1
    k1=h*f(t(i-1),w(i-1))
    k2=h*f(t(i-1)+h/2.0,w(i-1)+k1/2.0)
    k3=h*f(t(i-1)+h/2.0,w(i-1)+k2/2.0)
    k4=h*f(t(i-1)+h,w(i-1)+k3)

        w(i)=w(i-1)+(k1+2*k2+2*k3+k4)/6.0
        t(i)=a+i*h
end do

do i=4,N
    tf=a+i*h
    wf=w(3)+h*(55*f(t(3),w(3))-59*f(t(2),w(2))+37*f(t(1),w(1))-9*f(t(0),w(0)))/24.0
    wf=w(3)+h*(9*f(tf,wf)+19*f(t(3),w(3))-5*f(t(2),w(2))+f(t(1),w(1)))/24.0


    
    write(*,*)tf,wf
    do j=0,2,1
        t(j)=t(j+1)
        w(j)=w(j+1)
    end do
t(3)=tf
w(3)=wf
end do

end program adamspredictorcorrector

real function f(t,y)
real,intent(in) :: t,y
f=y-t*t+1.0
end function