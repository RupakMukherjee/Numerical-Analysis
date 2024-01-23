program CompSimp
implicit none

integer :: n,i,j
real :: a,b,int,h,XI0,XI1,XI2,X,XI
real,external :: f

n=8.0
a=0
b=4

h=(b-a)/n

XI0=f(a)+f(b)
XI1=0
XI2=0


do i=1,n-1
    X=a+i*h

    if (MOD(i,2)/=0) then
        XI1=XI1+f(X)
    
    else 
        
        XI2=XI2+f(X)
    end if


    
end do

XI=h*(XI0+2*XI2+4*XI1)/3.0

print*,"The approximation for the integral is :"
write(*,*) XI
    

end program

real function f(x)
real,intent(in) :: X
f=exp(x)
return
end function