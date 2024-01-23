program Romberg
    implicit none

    integer :: i,j,n,k
    reaL :: a,b,h,sum,sumf
    real,dimension(:,:),allocatable  :: R
    real,external :: f


    a=0
    b=1
    n=4
    h=b-a
    

    allocate(R(1:n,1:n))

    R(1,1)=(h/2.0)*(f(a)+f(b))

    print*, "Value of R(1,1) is :", R(1,1)


    do i=2,n,1
       do k=1,2**(i-2)
        sumf=sum+f((a + (k-0.5)*h))
       end do
    end do

    
    do i=2,n
        
        R(2,1)=(R(1,1)+h*sumf)/2.0

        do j=2,i
            write(*,"(a,i0,a)") "Value for R(2,", j, "):"
            R(2,j)=R(2,j-1)+(R(2,j-1)-R(1,j-1))/(4**(j-1)-1)
             write(*,*) R(2,j)
        end do
             h=h/2.0

        do j=1,i
            R(1,j)=R(2,j)
        end do  
        
          
    end do
    
    end program Romberg

real function f(x)
real,intent(in) :: X
f=sin(x)
!f=x**2+4
return
end function