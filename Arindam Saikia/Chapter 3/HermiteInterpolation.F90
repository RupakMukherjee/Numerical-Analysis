program Hermite
implicit none
    
!real:: f,f1
integer :: i, j, n,k
real, dimension(:), allocatable :: x,z,f,f1
real, dimension(:,:), allocatable :: Q
real :: H,l
    

!print *, "Give the value of n:"
!read *, n

n=2

allocate(x(0:n))
allocate(z(0:n))
allocate(f(0:n))
allocate(f1(0:n))
allocate(Q(0:2*n+1,0:2*n+1))

x(0)=1.3
x(1)=1.6
x(2)=1.9

f(0)=0.6200860
f(1)=0.4554022
f(2)=0.2818186

f1(0)=-0.5220232
f1(1)=-0.5698959
f1(2)=-0.5811571

!do i = 1, n+1
!    print *, "Give the value of x(",i-1,"):"
!    read *, x(i-1)
!end do

    z(0)=x(0)
    z(1)=x(0)
    Q(0,0)=f(0)
    Q(1,0)=f(0)
    Q(1,1)=f1(0)

do i = 1,n
    z(2*i)=x(i)
    z(2*i+1)=x(i)
    Q(2*i,0)=f(i)
   Q(2*i+1,0)=f(i)
   Q(2*i+1,1)=f1(i)
    
    Q(2*i,1)=(Q(2*i,0)-Q(2*i-1,0))/(z(2*i)-z(2*i-1))
    

end do

    do i=2,2*n+1
        do j=2,i

            Q(i,j)=(Q(i,j-1)-Q(i-1,j-1))/(z(i)-z(i-j))

        end do
    end do


    !write(*,*)Q(0,0)
    !write(*,*)Q(1,1)
    !write(*,*)Q(2,2)
    !write(*,*)Q(3,3)
    !write(*,*)Q(4,4)
    !write(*,*)Q(5,5)

!FIX :
    
! for odd : upto given sequence

! for even : upto (x-x(n-1))^2 term only

H=Q(0,0)+Q(1,1)*(1.5-1.3)+Q(2,2)*(1.5-1.3)*(1.5-1.3)&
+Q(3,3)*(1.5-1.3)*(1.5-1.3)*(1.5-1.6)+Q(4,4)*(1.5-1.3)&
*(1.5-1.3)*(1.5-1.3)+Q(5,5)*(1.5-1.3)*(1.5-1.3)*(1.5-1.3)*(1.5-1.9)

write(*,*)H

!do i=1,2*n+1
!write(*,*)Q(i,i)
!end do

    
end program
    
    !real function  f(x)
    !real,intent(in) :: x
    !f=exp(x)
    !return
    !end function
   ! 
   ! real function f1(x)
   ! real,intent(in) :: x
   ! f1=exp(x)
   ! return
   ! end function