program SORmethod
    implicit none
  
    integer,parameter :: n=3
    integer :: i, j, k, iter, iter_max
    real :: a(n,n), b(n), xo(n), x(n), omega, tol, error,sum,sum1
  
  
    A(1,1)=4.0
    A(1,2)=3.0
    A(1,3)=0.0

    A(2,1)=3.0
    A(2,2)=4.0
    A(2,3)=-1.0

    A(3,1)=0.0
    A(3,2)=-1.0
    A(3,3)=4.0
    
    b(1)=24.0
    b(2)=30.0
    b(3)=-24.0
    
    xo(1)=1.0
    xo(2)=1.0
    xo(3)=1.0

    tol=10E-5
    iter_max=20

    k = 1
    omega=1.25
    do while (k <= iter_max)

      do i = 1, n
        sum = 0
        do j = 1, i-1
            sum = sum + A(i,j)*x(j)
        end do
      
        sum1=0
        do j = i+1, n
            sum1 = sum1 + A(i,j)*xo(j)
        end do
        x(i) = (1.0 - omega) * xo(i) + omega*(b(i) - sum-sum1)/A(i,i)
      end do
  
      !if (x(i)-xo(i) < tol) then
      if (all(abs(x-xo) < tol)) then 
      print *, "Solution:"
        do i = 1, n
          print *, x(i)
        end do
        stop
      end if
  
      k = k + 1
      xo=x
  
      
    end do

    print *, "Maximum number of iterations exceeded"
    stop
  
  end program SORmethod

  !Desired Output : 3.0000498  ,  4.0002586   ,   -5.0003486
  