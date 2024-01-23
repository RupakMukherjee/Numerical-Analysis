program NM
  implicit none
  
  integer, parameter :: n = 3
  integer :: iter_max, i, k, j0
  real :: tol
  real, dimension(:), allocatable :: x, y, f
  real, dimension(:,:), allocatable :: J
  
  allocate(x(n), y(n), f(n), J(n,n))
  
  x(1) = 0.1
  x(2) = 0.1
  x(3) = -0.1
  
  tol = 10e-4
  
  iter_max = 1000
  
  ! Step 1
  k = 1
  
  ! Step 2
  do while (k <= iter_max)
    
    ! Step 3
    f(1) = 3*x(1) - cos(x(2)*x(3)) - 0.5
    f(2) = x(1)**2 - 81*(x(2)+0.1)**2 + sin(x(3)) + 1.06
    f(3) = exp(-x(1)*x(2)) + 20*x(3) + (10*3.14 - 3)/3.0
    
    J = reshape([3.0, x(3)*sin(x(2)*x(3)), x(2)*sin(x(2)*x(3)), &
                 2*x(1), -162*(x(2)+0.1), cos(x(3)), &
                 -x(2)*exp(-x(1)*x(2)), -x(1)*exp(-x(1)*x(2)), 20.0], [n, n])
    
    ! Step 4
    call gaussianwithbs(J, y, f)
    
    ! Step 5
    x = x + y
    
    ! Step 6
    if (norm2(y) < tol) then
      print *, "Solution converged to: ", x
      stop
    end if
    
    ! Step 7
    k = k + 1
    
  end do
  
  ! Step 8 :
  print *, "Maximum number of iterations exceeded."
  stop
  
contains

subroutine gaussianwithbs(J, y, f)
  implicit none
  
  integer, parameter :: n = 3
  real, dimension(n, n), intent(in) :: J
  real, dimension(n), intent(inout) :: y
  real, dimension(n), intent(in) :: f
  
  real, dimension(n, n+1) :: a, O
  real, dimension(n) :: E, O2, x
  real :: sumf
  integer :: i, z, p0
  
  a(:, 1:n) = J
  a(:, n+1) = f
  
  do i = 1, n-1
    ! partial pivoting
    p0 = i
    do z = i+1, n
      if (abs(a(p0, i)) < abs(a(z, i))) p0 = z
    end do
    if (p0 /= i) then
      O(i,:) = a(p0, :)
      a(p0, :) = a(i, :)
      a(i, :) = O(i,:)
    end if
    
    do z = i+1, n
      a(z, i:n+1) = a(z, i:n+1) - a(z, i)/a(i, i) * a(i, i:n+1)
    end do
  end do
  
  ! back substitution
  end subroutine gaussianwithbs

end program NM