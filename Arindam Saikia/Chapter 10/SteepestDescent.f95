program steepestdescent
implicit none
  
    ! Parameters
    integer, parameter :: n = 2 
    integer, parameter :: max_iterations = 100 
    real, parameter :: tolerance = 1e-6 
  
    ! Variables
    integer :: k 
    real :: objfun,alpha,g0,alpha0, h1, h2, h3,alpha1, alpha2, alpha3 ,z(n), z0 ,g1, g2, g3,x(n) 
    
      
  
    ! Initial approximation
    x = (/ 1.0, 1.0 /)
  
    !Step 1
    k = 1

    !Step 2
    do while (k <= max_iterations)
  
      !Step 3
      g1 = objfun(x)
      call gradfun(x, z)
      z0 = norm2(z)
  
      !Step 4
      if (z0 == 0.0) then
        write(*,*) "Zero gradient"
        write(*,*) "x =", x
        write(*,*) "f(x) =", g1
        stop
      end if
  
      !Step 5
      z = z / z0
      alpha1 = 0.0
      alpha3 = 1.0
      g3 = objfun(x - alpha3 * z)
  
      !Step 6
      do while (g3 >= g1)
        alpha3 = alpha3 / 2.0
        g3 = objfun(x - alpha3 * z)
        if (alpha3 < tolerance / 2.0) then
          write(*,*) "No likely improvement"
          write(*,*) "x =", x
          write(*,*) "f(x) =", g1
          stop
        end if
      end do
  
      !Step 9
      alpha2 = alpha3 / 2.0
      g2 = objfun(x - alpha2 * z)

      !Step 10
      h1 = (g2 - g1) / alpha2
      h2 = (g3 - g2) / (alpha3 - alpha2)
      h3 = (h2 - h1) / alpha3
      alpha0 = 0.5 * (alpha2 - h1 / h3)
      g0 = objfun(x - alpha0 * z)
  
      ! Find best step size
      if (g0 < g1) then
        alpha = alpha0
      else
        alpha = alpha3
      end if
  
      !Step 13
      x = x - alpha * z
      if (abs(g0 - g1) < tolerance) then
        write(*,*) "Success: solution found"
        write(*,*) "x =", x
        write(*,*) "f(x) =", g0
        stop
      end if
  
      !Step 15
      k = k + 1
  
    end do
  
    !Step 16
    write(*,*) "Failure: maximum iterations exceeded"
    write(*,*) "x =", x
    write(*,*) "f(x) =", g1
  
  contains
  subroutine gradfun(x, grad)
  implicit none
  !integer,parameter, intent(in) :: n=2
  real*8, dimension(n), intent(in) :: x
  real*8, dimension(n), intent(out) :: grad
  ! calculate the gradient of the function
  ! here is an example implementation
  grad(1) = 2.0d0 * x(1)
  grad(2) = 4.0d0 * x(2)**3
end subroutine gradfun

  
end program steepestdescent

function objfun(x) result(f)
    implicit none
    real, dimension(:), intent(in) :: x
    real :: f
    
    f = x(1)**2 + 2.0d0 * x(2)**2
    
  end function objfun