program nonlinear_shooting
  implicit none
  real, external :: p, q, r, dp, dq, dr
  real, dimension(:,:), allocatable :: u, v, k, kk
  real, dimension(:), allocatable :: W1, W2
  real :: a, b, h, alpha, beta, tol
  real :: y0, y1, yp0, yp1, y
  real :: x
  integer :: i, N, max_iter, iter
  logical :: converged

  N = 10
  tol = 1.0E-6
  max_iter = 100

  a = 1.0
  b = 2.0
  alpha = 1.0
  beta = 2.0

  h = (b - a) / N
  allocate(u(0:1, 0:N))
  allocate(v(0:1, 0:N))
  allocate(k(0:1, 0:1))
  allocate(kk(0:1, 0:1))
  allocate(W1(1:N))
  allocate(W2(1:N))

  ! Initial guesses for y0 and y1
  y0 = alpha
  y1 = alpha + h

  iter = 0
  converged = .false.

  do while (.not. converged)
      iter = iter + 1

      yp0 = (y1 - y0) / h
      yp1 = yp0 + h * ((p(a) * y0 + q(a) * yp0 + r(a)) - alpha) / dp(a)

      u(0, 0) = y0
      u(1, 0) = yp0
      v(0, 0) = 0.0
      v(1, 0) = 1.0

      do i = 0, N-1
          x = a + i * h
          k(0, 0) = h * u(1, i)
          k(0, 1) = h * (p(x) * u(1, i) + q(x) * u(0, i) + r(x))
          k(1, 0) = h * (u(1, i) + 0.5 * k(0, 1))
          k(1, 1) = h * (p(x + h/2) * (u(1, i) + 0.5 * k(0, 1)) + q(x + h/2) * (u(0, i) + 0.5 * k(0, 0)) + r(x + h/2))
          u(0, i+1) = u(0, i) + (1.0/6.0) * (k(0, 0) + 2 * k(1, 0))
          u(1, i+1) = u(1, i) + (1.0/6.0) * (k(0, 1) + 2 * k(1, 1))

          kk(0, 0) = h * v(1, i)
          kk(0, 1) = h * (p(x) * v(1, i) + q(x) * v(0, i))
          kk(1, 0) = h * (v(1, i) + 0.5 * kk(0, 1))
          kk(1, 1) = h * (p(x + h/2) * (v(1, i) + 0.5 * kk(0, 1)) + q(x + h/2) * (v(0, i) + 0.5 * kk(0, 0)))
          v(0, i+1) = v(0, i) + (1.0/6.0) * (kk(0, 0) + 2 * kk(1, 0))
          v(1, i+1) = v(1, i) + (1.0/6.0) * (kk(0, 1) + 2 * kk(1, 1))
      end do

      y = u(0, N) + (beta - u(0, N)) * v(0, N) / v(0, N)

      if (abs(y - y1) < tol .and. abs(yp1 - yp0) < tol) then
          converged = .true.
      elseif (iter > max_iter) then
          converged = .true.
          print*, "Maximum iterations reached. Solution did not converge."
      else
          y0 = y1
          y1 = y
      end if
  end do

  print*, "x:", a, "y(x):", y0, "y'(x):", yp0

  do i = 1, N
      W1(i) = u(0, i) + y0 * v(0, i)
      W2(i) = u(1, i) + yp0 * v(1, i)
      x = a + i * h
      print*, "x:", x, "y(x):", W1(i), "y'(x):", W2(i)
  end do

  stop
end program

real function p(x)
  implicit none
  real :: x
  p = -2.0 / x
  return
end function

real function q(x)
  implicit none
  real :: x
  q = 2.0 / (x**2)
  return
end function

real function r(x)
  implicit none
  real :: x
  r = (sin(log(x))) / (x**2)
  return
end function

real function dp(x)
  implicit none
  real :: x
  dp = 2.0 / (x**2)
  return
end function

real function dq(x)
  implicit none
  real :: x
  dq = -4.0 / (x**3)
  return
end function

real function dr(x)
  implicit none
  real :: x
  dr = (-2.0 * cos(log(x)) - sin(log(x))) / (x**2)
  return
end function
