program linear_shooting
  implicit none
  real, external :: p, q, r
  real, dimension(:,:), allocatable :: w, u, v, k, kk
  real, dimension(:), allocatable :: y, W1, W2
  real :: a, b, h, alpha, beta, x
  integer :: i, N
  
  N = 10
  a = 1.0
  b = 2.0
  alpha = 1.0
  beta = 2.0
  
  h = (b - a) / N
  
  allocate(u(1:2, 0:N))
  allocate(v(1:2, 0:N))
  allocate(w(1:2, 0:N))
  allocate(k(1:4, 1:2))
  allocate(kk(1:4, 1:2))
  allocate(y(0:N))
  allocate(W1(1:N))
  allocate(W2(1:N))
  
  y(0) = alpha
  y(N) = beta
  
  u(1, 0) = alpha
  u(2, 0) = 0.0
  v(1, 0) = 0.0
  v(2, 0) = 1.0
  
  do i = 0, N-1
      x = a + i * h
      
      k(1, 1) = h * u(2, i)
      k(1, 2) = h * (p(x) * u(2, i) + q(x) * u(1, i) + r(x))
      
      k(2, 1) = h * (u(2, i) + 0.5 * k(1, 2))
      k(2, 2) = h * (p(x + h/2) * (u(2, i) + 0.5 * k(1, 2)) + q(x + h/2) * (u(1, i) + 0.5 * k(1, 1)) + r(x + h/2))
      
      k(3, 1) = h * (u(2, i) + 0.5 * k(2, 2))
      k(3, 2) = h * (p(x + h/2) * (u(2, i) + 0.5 * k(2, 2)) + q(x + h/2) * (u(1, i) + 0.5 * k(2, 1)) + r(x + h/2))
      
      k(4, 1) = h * (u(2, i) + k(3, 2))
      k(4, 2) = h * (p(x + h) * (u(2, i) + k(3, 2)) + q(x + h) * (u(1, i) + k(3, 1)) + r(x + h))
      
      u(1, i+1) = u(1, i) + (1.0/6.0) * (k(1, 1) + 2 * k(2, 1) + 2 * k(3, 1) + k(4, 1))
      u(2, i+1) = u(2, i) + (1.0/6.0) * (k(1, 2) + 2 * k(2, 2) + 2 * k(3, 2) + k(4, 2))
      
      kk(1, 1) = h * v(2, i)
      kk(1, 2) = h * (p(x) * v(2, i) + q(x) * v(1, i))
      
      kk(2, 1) = h * (v(2, i) + 0.5 * kk(1, 2))
      kk(2, 2) = h * (p(x + h/2) * (v(2, i) + 0.5 * kk(1, 2)) + q(x + h/2) * (v(1, i) + 0.5 * kk(1, 1)))
      
      kk(3, 1) = h * (v(2, i) + 0.5 * kk(2, 2))
      kk(3, 2) = h * (p(x + h/2) * (v(2, i) + 0.5 * kk(2, 2)) + q(x + h/2) * (v(1, i) + 0.5 * kk(2, 1)))
      
      kk(4, 1) = h * (v(2, i) + kk(3, 2))
      kk(4, 2) = h * (p(x + h) * (v(2, i) + kk(3, 2)) + q(x + h) * (v(1, i) + kk(3, 1)))
      
      v(1, i+1) = v(1, i) + (1.0/6.0) * (kk(1, 1) + 2 * kk(2, 1) + 2 * kk(3, 1) + kk(4, 1))
      v(2, i+1) = v(2, i) + (1.0/6.0) * (kk(1, 2) + 2 * kk(2, 2) + 2 * kk(3, 2) + kk(4, 2))
  end do
  
  w(1, 0) = alpha
  w(2, 0) = (beta - u(1, N)) / (v(1, N))
  
  print*, 'x:', a, 'y(x):', w(1, 0), 'y''(x):', w(2, 0)
  
  do i = 1, N
      W1(i) = u(1, i) + w(2, 0) * v(1, i)
      W2(i) = u(2, i) + w(2, 0) * v(2, i)
      x = a + i * h
      
      print*, 'x:', x, 'y(x):', W1(i), 'y''(x):', W2(i)
  end do
  
  stop
end program

real function p(g)
  implicit none
  real :: g
  p = -2.0 / g
  return
end function

real function q(g)
  implicit none
  real :: g
  q = 2.0 / (g**2)
  return
end function

real function r(g)
  implicit none
  real :: g
  r = (sin(log(g))) / (g**2)
  return
end function
