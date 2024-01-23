program Householders
    implicit none
  
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer :: i, j, k, l, n
    real(dp) :: alpha, rsq, prod
    real(dp), dimension(:,:), allocatable :: A
    real(dp), dimension(:), allocatable :: v, u, z
  
    ! Read the size of the matrix
    write(*,*) 'Enter the size of the matrix:'
    read(*,*) n
  
    ! Allocate memory for the matrix and vectors
    allocate(A(n,n), v(n), u(n), z(n))
  
    ! Read the matrix A from standard input
    do j = 1, n
       write(*,*) 'Enter the elements of row', j, ' of the matrix A:'
       read(*,*) (A(j,i), i=1,n)
    end do
  
    ! Perform Householder reduction
    do k = 1, n-2
       ! Compute the Householder vector v
       alpha = dot_product(A(k+1:n,k), A(k+1:n,k))
       v(k+1:n) = A(k+1:n,k)
       v(k+1) = v(k+1) + sign(sqrt(alpha), v(k+1))
       rsq = dot_product(v(k+1:n), v(k+1:n))
       v(k+1:n) = v(k+1:n) / (sqrt(2.0_dp * rsq))
  
       ! Compute the matrix H and apply it to A
       u(k+1:n) = dot_product(v(k+1:n), A(k+1:n,k))

       z(k+1:n) = 2.0_dp * u(k+1:n) / dot_product(v(k+1:n), v(k+1:n))
       do l = k+1, n
          A(l,k:n) = A(l,k:n) - v(l)*z(k:n) - z(l)*v(k:n)
       end do
    end do
  
    ! Print the resulting matrix A
    write(*,*) 'The resulting matrix A is:'
    do j = 1, n
       write(*,*) (A(j,i), i=1,n)
    end do
  
    ! Deallocate memory
    deallocate(A, v, u, z)
  
  end program Householders
  