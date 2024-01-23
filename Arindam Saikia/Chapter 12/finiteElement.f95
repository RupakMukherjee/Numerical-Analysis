program finite
  implicit none

  integer, parameter :: ne = 100   
  integer, parameter :: nn = 101   
  integer, parameter :: nq = 10    
  real(kind=8), parameter :: pi = 3.14

  ! Variables
  real(kind=8) :: a, b, h, x, y, x0, x1
  real(kind=8), dimension(nn) :: xval, yval, u
  real(kind=8), dimension(nq) :: q, w
  real(kind=8), dimension(nn, nn) :: aMatrix
  real(kind=8), dimension(nn) :: f, s

  ! Subroutine declarations
  abstract interface
    subroutine gaussQuad_interface(n, a, b, x, w)
      integer, intent(in) :: n
      real(kind=8), intent(in) :: a, b
      real(kind=8), dimension(n), intent(out) :: x, w
    end subroutine gaussQuad_interface
  end interface

  subroutine gaussQuad(n, a, b, x, w)
    integer, intent(in) :: n
    real(kind=8), intent(in) :: a, b
    real(kind=8), dimension(n), intent(out) :: x, w
  end subroutine gaussQuad

  ! Main program
  integer :: i, ie, iq, jq, k  ! Declaration added

  h = (2.0 * pi) / (ne - 1)

  ! Generate the mesh
  do i = 1, nn
    xval(i) = (i - 1) * h
    yval(i) = sin(xval(i))
  end do

  ! Initialize the stiffness matrix and the load vector
  aMatrix = 0.0
  f = 0.0

  ! Loop over the elements
  do ie = 1, ne-1
    x0 = xval(ie)
    x1 = xval(ie+1)

    ! Compute the element stiffness matrix
    call gaussQuad(nq, -1.0, 1.0, q, w)
    do iq = 1, nq
      x = x0 + ((x1 - x0) / 2.0) * (q(iq) + 1.0)
      y = sin(x)

      do jq = 1, nq
        a = (x1 - x0) / 2.0
        b = (x1 - x0) / 2.0 * q(jq)
        s(iq) = a * b
      end do

      do jq = 1, nq
        do k = 1, nn
          aMatrix(k, k) = aMatrix(k, k) + s(iq) * s(jq)
        end do

        f(iq) = f(iq) + s(iq) * y * s(jq)
      end do
    end do

    ! Assemble the element stiffness matrix and the load vector
    do i = 1, nn
      do k = 1, nn  ! Changed j to k
        aMatrix(i, k) = aMatrix(i, k) + aMatrix(i, k)
      end do

      f(i) = f(i) + f(i)
    end do
  end do

  ! Solve the system of equations
  call solveSystem(nn, aMatrix, f, u)

  ! Output the solution
  do i = 1, nn
    write(*,*) xval(i), u(i)
  end do

contains

  ! Subroutine to solve the system of equations using Gaussian elimination
  subroutine solveSystem(n, a, b, x)
    integer, intent(in) :: n
    real(kind=8), dimension(n, n), intent(in) :: a
    real(kind=8), dimension(n), intent(in) :: b
    real(kind=8), dimension(n), intent(out) :: x

    ! Local variables
    real(kind=8), dimension(n, n+1) :: augMatrix
    integer :: i, j, k
    real(kind=8) :: ratio

    ! Augment the coefficient matrix with the load vector
    augMatrix(:, 1:n) = a
    augMatrix(:, n+1) = b

    ! Gaussian elimination
    do k = 1, n-1
      do i = k+1, n
        ratio = augMatrix(i, k) / augMatrix(k, k)
        do j = k+1, n+1
          augMatrix(i, j) = augMatrix(i, j) - ratio * augMatrix(k, j)
        end do
      end do
    end do

    ! Back substitution
    x(n) = augMatrix(n, n+1) / augMatrix(n, n)
    do i = n-1, 1, -1
      x(i) = augMatrix(i, n+1)
      do j = i+1, n
        x(i) = x(i) - augMatrix(i, j) * x(j)
      end do
      x(i) = x(i) / augMatrix(i, i)
    end do
  end subroutine solveSystem

end program finite
