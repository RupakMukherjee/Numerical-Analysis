program jacobiiterative
    implicit none
    
    integer,parameter :: n=4
    integer :: k, i, j, iter_max
    real :: tol
    real, dimension(:,:), allocatable :: A
    real, dimension(:), allocatable :: b, x, x_old
    
    !print *, "Enter the number of equations and unknowns (n):"
    !read *, n
    
    allocate(A(n,n))
    allocate(b(n))
    allocate(x(n))
    allocate(x_old(n))
    
    !print *, "Enter the entries of the coefficient matrix A:"
    !do i = 1, n
    !    do j = 1, n
    !        read *, A(i,j)
    !    end do
    !end do
    
    A(1,1)=10
    A(1,2)=-1
    A(1,3)=2
    A(1,4)=0
    A(2,1)=-1
    A(2,2)=11
    A(2,3)=-1
    A(2,4)=3
    A(3,1)=2
    A(3,2)=-1
    A(3,3)=10
    A(3,4)=-1
    A(4,1)=8
    A(4,2)=3
    A(4,3)=-1
    A(4,4)=8

    !print *, "Enter the entries of the right-hand side vector b:"
    !do i = 1, n
    !    read *, b(i)
    !end do

    b(1)=6
    b(2)=25
    b(3)=-11
    b(4)=15
   
    !print *, "Enter the initial approximation X0:"
    !do i = 1, n
    !    read *, x_old(i)
    !end do
    
    x_old(1)=0
    x_old(2)=0
    x_old(3)=0
    x_old(4)=0

    !print *, "Enter the tolerance (TOL):"
    !read *, tol
    
    tol=10e-5

    !print *, "Enter the maximum number of iterations (iter_max):"
    !read *, iter_max
    
    iter_max=1000

    k = 1
    x(1)=0
    x(2)=0
    x(3)=0
    x(4)=0
    do while (k <= iter_max)
        do i = 1, n
            do j = 1, n
                if (j /= i) then
                    x(i) = x(i) + A(i,j)*x_old(j)
                end if
            end do
            x(i) = (-x(i)+b(i))/A(i,i)
        end do
    
        
        
        if (abs(x(i) - x_old(i)) < tol) then
            print *, "Solution found:"
            !do i = 1, n
            
                print *, x(1)
                print *, x(2)
                print *, x(3)
                print *, x(4)

                !print *, x(i)
            stop
            !end do
        end if
        
        
        k = k + 1
        x_old = x

    end do
    
    print *, "Maximum number of iterations exceeded."
    stop
    
    end program jacobiiterative
    