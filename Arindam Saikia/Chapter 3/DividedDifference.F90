program ndd
    implicit none
    
    integer :: i,j
    integer,parameter :: n=4
    real :: f(0:n+1),x(0:n+1),Q(0:n+1,0:n+1)

    do i = 0, n
        write(*,*) 'Enter the value of x(', i, '):'
        read(*,*) x(i)
        
        write(*,*) 'Enter the value of f(', i, '):'
        read(*,*) f(i)
        
        Q(i, 0) = f(i)
    end do

    do j=1,n
        do i=j,n
            Q(i,j)=(Q(i, j-1) - Q(i-1, j-1)) / (x(i) - x(i-j))
        end do
    end do
    
    write(*,*) 'The divided-difference coefficients are:'
    do i = 0, n
        write(*,*) Q(i, i)
    end do
    
end program ndd
