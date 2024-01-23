program gaussianwithbs
    implicit none

    integer,parameter :: n=2 , p=1
    integer :: i,j
    real,dimension(1:n,1:n) :: m
    real,dimension(1:n, 1:n+1) :: a
    real,dimension(1:n) :: E,O,O2
    real :: sumf
    real,dimension(1:n) :: x


    a(1,1)=1
    a(1,2)=2
    a(2,1)=2
    a(2,2)=1
    a(1,3)=5
    a(2,3)=4

    E(1)=a(1,1)*x(1)+a(1,2)*x(2)-a(1,3)
    E(2)=a(2,1)*x(1)+a(2,2)*x(2)-a(2,3)



    do i=1,n-1,1
        !if (i<=p .and. p<=n .and. a(p,i) /= 0) then
        
        !else 
        !    print*,"No unique solution exists"
        !end if
        
        if (p/=i) then

            O(p)=E(P)
            E(p)=E(i)
            E(i)=O(p)
               
              do j=i+1,n
                m(i,j)=a(j,i)/a(i,i)


                E(j)=(E(j)-m(j,i)*E(i))
              
            
              end do


        end if


    end do

    if (a(n,n)==0) then
        print*,"No unique solution exists:"
        stop
    end if

    x(n)=a(n,n+1)/a(n,n)

    do i=n-1,1,-1
        sumf=0
        do j=i+1,n
            sumf=sumf+a(i,j)*x(j)
           end do
        x(i)=(a(i,n+1)-sumf)/a(i,i)
    end do 

    do i=1,n
    write(*,*)x(i)
    end do
    
end program gaussianwithbs