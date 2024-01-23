program fastfourier
    
    implicit none
    real, parameter::pi=3.142857
    integer,parameter::m=2,p=1
    complex y(0:2*m-1),zeta,zai(0:m+m),c(0:2*m-1),n
    integer l,i,j,k,mm,q,qq,kkk
    real a(0:m),b(1:m-1),k1,k2,kk(0:p)

    ! INPUTS -
    y(0)=(1,0)
    y(1)=(-1,0)
    y(2)=(1,0)
    y(3)=(-1,0)


    ! Step 1
    mm = m
    q = p
    zeta = cmplx(cos(pi/M), sin(pi/M))

    ! Step 2
    do j = 0, 2*m-1 ,1
        c(j) = y(j)
    end do

    ! Step 3
    do j = 1,mm,1
        zai(j) = zeta ** j
        zai(j+M) = -zai(j)
    end do
    
    ! Step 4
    K = 0
    zai(0) = 1

    ! Step 5
    do L = 1, p+1,1 
        ! Step 6
        do while (k < 2*m-1)
            ! Step 7
            do j = 1, mm,1
                ! Step 8
                kkk = k
                do i=p,0,-1
                    if(mod(kkk,2)==0) then
                        kk(p)=0
                    else
                        kk(p)=1
                    endif
                    kkk=kkk/2
                enddo
                k1=0
                k2=0

                do i=q,p,1
                    k1=k1+kk(i)*2**(i-q)
                    k2=kk(i)*2**(p-(i-q))
                end do

            ! Step 9
                n=c(k)+mm*zeta*k2
                c(k+m)=c(k)-n
                c(k)=c(k)+n

            ! Step 10
                K = K + 1
            end do
            
            ! Step 11
            k=k+mm
        end do

        ! Step 12
        k=0
        mm=mm/2
        q=q-1
    end do
    
    ! Step 13
    do while (K < 2*m-1)
        ! Step 14
        kkk=k
        do i=p,0,-1
            if(mod(kkk,2)==0) then
                kk(p)=0
            else
                kk(p)=1
            endif
            kkk=kkk/2

        enddo
        j=0
        do i=0,p,1
            j=j+kk(i)*2**(p-i)

        enddo
        
         ! Step 15
        if (j>k) then
            qq=c(k)
            c(k)=c(j)
            c(j)=qq
        end if
        
        ! Step 16
        k=k+1
    end do

    ! Step 17
    a(0)=c(0)/m
    a(m)=real(cmplx(cos(pi*m),-sin(pi*m))*cmplx(c(m))/cmplx(m,0))

    ! Step 18
    
    do j=1,m-1
        a(j)=real(cmplx(cos(pi*j),-sin(pi*j))*cmplx(c(j))/cmplx(m,0))
        b(j)=aimag(cmplx(cos(pi*j),-sin(pi*j))*cmplx(c(j))/cmplx(m,0))

    enddo

    ! Step 19
    print*,c
    print*,"---"
    print*,a
    print*,"---"
    print*,b
    stop


end program fastfourier