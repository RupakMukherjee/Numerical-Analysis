program iterativeRefine
    implicit none
    integer i, j, k, NN, kk, nnn
    integer, parameter :: n = 4
    real a(n, n), b(n), xo(n), tol, x(n), dum, dums, m(n, n - 1), y(n), r(n), cond, yo(n), t, xx(n)
    NN = 10
    nnn = 1
    t = 1.0 / 10.0

    a(1, 1) = 10
    a(1, 2) = -1
    a(1, 3) = 2
    a(1, 4) = 0
    a(2, 1) = -1
    a(2, 2) = 11
    a(2, 3) = -1
    a(2, 4) = 3
    a(3, 1) = 2
    a(3, 2) = -1
    a(3, 3) = 10
    a(3, 4) = -1
    a(4, 1) = 0
    a(4, 2) = 3
    a(4, 3) = -1
    a(4, 4) = 8

    b(1) = 6
    b(2) = 25
    b(3) = -11
    b(4) = 15

    xo(1) = 0
    xo(2) = 0
    xo(3) = 0
    xo(4) = 0

    tol = 0.001

    k = 1
    do while (k <= NN)
        do i = 1, n, 1
            dum = 0
            dums = 0
            do j = 1, i - 1, 1
                dums = dums + a(i, j) * x(j)
            enddo
            do j = i + 1, n, 1
                dum = dum + a(i, j) * xo(j)
            enddo
            x(i) = (1.0 / a(i, i)) * (-dums - dum + b(i))
        enddo

        if (all(abs(x - xo) < tol)) then
            exit
        endif

        k = k + 1

        xo = x
    enddo

    k = 1
    do while (k <= nnn)
        do i = 1, n, 1
            dum = 0
            do j = 1, n, 1
                dum = dum + a(i, j) * x(j)
            enddo
            r(i) = b(i) - dum
        enddo

        yo(1) = 0
        yo(2) = 0
        yo(3) = 0
        yo(4) = 0

    kk=1
    do while(kk<=NN)
        do i=1,n,1
            dum=0
            dums=0

            do j=1,i-1,1
                dums=dums+a(i,j)*y(j)
            enddo
            do j=i+1,n,1
                dum=dum+a(i,j)*yo(j)
            enddo
            y(i)=(1/a(i,i))*(-dums-dum+r(i))
        enddo

            if(abs(y(1)-yo(1))<tol .and. abs(y(2)-yo(2))<tol .and. abs(y(3)-yo(3))<tol .and. abs(y(4)-yo(4))<tol) then

                exit
            endif

            kk=k+1

            do i=1,n,1
                yo(i)=y(i)
            enddo
    enddo
    do i=1,n,1
        xx(i)=x(i)+y(i)
    enddo
    if(k==1) then
            cond=(y(1)/xx(1))*10**t
    endif
    if(abs(x(1)-xx(1))<tol) then
        write(*,*) 'xx'
        write(*,*) xx
        !write(*,*) 'cond:', cond
        stop
    endif

    k=k+1
    do i=1,n,1
        x(i)=xx(i)
    enddo
    enddo


!Step 10 :
print*,"Maximum number of iterations exceeded"
print*,COND
stop

    
end program iterativeRefine

