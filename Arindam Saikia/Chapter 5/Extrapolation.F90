program extrapolation
    implicit none
    real :: a, b, alpha, tol, hmax, hmin, h,TO,WO
    real :: t, w, hk, yk, v, w1, w2, w3, q(8,8)
    integer :: i, j, k, nk(8), NFLAG, FLAG
    real,dimension(1:8) :: x,y
    real,external :: f
    
    nk = (/2, 4, 6, 8, 12, 16, 24, 32/)

    a=0
    b=2
    alpha=0.5
    tol=10E-4
    hmax=0.25
    hmin=0.01

    TO=a
    WO=alpha
    h=hmax
    FLAG=1

    do i=1,7,1
        do j=1,i,1
            Q(i,j)=(NK(i+1)/NK(j))**2

        end do
    end do

    do while (FLAG==1)
        k=1
        NFLAG=0

        do while (k<=8 .and. NFLAG==0)
            HK=h/NK(k)
            T=TO
            W2=WO
            W3=W2+HK*f(T,W2)
            T=TO+HK
            
            do j=1,NK(k)-1,1
            W1=W2
            W2=W3
            W3=W1+2.0*HK*f(T,W2)
            T=TO+(j+1)*HK
            end do

            y(k)=(W3+W2+HK*f(T,W3))/2.0

                if (K>=2) then
                j=K
                v=y(1)

                do while (j>=2)
                
                    y(j-1)=y(j)+(y(j)-y(j-1))/(Q(k-1,j-1)-1)
                    j=j-1
                    
                    end do
                
                        if (abs(y(1)-v)<=tol) then
                    NFLAG=1
                        end if
                end if  

                  k=k+1
        end do

        k=k-1
        
        if (NFLAG==0) then
            h=h/2.0
            if (h<hmin) then
                print*,"hmin exceeded"
                FLAG=0
            end if
        
        else 
            WO=y(1)
            TO=TO+h
            Write(*,*)TO,WO,h

            if (TO>=b) then
                FLAG=0

            else if (TO+h>b) then
                h=b-TO

            else if (k<=3 .and. h<0.5*hmax)then
                h=2.0*h
            end if
        end if
        
    end do





end program extrapolation

real function f(t,y)
real,intent(in) :: t,y
f=y-t*t+1.0
end function