program adamsvariablestepsize
    implicit none
    
    integer :: i,j,N
    real :: a,b,alpha,h,k1,k2,k3,k4,wf,tf,FLAG,LAST,sigma,tt,q,NFLAG,hmax,hmin,Tol,WC,WP
    real, dimension(0:3):: t,w
    real,external :: f

    a=0
    b=2
    hmax=0.25
    hmin=0.01
    alpha=0.5
    Tol=10E-4
    t(0)= a
    w(0)= alpha
    h=hmax
    FLAG=1
    LAST=0

    write(*,*) t(0), w(0)

    call rk4(h, w, t, 0, 3)

    NFLAG=1
    i=4
    tt=t(3)+h

    do while(FLAG==1)
        
        WP=w(i-1)+(h/24)*(55*f(t(i-1),w(i-1))-59*f(t(i-2),w(i-2))+37*f(t(i-3),w(i-3))-9*f(t(i-4),w(i-4)))
        WC=w(i-1)+(h/24)*(9*f(tt,WP)+19*f(t(i-1),w(i-1))-5*f(t(i-2),w(i-2))+f(t(i-3),w(i-3)))
        sigma=19*abs(WC-WP)/(270*h)

        if(sigma<=Tol) then
        w(i)=WC
        t(i)=tt
            if(NFLAG==1) then
               do j=i-3,i,1
               write(*,*) t(j),w(j)
            enddo
            else
               write(*,*) t(i),w(j)
            endif

            if(LAST==1) then
                FLAG=0
            else
                i=i+1
                NFLAG=0
                if(sigma<=0.1*Tol .or. (t(i-1)+h)>b) then
                    q=(Tol/(2*sigma))**(1/4)
                    if(q>4) then
                        h=4*h
                    else
                        h=q*h
                    endif

                    if(h>hmax) then
                        h=hmax
                    endif
                    if((t(i-1)+4*h)>b) then
                        h=(b-t(i-1))/4
                        LAST=1
                    endif

                    call rk4(h, w, t, i-1, i+2)
                    NFLAG=1
                    i=i+3
                endif
            endif
        endif

        q=(Tol/(2*sigma))**(1/4)

        if(q<0.1) then
            h=0.1*h
        else
            h=q*h
        endif

        if(h<hmin) then
            FLAG=0
            write(*,*)'hmin exceeded'
        else
            if(NFLAG==1) then
                i=i-3
                call rk4(h,w,t,i-1,i+2)
                i=i+3
                NFLAG=1
            endif
        endif
        tt=t(i-1)+h
    enddo
    write(*,*)'Process ended'

end program adamsvariablestepsize

real function f(t,w)
real, intent(in):: t, w
f=w-t**2+1
end function

subroutine rk4(h, v, x, m, n)

    implicit none
    real h,k1,k2,k3,k4,f
    integer j,m,n
    real, dimension(m:n):: v
    real, dimension(m:n):: x
    intent(in) h,m,n
    intent(inout) v,x

    do j=1,3,1

    k1=h*f(x(j-1),v(j-1))
    k2=h*f(x(j-1)+h/2,v(j-1)+k1/2)
    k3=h*f(x(j-1)+h/2,v(j-1)+k2/2)
    k4=h*f(x(j-1)+h,v(j-1)+k3)

    v(j)=v(j-1)+(k1+2*k2+2*k3+k4)/6
    x(j)=x(0)+j*h

    enddo



end subroutine rk4