program waveeqnfinitediff
implicit none

integer,parameter::n=20,m=10
real l,t,alpha,w(0:m,0:n),h,k,lamda,f,g,x
integer i,j


print*,'endpoint l'
read(*,*) l

print*, 'max time T'
read(*,*) t

print*, "alpha"
read(*,*) alpha

h=l/m
k=t/n
lamda=k*alpha/h

    do j=1,n,1
        w(0,j)=0
        w(m,j)=0
    enddo
    w(0,0)=f(0.0)
    w(m,0)=f(l)

    do i=1,m-1
        w(i,0)=f(i*h)
        w(i,1)=(1-lamda**2)*f(i*h)+((lamda**2)/2.0)*(f((i+1)*h)+f((i-1)*h))+k*g(i*h)
    enddo

    do j=1,n-1
        do i=1,m-1
            w(i,j+1)=2*(1-lamda**2)*w(i,j)+(lamda**2)*(w(i+1,j)+w(i-1,j))-w(i,j-1)
        enddo
    enddo

    do i=0,m
        x=i*h

        do j=0,n
            t=j*k
            if(j==20) then
            print*, x,t,w(i,j)
            endif
        enddo
    enddo

    end program waveeqnfinitediff


real function f(x)
implicit none
real x
f=sin(3.14*x)
endfunction
real function g(x)
implicit none
real x
g=0.0
endfunction
