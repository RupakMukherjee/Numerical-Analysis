program symmetric_power_method
  implicit none
  integer,parameter::n=3
  real a(n,n),x(n,1),tol,nn,mu(1,1),er,y(n,1),dum(n),xt(1,n)
  integer i,j,k,p

  tol=10e-8
  nn=100
  a(1,1)=4
  a(1,2)=-1
  a(1,3)=1
  a(2,1)=-1
  a(2,2)=3
  a(2,3)=-2
  a(3,1)=1
  a(3,2)=-2
  a(3,3)=3
  x(1,1)=1
  x(2,1)=0
  x(3,1)=0
    
  k=1
  do i=1,n
      x(i,1)=x(i,1)/norm2(x(1:n,1))
  enddo


  do while(k<=nn)
      call mm(a,x,y,n,n,1)
      !do i=1,n
          !xt(1,i)=x(i,1)
      !enddo

      !call mm(xt,y,mu,1,n,1)
      call mm(x,y,mu,1,n,1)

      if(norm2(y)==0) then
          write(*,*) 'eigenvector',x
          write(*,*) 'a has the eigenvalue 0, select a new vector x and restart'
          stop
      endif
      do i=1,n,1
          dum(i)=x(i,1)-(y(i,1)/norm2(y))
      enddo
      er=norm2(dum)
      do i=1,n,1
          x(i,1)=y(i,1)/norm2(y)
      enddo
      if (er<tol) then
          write(*,*) 'mu=',mu,'x=',x
          stop
      endif
      k=k+1
  enddo
  write(*,*)'max no. of iterations exceeded'

end program symmetric_power_method

subroutine mm(a,b,c,n,m,p)
  implicit none
  integer n,m,p
  real a(n,m),b(m,p),c(n,p),dum
  integer i,j,k

  do i=1,n,1
      do j=1,p,1
          dum=0
          do k=1,m,1
              dum=dum+a(i,k)*b(k,j)
          enddo
          c(i,j)=dum
      enddo
  enddo

endsubroutine
  