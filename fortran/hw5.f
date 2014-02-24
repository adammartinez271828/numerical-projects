c  Adam Martinez
c  2 FEB 2010
c  3450:730 ANSPDE
c  Homework 4

c  Use the tridiagonal solver to solve Ax=b.

      implicit none
      double precision L(1001),M(1001),R(1001),b(1001),x(1001)
	integer n,nmax

	nmax = 1001

	do n = 2,nmax-1
	  L(n) =  1.0d0
	  M(n) = -2.0d0
	  R(n) =  1.0d0
	  b(n) = dfloat(n-1)
	end do
	L(1) = 0.0d0
	L(nmax) = 0.0d0
	M(1) = 1.0d0
	M(nmax) = 1.0d0
	R(1) = 0.0d0
	R(nmax) = 0.0d0
	b(1) = 0.0d0
	b(nmax) = 1000.0d0

	do n = 1,nmax
	  write(2,*) L(n),M(n),R(n),b(n)
	end do

	call tridiag(L,M,R,b,x,nmax)

	do n = 1,nmax
	write(6,*) '[',x(n),']'
	end do

	open(1,FILE='output.txt')
	
	do n = 1,nmax
	  write(1,*) x(n)
	end do

	stop
	end

      SUBROUTINE tridiag(a,b,c,r,u,n)
      implicit none
      INTEGER n
      double precision a(n),b(n),c(n),r(n),u(n)
      INTEGER j
      double precision bet,gam(n)

c     a is AL, b is AM, c is AR
c     r is rhs, u is output, n is dimension
      if(b(1).eq.0.)pause 'tridiag: rewrite equations'
      bet=b(1)
      u(1)=r(1)/bet
      do 11 j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if(bet.eq.0.)pause 'tridag failed'
        u(j)=(r(j)-a(j)*u(j-1))/bet
11    continue
      do 12 j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
12    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +%6V+j)D2.