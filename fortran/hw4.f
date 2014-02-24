c  Adam Martinez
c  2 FEB 2010
c  3450:730 ANSPDE
c  Homework 4

c  Use the tridiagonal solver to solve Ax=b.

      implicit none
      double precision L(4),M(4),R(4),b(4),x(4)
	integer n,nmax

	nmax = 4

	do n = 1,nmax
	L(n) = -1.0d0
	M(n) = 2.0d0
	R(n) = -1.0d0
	end do
	L(1) = 0
	R(nmax) = 0

	b(1) = 1.0d0
	b(2) = 0.0d0
	b(3) = 0.0d0
	b(4) = 1.0d0

	call tridiag(L,M,R,b,x,nmax)

	do n = 1,nmax
	write(6,*) '[',x(n),']'
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
             write(6,*) a
             write(6,*) b
             write(6,*) c
             write(6,*) r
             write(6,*) n
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
