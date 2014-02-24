c  Adam Martinez
c  23 FEB 2010
c  3450:730 ANSPDE
c  Computer Project 1

c  Numerically solve the diffusion equation, Ut = D(x)Uxx + E(x)Ux + S(x,t)
c  U (xmin,t) = f(t)
c  Ux(xmax,t) = g(t)
c  U (   x,0) = h(x)

      implicit none
      double precision U(3001),x(3001),t(4001)
	double precision L(3001),M(3001),R(3001),bvec(3001)
      double precision rt,bt,tmin,tmax,xmin,xmax,dt,dx
	double precision f,g,h,D,E,S,Uxx,gm,gc,gp
      integer nmax,n,imax,i

c  xmin,xmax - endpoints of spacial domain
c  tmin,tmax - endpoints of time domain
c  dt,dx - finite differences of domains
c  f,g - left and right	boundaries
c  h - initial condition
c  N - number of points	in time grid
c  I - number of points	in spacial grid

      tmin = 0.0d0
      tmax = 2.0d0
      xmin = 1.0d0
      xmax = 4.0d0

	dx = 1.0d-2
	dt = 2.0d0*dx
	open(9,file = 'cp1a.txt')

c	dx = 1.0d-2
c	dt = dx
c	open(9,file = 'cp1b.txt')

c	dx = 1.0d-2
c	dt = dx/2.0d0
c	open(9,file = 'cp1c.txt')

c	dx = 1.0d-3
c	dt = 2.0d0*dx
c	open(9,file = 'cp1d.txt')

c	dx = 1.0d-3
c	dt = dx
c	open(9,file = 'cp1e.txt')

c	dx = 1.0d-3
c	dt = dx/2.0d0
c	open(9,file = 'cp1f.txt')

      nmax = (tmax-tmin)/dt + 1.5d0
      imax = (xmax-xmin)/dx + 1.5d0

      write(6,*) 'dt = ',dt
      write(6,*) 'dx = ',dx
      write(6,*) 'nmax = ',nmax
      write(6,*) 'imax = ',imax

c  Generate t and x meshes

      do n = 1,nmax
	t(n) = tmin + dfloat(n-1)*dt
      end do

      do i = 1,imax
	x(i) = xmin + dfloat(i-1)*dx
      end do

c  Apply initial condition

	do i = 1,imax
	U(i) = h(x(i))
      end do

c  Interior point calculation

	do n = 1,nmax-1

c  Construct Ax=b system for step n->n+1
c  First row
	  L(1) = 0.0d0
	  M(1) = 1.0d0
	  R(1) = 0.0d0
	  bvec(1) = f(t(n+1))

c  Interior rows
c  gm, gc, gp are dependent on x(i), so must be inside the interior loop.
	  do i = 2,imax-1
        rt = D(x(i))*dt/(dx**2)
c	  write(6,*) rt
	  bt = E(x(i))*dt/(2.0d0*dx)
	  gm = (rt-bt)/2.0d0
	  gc = rt
	  gp = (rt+bt)/2.0d0
c	  write(6,*) gm, gc, gp

	  L(i) = -gm
	  M(i) = 1.d0+gc
	  R(i) = -gp
	  bvec(i) = gm*U(i-1) + (1-gc)*U(i)
	  bvec(i) = bvec(i) + gp*U(i+1) + dt*S(x(i),t(n)+dt/2.0d0)
	  end do

c  Last row
        rt = D(x(imax))*dt/(dx**2)
	  bt = E(x(imax))*dt/(2.0d0*dx)
	  gm = (rt-bt)/2.0d0
	  gc = rt
	  gp = (rt+bt)/2.0d0

	  i = imax
	  L(i) = -gm-gp
	  M(i) = 1+gc
	  R(i) = 0.0d0
	  bvec(i) = 2.0d0*dx*gp*(g(t(n+1))+g(t(n)))
	  bvec(i) = bvec(i) + (gm+gp)*U(i-1) + (1-gc)*U(i) 
	  bvec(i) = bvec(i) + dt*S(x(i),t(n)+dt/2.0d0)
c  Matrix has been constructed
c  Call tridiag and load data into U(i)

	call tridiag(L,M,R,bvec,U,imax)

	end do

c  Output

	do i = 1,imax
	write(9,*) x(i),U(i)
	end do
	close(9)

	do i = 1,10
	write(6,*) x(i),U(i)
	end do

      stop
      end

      function f(t)
      implicit none
      double precision f,t
      f = 1.0d0 + dsin(t)
      return
      end

      function g(t)
      implicit none
      double precision g,t
      g = -5.0d0 - exp(-t)
      return
      end

      function h(x)
      implicit none
      double precision h,x
      h = dexp(-50.0d0*(x-2.0d0)**2)
      return
      end

	function D(x)
	implicit none
	double precision D,x
	D = 4.0d0 + dsin(7.0d0*x)
	return
	end

	function E(x)
	implicit none
	double precision E,x
	E = 1.0d-3
	return
	end

	function S(x,t)
	implicit none
	double precision S,x,t,pi
	pi = 4.0d0*datan(1.0d0)
	S = 8000.0d0 * dsin( 4.0d0*pi*(x-1.0d0) ) * dexp(-t)
	return
	end

      SUBROUTINE tridiag(a,b,c,r,u,n)
      implicit none
      INTEGER n
      double precision a(n),b(n),c(n),r(n),u(n)
      INTEGER j
      double precision bet,gam(n)

c     a is AL, b is AM, c is AR
c     r is rhs, u is output, n is dimension

      bet=b(1)
      u(1)=r(1)/bet
      do 11 j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        u(j)=(r(j)-a(j)*u(j-1))/bet
11    continue
      do 12 j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
12    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +%6V+j)D2.