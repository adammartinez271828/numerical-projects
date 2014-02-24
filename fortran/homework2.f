c  Adam Martinez
c  19 JAN 2010
c  3450:730 ANSPDE
c  Homework 2

c  Numerically solve the heat equation, Ut = DUxx
c  U(	0,t) = f(t)
c  U(xmax,t) = g(t)
c  U(	x,0) = h(x)

      implicit none
      double precision U(1001,100001),x(1001),t(100001)
      double precision D,tmin,tmax,xmin,xmax,dt,dx,f,g,h,r,Uxx
      integer nmax,n,imax,i

c  D - thermal conductivity
c  xmin,xmax - endpoints of spacial domain
c  tmin,tmax - endpoints of time domain
c  dt,dx - finite differences of domains
c  f,g - left and right	boundaries
c  h - initial condition
c  r - discretization parameter
c  N - number of points	in time grid
c  I - number of points	in spacial grid

      D = 3.0d-2

      tmin = 0.0d0
      tmax = 1.0d0
      xmin = 0.0d0
      xmax = 1.0d0

      dt = 1.0d-5
      dx = 1.0d-3

      nmax = (tmax-tmin)/dt + 1.5d0
      imax = (xmax-xmin)/dx + 1.5d0

      r = D*dt/(dx**2)

      write(6,*) 'dt = ',dt
      write(6,*) 'dx = ',dx
      write(6,*) 'r  = ',r
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
	U(i,1) = h(x(i))
      end do

c  Interior point calculation

      do n = 1,nmax-1
        if (n==(n-n%10000)) then
          write(6,*) n
        end if
	U(1,n+1) = f(t(n+1))
	do i = 2,imax-1
	  Uxx = U(i-1,n) - 2*U(i,n) + U(i+1,n)
	  U(i,n+1) = U(i,n) + r*Uxx
	end do
	u(imax,n+1) = g(t(n+1))
      end do

c  Output

      open(9,file = 'hw2output.txt')
      do i = 1,imax
	write(9,*) x(i),U(i,nmax)
      end do
      close(9)

      do i = 1,10
	write(6,*) t(i),U(i,nmax)
      end do

      stop
      end

      function f(t)
      implicit none
      double precision f,t
      f = 3.0d0
      return
      end

      function g(t)
      implicit none
      double precision g,t
      g = 1.0d0
      return
      end

      function h(x)
      implicit none
      double precision h,x
      h = 0.0d0
      return
      end
