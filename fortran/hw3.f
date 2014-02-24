c  Adam Martinez
c  04 MAR 2010
c  3450:730 ANSPDE
c  Computer Project 2

c  Numerically solve the heat equation, Ut = alpha*Uxx + sin(U) + (Ux)^2
c  U(	  0,t) = f(t)
c  U(xmax,t) = g(t)
c  U(	  x,0) = h(x)

      implicit none
      double precision U(1001),x(1001),t(10001),
     .tmin,tmax,xmin,xmax,dt,dx,
     .f,g,h,rt,Uxx,Ux,
     .alpha,err,tol
      integer nmax,n,imax,i,mmax,m

	character*1 cflag
	character*13 flname
	character*3 chan3
	integer np

c  D - thermal conductivity
c  xmin,xmax - endpoints of spacial domain
c  tmin,tmax - endpoints of time domain
c  dt,dx - finite differences of domains
c  f,g - left and right	boundaries
c  h - initial condition
c  r - discretization parameter
c  N - number of points	in time grid
c  I - number of points	in spacial grid

	write(6,*) 'Enter ID flag (1 character)'
	read(5,77) cflag
77	format(a1)

	np = 10000

      alpha = 5.0d-1

      tmin = 0.d0
      tmax = 1.d0
      xmin = 0.d0
      xmax = 1.d0

      dt = 1.d-6
      dx = 1.d-3

      nmax = (tmax-tmin)/dt + 1.5d0
      imax = (xmax-xmin)/dx + 1.5d0

      rt = alpha*dt/(dx**2)

      write(6,*) 'dt = ',dt
      write(6,*) 'dx = ',dx
      write(6,*) 'rt  = ',rt
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

      t = 0
	call ooutfn('smp2',chan3(t),cflag,1)
	do i = 1,imax
	  U(i,0) = h(x(i))
	  write(1,*) x(i),U(i,0)
      end do

c  Interior point calculation
c  U(i,0) is U_i^(k)
c  U(i,1) is U_i^(k+1)

      do n = 1,nmax-1
	  U(1,1) = f(t(n+1))

	  do i = 2,imax-1
	    Uxx = U(i-1,0) - 2*U(i,0) + U(i+1,0)
	    U(i,1) = U(i,0) + r*Uxx
	  end do

	  U(imax,1) = g(t(n+1))

	  if (np*((n+1)/np) .eq. n+1) then
	    call ooutfn('smp2',chan3((n+1)/np),cflag,1)
	    write(6,*) n+1,'    ',chan3((n+1)/np)

	    do i = 1,imax
	      write(1,*) x(i),U(i,1)
	    end do

	    close(1)
	  end if

	  do i = 1,imax
	    U(i,0) = U(i,1)
	  end do

	end do

c  Output

c	open(9,file = 'hw2output.txt')
c	do i = 1,imax
c	  write(9,*) x(i),U(i,1)
c	end do
c	close(9)

c	do i = 1,10
c	  write(6,*) x(i),U(i,1)
c	end do

      stop
      end

      function f(t)
      implicit none
      double precision f,t
        f = 1.d0
      return
      end

      function g(t)
      implicit none
      double precision g,t
        g = -3.d0
      return
      end

      function h(x)
      implicit none
      double precision h,x
        h = 3.d0*exp(-500.d0*(x-1.d0/3.d0)**2)+1.d0-4.d0*x
      return
      end

	subroutine ooutfn(fpr,fid,cflag,lun)
c  usage
c  type *,' Enter file id'
c  accept (2a), id
c  call ooutfn(id,'i',1) creates output file i'id'.dat
c  taken from Nazanin Imani

	character*3 fid
	character*1 cflag
	character*4 fpr
	integer lun
	character*20 flname

	flname = fpr // cflag // fid
	open(unit=lun,file=flname,err=150)
	return
150	write(6,*),'Error opening file.'
	return
	end

	function chan3(M)
c  Stephen Cardarelli, Oct 2003
c  Provides time step files 000 to 999
	implicit none
	character*3 chan3
	integer i,j,k,M,countout

	countout = M
	i = countout/100
	countout = countout - (i*100)
	j = countout/10
	k = countout - (j*10)

	chan3 = char(i+48)//char(j+48)//char(k+48)
	return
	end