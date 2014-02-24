c  Adam Martinez
c  04 MAR 2010
c  3450:730 ANSPDE
c  Computer Project 2

c  Numerically solve the heat equation, Ut = alpha*Uxx + sin(U) + (Ux)^2
c  U(	  0,t) = f(t)
c  U(xmax,t) = g(t)
c  U(	  x,0) = h(x)
c  0 < x < 1
c  0 < t < inf

      implicit none
      double precision U(1001),v(1001),x(1001),t(100001),
     .b(1001),L(1001),M(1001),R(1001),output(1001),
     .tmin,tmax,xmin,xmax,dt,dx,
     .f,g,h,rt,bt,Uxx,Ux,
     .alpha,err,tol
      integer nmax,n,imax,i,mmax,mval

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

	np = 1

      alpha = 5.0d-1

      dt = 1.d-6
      dx = 1.d-3

      tmin = 0.d0
      tmax = 2*dt
c	tmax = 1.d-1
      xmin = 0.d0
      xmax = 1.d0

      nmax = (tmax-tmin)/dt + 1.5d0
      imax = (xmax-xmin)/dx + 1.5d0
	mmax = 50

c      rt = alpha*dt/(dx**2)
	rt = 5.d-1
	bt = 1.d0*dt/(2.d0*dx)

      write(6,*) 'dt   = ',dt
      write(6,*) 'dx   = ',dx
      write(6,*) 'rt   = ',rt
	write(6,*) 'bt   = ', bt
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

      n = 0
	call ooutfn('smp2',chan3(n),cflag,1)
	do i = 1,imax
	  U(i) = h(x(i))
	  write(1,*) x(i),U(i)
      end do
	tol = 1.d-6
	
c  Construct the matrix
c  First row
	  i = 1
	  L(i) = 0.d0
	  M(i) = 1.d0
	  R(i) = 0.d0

c  Interior rows
	  do i = 2,imax-1
	    L(i) = rt/2.d0
c  Unstable version
	    M(i) = -rt
c  Stable version
c	    M(i) = -(1.d0+rt)
	    R(i) = rt/2.d0
	  end do

c  Last row
	  i = imax
	  L(i) = 0.d0
	  M(i) = 1.d0
	  R(i) = 0.d0

c  Primary loop

	write(6,*) 'Beginning loop...'
	open(unit=2,file="ErrorOutput.txt")
	write(2,*) "Error Output"
      do n = 1,nmax-1
	write(2,*) 'Iteration ', n
	  do i = 1,imax
	    v(i) = u(i)
	  end do

	  err = 2.d0*tol
	  mval = 0

	  do while((err.gt.tol).and.(mval.le.mmax))

	    i = 1
	    b(i) = f(x(i))
	    do i = 2,imax-1
c  Unstable version
	      b(i) = - (rt/2.d0*U(i-1) - rt*U(i) + rt/2.d0*U(i+1))
     .             - dt*dsin( (v(i)+U(i))/2.d0 ) 
     .             - dt/(16*dx**2)*(v(i+1)-v(i-1)+U(i+1)-U(i-1))**2

c  Stable version
c	      b(i) = - (rt/2.d0*U(i-1) + (1-rt)*U(i) + rt/2.d0*U(i+1))
c     .             - dt*dsin( (v(i)+U(i))/2.d0 ) 
c     .             - dt/(16*dx**2)*(v(i+1)-v(i-1)+U(i+1)-U(i-1))**2
	    end do
	    i = imax
	    b(i) = g(x(i))

c  Solve the system
	    call tridiag(L,M,R,b,output,imax)

c  Calculate error and update for next iteration
	    err = 0
	    do i = 1,imax
	      err = err + (v(i)-output(i))**2
	    end do
	    err = sqrt(err)
	    write(6,*) 'err = ', err

	    do i = 1,imax
	      v(i) = output(i)
	    end do

	  mval = mval + 1

	write(2,*) err, tol
	write(2,*) mval, mmax

	  end do

	  do i = 1,imax
	    U(i) = v(i)
	  end do

c  Output
c	  if (np*((n+1)/np) .eq. n+1) then
	    call ooutfn('smp2',chan3((n+1)/np),cflag,1)
	    write(6,*) n+1,'    ',chan3((n+1)/np)

	    do i = 1,imax
	      write(1,*) x(i),U(i)
	    end do

	    close(1)
c	  end if

	end do

c  Output

c	open(9,file = 'hw2output.txt')
c	do i = 1,imax
c	  write(9,*) x(i),U(1)
c	end do
c	close(9)

c	do i = 1,10
c	  write(6,*) x(i),U(1)
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