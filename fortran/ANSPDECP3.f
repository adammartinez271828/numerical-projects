c  Adam Martinez
c  20 APR 2010
c  3450:730 ANSPDE
c  Computer Project 3

c  Numerically solve the diffential equation, Ut = alpha*(Uxx+Uyy) + AUx + BUy + CU + S(x,y,t)
c  U(x,0,t) = G(x,t)
c  U(x,2,t) = Q(x,t)
c  U(0,y,t) + beta*Ux(0,y,t) = F(y,t)
c  U(1,y,t) = H(y,t)
c  U(x,y,0) = P(x,y)
c  0 < x < 1
c  0 < y < 2
c  0 < t < inf

      implicit none
      double precision U(61,151),V(61,151),W(61,151),
     .x(61),y(151),t(1001),AAL(61),AAM(61),AAR(61),
     .BAL(151),BAM(151),BAR(151),AB(61),BB(151),
     .tmin,tmax,xmin,xmax,dt,dx,dy,
     .G,Q,F,H,P,rx,ry,bx,bt,ct,alpha,beta
      integer nmax,n,imax,i,jmax,j

	character*1 cflag
	character*13 flname
	character*3 chan3
	integer np

c  xmin,xmax - endpoints of x domain
c  ymin,ymax - endpoints of y domain
c  tmin,tmax - endpoints of time domain
c  dt,dx,dy - finite differences of domains
c  F,H - left and right	boundary functions
c  Q,G - top and bottom boundary functions
c  P - initial condition
c  nmax - number of points in time grid
c  imax - number of points in x grid
c  jmax - number of points in y grid

	write(6,*) 'Enter ID flag (1 character)'
	read(5,77) cflag
77	format(a1)

	np = 1000

      alpha = 5.d0
	beta = 3.d0
	A = 1.d0
	B = 5.d-1
	C = -3.d0

      dt = 1.d-3
      dx = (1.d0-0.d0)/6.d1
	dy = (2.d0-0.d0)/1.5d1

      tmin = 0.d0
	tmax = 1.d0
      xmin = 0.d0
      xmax = 1.d0
	ymin = 0.d0
	ymax = 2.d0

      nmax = (tmax-tmin)/dt + 1.5d0
      imax = (xmax-xmin)/dx + 1.5d0

	rx = alpha*dt/dx**2
	ry = alpha*dt/dy**2
	bx = dt/(2.d0*dx)
	by = dt/(2.d0*dy)
	ct = C*dt/2.d0

      write(6,*) 'dt   = ',dt
      write(6,*) 'dx   = ',dx
	write(6,*) 'dy   = ',dy
      write(6,*) 'rx   = ',rx
      write(6,*) 'ry   = ',ry
      write(6,*) 'bx   = ',bx
      write(6,*) 'by   = ',by
	write(6,*) 'ct   = ',ct
      write(6,*) 'nmax = ',nmax
      write(6,*) 'imax = ',imax
	write(6,*) 'jmax = ',jmax

c  Generate t, x, y meshes

      do n = 1,nmax
	  t(n) = tmin + dfloat(n-1)*dt
      end do

      do i = 1,imax
	  x(i) = xmin + dfloat(i-1)*dx
      end do

	do j = 1,jmax
	  y(j) = jmin + dfloat(j-1)*dy
	end do

c  Apply initial condition

      n = 0
	call ooutfn('smp2',chan3(n),cflag,1)
	do i = 1,imax
	  do j = 1,jmax
	    U(i,j) = P(x(i),y(j))
	    write(1,*) x(i),y(j),U(i,j)
	  end do
      end do
	
c  Construct the matrix
c  First row
c	  i = 1
c	  L(i) = 0.d0
c	  M(i) = 1.d0
c	  R(i) = 0.d0
c
c  Interior rows
c	  do i = 2,imax-1
c	    L(i) = rt/2.d0
c  Unstable version
c	    M(i) = -rt
c  Stable version
c	    M(i) = -(1.d0+rt)
c	    R(i) = rt/2.d0
c	  end do
c
c  Last row
c	  i = imax
c	  L(i) = 0.d0
c	  M(i) = 1.d0
c	  R(i) = 0.d0
c
c  Primary loop
c
c	write(6,*) 'Beginning loop...'
c     do n = 1,nmax-1
c	write(6,*) 'Iteration ', n
c	  do i = 1,imax
c	    v(i) = u(i)
c	  end do
c
c	  err = 2.d0*tol
c	  mval = 0
c
c	  do while((err.gt.tol).and.(mval.le.mmax))
c
c	    i = 1
c	    b(i) = f(x(i))
c	    do i = 2,imax-1
c  Unstable version
c	      b(i) = v(i) - U(i) - dt*dsin(5.d-1*(v(i)+U(i)) 
c    .             - 2.5d-1*bt**2
c    .             * (v(i+1)-v(i-1)+U(i+1)-U(i-1))**2

c  Stable version
c	      b(i) = - (rt/2.d0*U(i-1) + (1-rt)*U(i) + rt/2.d0*U(i+1))
c     .             - dt*dsin( (v(i)+U(i))/2.d0 ) 
c     .             - dt/(16*dx**2)*(v(i+1)-v(i-1)+U(i+1)-U(i-1))**2
c	    end do
c	    i = imax
c	    b(i) = g(x(i))

c  Solve the system
c	    call tridiag(L,M,R,b,output,imax)

c  Calculate error and update for next iteration
c	    err = 0
c	    do i = 1,imax
c	      err = err + (v(i)-output(i))**2
c	    end do
c	    err = sqrt(err)
c	    write(6,*) 'err = ', err

c	    do i = 1,imax
c	      v(i) = output(i)
c	    end do

c	  mval = mval + 1

c	  end do

c	  do i = 1,imax
c	    U(i) = v(i)
c	  end do

c  Output
c	  if (np*((n+1)/np) .eq. n+1) then
c	    call ooutfn('smp2',chan3((n+1)/np),cflag,1)
c	    write(6,*) n+1,'    ',chan3((n+1)/np)

c	    do i = 1,imax
c	      write(1,*) x(i),U(i)
c	    end do

c	    close(1)
c	  end if

c	end do

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

	function Q(x,t)
      implicit none
      double precision Q,x,t
        Q = x + dasin(t)/1.d1
      return
      end

	function G(x,t)
      implicit none
      double precision G,x,t
        F = x + dasin(t)/1.d1
      return
      end

      function F(y,t)
      implicit none
      double precision F,y,t
        F = dasin(t)/1.d1
      return
      end

      function H(y,t)
      implicit none
      double precision H,y,t
        H = 1.d0 - y(2.d0-y) + dasin(t)/1.d1
      return
      end

      function P(x,y)
      implicit none
      double precision h,x,y
	  P = 0
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