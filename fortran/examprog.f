c  Adam Martinez
c  03 MAY 2010
c  3450:730 ANSPDE
c  Final Exam

c  Numerically solve the diffential equation, Ut = U^2*Ux - 20Ux
c  U( 0,t) = 0
c  U(10,t) = 0
c  U( x,0) = P(x,y)
c  0 < x < 1
c  0 < y < 2
c  0 < t < inf

      implicit none
      double precision U(10001),V(10001),W(10001),
     .x(10001),t(50001),tmin,tmax,xmin,xmax,dt,dx,F,P
      integer nmax,n,imax,i

	character*1 cflag
	character*13 flname
	character*3 chan3
	integer np

c  xmin,xmax - endpoints of x domain
c  tmin,tmax - endpoints of time domain
c  dt,dx - finite differences of domains
c  P - initial condition
c  nmax - number of points in time grid
c  imax - number of points in x grid

	write(6,*) 'Enter ID flag (1 character)'
	read(5,77) cflag
77	format(a1)

	np = 10000

      dx = 1.d-3
      dt = dx/1.d3

      tmin = 0.d0
	tmax = 5.d-2
	xmin = 0.d0
	xmax = 1.d1

      nmax = (tmax-tmin)/dt + 1.5d0
      imax = (xmax-xmin)/dx + 1.5d0

      write(6,*) 'dt   = ',dt
      write(6,*) 'dx   = ',dx
      write(6,*) 'nmax = ',nmax
      write(6,*) 'imax = ',imax

c  Generate t, x meshes

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
	  U(i) = P(x(i))
        write(1,*) x(i),U(i)
      end do
	
c  Primary loop
c	write(6,*) 'Beginning loop...'
	do n = 1,nmax-1
c  First stage
        do i = 1,imax
	    v(i) = u(i)-dt/dx*(F(U(i+1)) - F(U(i)))
	  end do

c  Second stage
	  do i = 1,imax
	    w(i) = (U(i)+V(i))/2.d0 
     . - dt/(2.d0*dx)*(F(V(i))-F(V(i-1)))
	  end do
	
c  Reset stages
	  do i = 1,imax
	    u(i) = w(i)
	  end do

c  Output
	  if (np*((n+1)/np) .eq. n+1) then
	    call ooutfn('smp2',chan3((n+1)/np),cflag,1)
	    write(6,*) n+1,'    ',chan3((n+1)/np)

	    do i = 1,imax
	      write(1,*) x(i),U(i)
	    end do

	    close(1)
	  end if

	end do

      stop
      end

      function P(x)
      implicit none
      double precision P,x
	  P = 3*exp(-5.d0*(x-3.d0)**2) + 2*exp(-1.d1*(x-5.d0)**2)
      return
      end

      function F(U)
      implicit none
      double precision F,U
	  F = 20.d0*U-U**3/3.d0
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