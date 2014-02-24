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
     .x(61),y(151),t(1001),
     .AAL(61),AAM(61),AAR(61),AB(61),
     .BAL(151),BAM(151),BAR(151),BB(151),
     .tmin,tmax,xmin,xmax,ymin,ymax,
     .dt,dx,dy,A,B,C,
     .G,Q,F,H,P,
     .rx,ry,bx,by,ct,alpha,beta,
     .xoutput(61),youtput(151),S
      integer nmax,n,imax,i,jmax,j

	character*1 cflag
	character*13 flname
	character*3 chan3
	integer np

c  U - function value at time t(n)
c  V - function value at time t(n+1/2)
c  W - function value at time t(n+1)
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
	dy = (2.d0-0.d0)/1.5d2

      tmin = 0.d0
	tmax = 1.d0
      xmin = 0.d0
      xmax = 1.d0
	ymin = 0.d0
	ymax = 2.d0

      nmax = (tmax-tmin)/dt + 1.5d0
      imax = (xmax-xmin)/dx + 1.5d0
	jmax = (ymax-ymin)/dy + 1.5d0

	rx = alpha*dt/dx**2
c	ry = alpha*dt/dy**2
	ry = 28.125d0
c	bx = dt/(2.d0*dx)
	bx = 3.d-2
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
	  y(j) = ymin + dfloat(j-1)*dy
	end do

	write(6,*) 't(nmax) = ',t(nmax)
	write(6,*) 'x(imax) = ',x(imax)
	write(6,*) 'y(jmax) = ',y(jmax)

c  Apply initial condition

      n = 0
	call ooutfn('smp2',chan3(n),cflag,1)
	do i = 1,imax
	  do j = 1,jmax
	    U(i,j) = P(x(i),y(j))
	    write(1,*) x(i),y(j),U(i,j)
	  end do
      end do
	
c  Construct the matricies
c  First row
	  i = 1
	  AAL(i) = 0.d0
	  AAM(i) = 1.d0 + rx - Ct/2.d0 
     . + (2.d0*dx/beta)*(-rx/2.d0 + A*bx/2.d0)
	  AAR(i) = -rx
	  j = 1
	  BAL(j) = 0.d0
	  BAM(j) = 1.d0
	  BAR(j) = 0.d0

c  Interior rows
	  do i = 2,imax-1
	    AAL(i) =      - rx/2.d0 + A*bx/2.d0
	    AAM(i) = 1.d0 + rx      -   Ct/2.d0
	    AAR(i) =      - rx/2.d0 - A*bx/2.d0
	  end do
	  do j = 2,jmax-1
	    BAL(j) =      - ry/2.d0 + B*by/2.d0
	    BAM(j) = 1.d0 + ry      -   Ct/2.d0
  	    BAR(j) =      - ry/2.d0 - B*by/2.d0
	  end do

c  Last row
	  i = imax
	  AAL(i) = 0.d0
	  AAM(i) = 1.d0
	  AAR(i) = 0.d0
	  j = jmax
	  BAL(j) = 0.d0
	  BAM(j) = 1.d0
	  BAR(j) = 0.d0

c  Primary loop

	write(6,*) 'Beginning loop...'
	do n = 1,nmax-1
c	write(6,*) 'Iteration ', n

c  Find V(i,j)
	j = 1
	  i = 1
	    V(i,j) = G(x(i),t(n) + dt/2.d0)
	  do i = 2,imax-1
	    V(i,j) = G(x(i),t(n) + dt/2.d0)
	  end do
	  i = imax
	    V(i,j) = (G(x(i),t(n)+dt/2.d0)+H(y(j),t(n)+dt/2.d0))/2.d0

	do j = 2,jmax-1
	  i = 1
	    AB(i) = U(i,j) + ry/2.d0*(U(i,j+1)-2*U(i,j)+U(i,j-1))
     . + B*by*(U(i,j+1)-U(i,j-1)) + Ct/2.d0*U(i,j) 
     . + dt/2.d0*S(x(i),y(j),t(n))
     . + 2*dt/beta*(-rx/2.d0 + A*bx/2.d0)*F(y(j),t(n)+dt/2.d0)

	  do i = 2,imax-1
	    AB(i) = U(i,j) + ry/2.d0*(U(i,j+1)-2*U(i,j)+U(i,j-1))
     . + B*by*(U(i,j+1)-U(i,j-1)) + Ct/2.d0*U(i,j) 
     . + dt/2.d0*S(x(i),y(j),t(n))
	end do

        i = imax
	    AB(i) = H(y(j),t(n)+dt/2.d0)

c  Solve the system
	  call tridiag(AAL,AAM,AAR,AB,xoutput,imax)

	  do i = 1,imax
	    V(i,j) = xoutput(i)
	  end do

	end do

	j = jmax
	  i = 1
	    V(i,j) = Q(x(i),t(n)+dt/2.d0)
	  do i = 2,imax-1
	    V(i,j) = Q(x(i),t(n)+dt/2.d0)
	  end do
	  i = imax
          V(i,j) = (Q(x(i),t(n)+dt/2.d0)+H(y(j),t(n)+dt/2.d0))/2.d0

c  Find W(i,j)

	i = 1
	  j = 1
	    BB(j) = G(x(i),t(n)+dt)

	  do j = 2,jmax-1
	    BB(j) = (V(i+1,j)+(2*dx/beta)*(V(i,j)-F(y(j),t(n)+dt/2.d0)))
     . * (rx/2.d0 - A*bx/2.d0)
     . + V(i,j)*(1 - rx + Ct/2.d0) + V(i+1,j)*(rx/2.d0 + A*bx/2.d0)
     . + dt/2.d0*S(x(i),y(j),t(n)+dt)
	  end do

	  j = jmax
	    BB(j) = Q(x(i),t(n)+dt)

c  Solve the system
	  call tridiag(BAL,BAM,BAR,BB,youtput,jmax)

	  do j = 1,jmax
	    W(i,j) = youtput(j)
	  end do

	do i = 2,imax-1
	  j = 1
	    BB(j) = G(x(i),t(n)+dt)

	  do j = 2,jmax-1
	    BB(j) = V(i-1,j)*(rx/2.d0 - A*bx/2.d0)
     . + V(i,j)*(1 - rx + Ct/2.d0) + V(i+1,j)*(rx/2.d0 + A*bx/2.d0)
     . + dt/2.d0*S(x(i),y(j),t(n)+dt)
	end do

        j = jmax
	    BB(j) = Q(x(i),t(n)+dt)

c  Solve the system
	  call tridiag(BAL,BAM,BAR,BB,youtput,jmax)

	  do j = 1,jmax
	    W(i,j) = youtput(j)
	  end do

	end do

	i = imax
	  j = 1
	    W(i,j) = (G(x(i),t(n)+dt) + H(y(j),t(n)+dt))/2.d0
	  do j = 2,jmax-1
	    W(i,j) = H(y(j),t(n)+dt)
	  end do
	  j = jmax
	    W(i,j) = (Q(x(i),t(n)+dt) + H(y(j),t(n)+dt))/2.d0	    

c  Copy W(i,j) into U(i,j) for next iteration.
	do i = 1,imax
	  do j = 1,jmax
	    U(i,j) = W(i,j)
	  end do
	end do

c  Output
	  if (np*((n+1)/np) .eq. n+1) then
	    call ooutfn('smp2',chan3((n+1)/np),cflag,1)
	    write(6,*) n+1,'    ',chan3((n+1)/np)

	    do i = 1,imax
		do j = 1,jmax
	        write(1,*) x(i),y(j),U(i,j)
	      end do
	    end do

	    close(1)
	  end if

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

	function Q(x,t)
      implicit none
      double precision Q,x,t
        Q = x + dasin(t)/1.d1
      return
      end

	function G(x,t)
      implicit none
      double precision G,x,t
        G = x + dasin(t)/1.d1
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
        H = 1.d0 - y*(2.d0-y) + dasin(t)/1.d1
      return
      end

      function P(x,y)
      implicit none
      double precision P,x,y
	  P = 0
      return
      end

	function S(x,y,t)
	implicit none
	double precision S,x,y,t
	  S = x*y*exp(-t)
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