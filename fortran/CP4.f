      program CP4

c  Adam Martinez
c  22 APR 2010
c  3450:730 ANSPDE
c  Computer Project 4

c  Solve ut + F(u)x = 0 st
c  u(x,tmin) = s(x)
c  u(xmin,t) = l(t)
c  uxx(xmax,t) = r(t)

      implicit none
      double precision u(1001,0:1),x(1001),t(24001),uip1n,
     .dx,dt,xmin,xmax,tmin,tmax,s,l,r,f,dfdu,alpha,beta,am,ap,k,
     .ri,rip1,rim1,urm,urp,ulm,ulp,a,b,c,fp,fm,phi
      integer i,imax,n,nmax,m

      character*1 cflag
      character*3 chan3
      integer np

      write(6,*) 'Enter id flag (1 character)'
      read(5,77) cflag
77    format(a1)

      dx   = 1.d-3
      dt   = 1.d-5
	xmin = 0.d0
	xmax = 1.d0
	tmin = 0.d0
	tmax = 2.d-1

	M = (tmax-tmin)/dt + 1.5d0

      np = (M-1)/5

      imax = (xmax-xmin)/dx + 1.5d0
      nmax = (tmax-tmin)/dt + 1.5d0

      write(6,*) 'dx   = ',dx
      write(6,*) 'dt   = ',dt
      write(6,*) 'xmin = ',xmin
      write(6,*) 'xmax = ',xmax
      write(6,*) 'tmin = ',tmin
      write(6,*) 'tmax = ',tmax
      write(6,*) 'imax = ',imax
      write(6,*) 'nmax = ',nmax
      write(6,*) 'M    = ',M
      write(6,*) 'np   = ',np

c   Construct spatial grids
      do i=1,imax
         x(i) = xmin + dfloat(i-1)*dx
      end do
      do n=1,nmax
         t(n) = tmin + dfloat(n-1)*dt
      end do

c   Apply IC
      n = 0
      call ooutfn('smp2',chan3(n),cflag,1)
      do i=1,imax
         u(i,0) = s(x(i))
         write(1,*) x(i),u(i,0)
      end do
      close(1)

c   main loop
c   u(i,0) is old value u_i^k
c   u(i,1) is new value u_i^{k+1}
      do n=1,nmax-1

c  Dirichlet condition at i = 1
	   i = 1
         u(i,1) = l(t(n+1))

c  Lax-Wendroff at i = 2
         i = 2
	      fp = u(i,0)
		ap = fp * (u(i,0)+u(i+1,0))/2.d0
		am = fp * (u(i,0)+u(i-1,0))/2.d0
		alpha = (dt/dx)/2.d0 * ( f(u(i+1,0)) - f(u(i-1,0)) )
		beta = ((dt/dx)**2)/2.d0 * (ap*(f(u(i+1,0))-f(u(i,0)))
     . - am*(f(u(i,0))-f(u(i-1,0))) )
            u(i,1) = u(i,0) - alpha + beta
         
c  MUSCL method for i = 3,imax-2
	   do i = 3,imax-2
		ri = (u(i,0)-u(i-1,0))/(u(i+1,0)-u(i,0))
		rip1 = (u(i+1,0)-u(i,0))/(u(i+2,0)-u(i+1,0))
		rim1 = (u(i-1,0)-u(i-2,0))/(u(i,0)-u(i-1,0))
		urm = u(i  ,0) - phi(ri  )*(u(i+1,0)-u(i  ,0))/2.d0
		urp = u(i+1,0) - phi(rip1)*(u(i+2,0)-u(i+1,0))/2.d0
		ulm = u(i-1,0) + phi(rim1)*(u(i  ,0)-u(i-1,0))/2.d0
		ulp = u(i  ,0) + phi(ri  )*(u(i+1,0)-u(i  ,0))/2.d0
		a = dfdu(u(i-1,0))
		b = dfdu(u(i  ,0))
		c = dfdu(u(i+1,0))
		am = max(abs(a),abs(b))
		ap = max(abs(b),abs(c))
		fp = (f(urp) + f(ulp) - ap*(urp-ulp))/2.d0
		fm = (f(urm) + f(ulm) - am*(urm-ulm))/2.d0
		u(i,1) = u(i,0) - dt/dx * (fp - fm)
	   end do

c  Lax-Wendroff at i = imax-1, imax
	   i = imax-1
	      fp = u(i,0)
		ap = fp * (u(i,0)+u(i+1,0))/2.d0
		am = fp * (u(i,0)+u(i-1,0))/2.d0
		alpha = (dt/dx)/2.d0 * ( f(u(i+1,0)) - f(u(i-1,0)) )
		beta = ((dt/dx)**2)/2.d0 * (ap*(f(u(i+1,0))-f(u(i,0)))
     . - am*(f(u(i,0))-f(u(i-1,0))) )
            u(i,1) = u(i,0) - alpha + beta
         	
	   i = imax
	   uip1n = 2.d0*u(i,0) - u(i-1,0)
		fp = u(i,0)
		ap = fp * (u(i,0)+uip1n)/2.d0
		am = fp * (u(i,0)+u(i-1,0))/2.d0
		alpha = (dt/dx)/2.d0 * ( f(uip1n) - f(u(i-1,0)) )
		beta = ((dt/dx)**2)/2.d0 * (ap*(f(uip1n)-f(u(i,0)))
     . - am*(f(u(i,0))-f(u(i-1,0))) )
            u(i,1) = u(i,0) - alpha + beta

         if (np*(n/np) .eq. n) then
            call ooutfn('smp2',chan3(n/np),cflag,1)
c           write(6,*) n,'   ',chan3(n/np),cflag
            do i=1,imax
               write(1,*) x(i),u(i,1)
            end do
            close(1)
         end if
         do i=1,imax
            u(i,0) = u(i,1)
         end do
      end do

      stop
      end

      function s(x)
      implicit none
      double precision s,x
        s = exp(-1.d2*(x-3.d-1)**2)
      return
      end

      function l(t)
      implicit none
      double precision l,t
        l = 0
      return
      end

	function r(t)
	implicit none
	double precision r,t
	  r = 0
	return
	end

	function f(u)
	implicit none
	double precision f,u
	  f = u**2/2.d0
	return
	end

	function dfdu(u)
	implicit none
	double precision dfdu,u
	  dfdu = u
	return
	end

	function phi(r)
	implicit none
	double precision phi,r
	  if (r.le.0) then
	    phi = 0
	  else
	    phi = 2.d0*r/(r**2+1.d0)
	  end if
	return
	end

      subroutine ooutfn(fpr,fid,cflag,lun)
c  usage
c      type *,' Enter file id'
c      accept (2a), id
c      call ooutfn(id,'i',1) creates output file i'id'.dat
c      taken from Nazanin Imani

      character*3 fid
      character*1 cflag
      character*4 fpr
      integer lun
      character*20 flname

      flname = fpr // cflag // fid
      open(unit=lun,file=flname,err=150)
      return
150     write(6,*) '*** ERROR IN OPENING FILE ***'
      return
      end

      function chan3(M)
c     Stephen Cardarelli, Oct 2003
c     provides time step files 000 to 999
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




      
