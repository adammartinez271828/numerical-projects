c  Adam Martinez
c  03 MAR 2010
c  3450:635 Optimization
c  Project #3

	implicit none
	double precision x(2),xp(2),alpha,c1,c2,o(2),
     .tol,err,gtol,gerr,LHS,RHS,g(2),p(2),grad(2)
	integer m,mmax

c  mmax is the maximal number of steps being allowed
c  tol is the desired precision

	write(6,*) "Input the starting point."
	read(5,*) x(1),x(2)
	o = x

	alpha = 1.d0
	c1 = 1.d-4
	c2 = 9.d-1
	tol = 1.d-12
	gtol = 1.d-12
	mmax = 50000
	xp = x

	m = 0
	err = 2*tol
	gerr = 2*gtol

	do while(((err.ge.tol).or.(gerr.ge.gtol)).and.(m.lt.mmax))

c  Compute LHS
	  call gradf(x,p)
	  p(1) = -p(1)
	  p(2) = -p(2)
	  call f(x + alpha*p,LHS)

c  Compute RHS
	  RHS = p(1)*(-p(1)) + (-p(2))*p(2)
	  call f(x+c1*alpha*RHS,RHS)

c  Check to see if step size is acceptable
	  if(LHS.le.RHS) then
c  If so, make step and prepare for next iteration
	    xp = x
	    x = x + alpha*p
	    RHS = LHS
	    call gradf(x,g)
	    call gradf(x,p)
	    p(1) = -p(1)
	    p(2) = -p(2)
	    m = m + 1
	    alpha = 1.d0
c  Calculate error of new point
	    err = sqrt((x(1)-xp(1))**2+(x(2)-xp(2))**2)
	    gerr = sqrt(g(1)**2+g(2)**2)
c	    write(6,*) 'Point: [',x,']'
c	    write(6,*) 'Delta(x): ',err
c	    write(6,*) 'Magnitude of Gradient: ',gerr

c  Otherwise, reduce alpha by half
	  else
	    alpha = alpha/2.d0
	  end if	

	end do

	open(unit=1,file="opproj3.txt")

	if(m.ge.mmax) then
	  write(1,*) 'Maximum number of iterations exceeded.'
	end if

	write(1,*) 'Steps taken: ',m
	write(1,*) 'Origin: [',o,']'
	write(1,*) 'Stopping point: [',x,']'
	write(1,*) 'Delta(x): ',err
	write(1,*) 'Magnitude of Gradient: ',gerr

	stop
	end

	subroutine f(x,output)
	implicit none
	double precision x(*),output
	output = 100*(x(2)-x(1)**2)**2 + (1.d0 - x(1))**2
	end subroutine f

	subroutine gradf(x,output)
	implicit none
	double precision x(*),output(*)
	output(1) = 400.d0*x(1)**3 + 2*x(1) - 400.d0*x(1)*x(2) - 2
	output(2) = 200.d0*x(2)-200.d0*(x(1)**2)
	end subroutine gradf