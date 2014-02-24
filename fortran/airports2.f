	program airports

c  Adam Martinez
c  27 JAN 2010
c  3450:635 Optimization
c  Project #2, 51 State tour

c  coord contains the lat/lon of each element
c  dmatrix contains the distances between each element
c  path contains an ordered list of the optimized math
c  index is a set of flags showing which elements have been placed in the path

	implicit none
	integer dlat,mlat,slat,dlon,mlon,slon,l,m,n,o
	integer path(51),index(51)
	integer location,closest
	integer rotateupflag,rotatedownflag
	double precision tourlength,lowest,temp1,temp2,temp3,temp4
	double precision coord(51,2),dmatrix(51,51),lat,lon,infdist,dist
	double precision earthrad,pi

c  Distance in kilometers

	earthrad = 6371.0d0
	pi = 3.141592653589793238462643d0
	infdist = 1.0d300

c  Load data from files and initialize coordinate and distance matricies.

	open(1,FILE='airports.txt')

	read(1,*)
	read(1,*)

	do n = 1,51
	read(1,*) dlat,mlat,slat,dlon,mlon,slon
	lat = dfloat(dlat) + dfloat(mlat)/60.0d0 + dfloat(slat)/360.0d0
	lon = dfloat(dlon) + dfloat(mlon)/60.0d0 + dfloat(slon)/360.0d0
	coord(n,1) = lat
	coord(n,2) = lon
	end do

	close(1)

	do n = 1,51
	do m = 1,n
	if(n == m) then
	dmatrix(n,m) = infdist
	else
	dmatrix(n,m) = dist(coord(n,1),coord(n,2),coord(m,1),coord(m,2))
	dmatrix(m,n) = dist(coord(n,1),coord(n,2),coord(m,1),coord(m,2))
	end if
	end do
	end do

c  Construct a quick tour.

c  Initialize last location to starting point

	do n = 1,49
	path(n) = 0
	index(n) = 0
	end do
	path(51) = 51
	index(51) = 1

c  Fill rest of tour

	do n = 1,50

	lowest = infdist
	location = -1
	closest = -1

c  Find nearest to current path
c  Choose an added point and check distance to all unadded points

	do m = 1,51
	if(index(m) == 1) then
	do l = 1,51
	if((index(l)==0).AND.(dmatrix(m,l)<lowest)) then
	lowest = dmatrix(m,l)
	location = m
	closest = l
	end if
	end do
	end if
	end do

c	write(6,*) 'Adding node ',closest
c	write(6,*) 'Above node ',location
c	write(6,*) 'With distance ',lowest

c  Add closest node to path by inserting before location

	rotateupflag = 0
	rotatedownflag = 0

c  Attempt to rotate up

	do m = location-1,1,-1
	if(path(m) == 0) then
c  Rotate up and insert above location
	do l = m,location-2
	path(l) = path(l+1)
	end do
	path(location-1) = closest
	index(closest) = 1
	rotatedownflag = 1
c	write(6,*) 'Inserted node up'
	exit
	end if
	end do

c  Otherwise rotate down
	if(rotatedownflag == 0) then
	do m = location+1,51
	if(path(m) == 0) then
c  Rotate down and insert at location
	do l = m,location+1,-1
	path(l) = path(l-1)
	end do
	path(location) = closest
	index(closest) = 1
	rotateupflag = 1
c	write(6,*) 'Inserted node down'
	exit
	end if
	end do
	end if

c	if((rotatedownflag == 0).AND.(rotateupflag == 0)) then
c	write(6,*) 'Error: Cannot insert node into path.'
c	end if

	end do

c  Swap algorithm

	do o = 1,20
	do n = 1,51
	do m = 1,51
	if(n /= m) then
	temp1 = dmatrix(path(n-1),path(n))+dmatrix(path(n),path(n+1))
	temp2 = dmatrix(path(m-1),path(m))+dmatrix(path(m),path(m+1))
	temp3 = dmatrix(path(n-1),path(m))+dmatrix(path(m),path(n+1))
	temp4 = dmatrix(path(m-1),path(n))+dmatrix(path(n),path(m+1))
	if((temp1+temp2-temp3-temp4)>0) then
	l  = path(n)
	path(n) = path(m)
	path(m) = l
	end if
	end if
	end do
	end do
	end do

c  Write out optimized path

	open(1,FILE='path.txt')
	write(1,*) 'Quick tour path.'
	do n = 1,51
	write(1,*) path(n)
	end do
	close(1)

	tourlength = 0
	do n = 1,49
	tourlength = tourlength + dmatrix(path(n),path(n+1))
	end do
	tourlength = tourlength + dmatrix(path(51),path(1))

	write(6,*) 'Optimal tour is ',tourlength,'km.'

	stop
	end

c  Great circle distance finding function.
c  Uses the Haversine formula for great circle distance.
c  Input is lat/lon of origin and lat/lon of destination in degrees (double precision).
c  Output is in kilometers.
c  Inaccurate for antipodal points.

	function dist(olat,olon,dlat,dlon)
	implicit none
	double precision dist,olat,olon,dlat,dlon
	double precision ola,olo,dla,dlo
	double precision pi,earthrad

	earthrad = 6371.0d0
	pi = 3.141592653589793238462643d0

	ola = olat*2*pi/360
	olo = olon*2*pi/360
	dla = dlat*2*pi/360
	dlo = dlon*2*pi/360

	dist = dsin((ola-dla)/2.0d0)**2
	dist = dist + dcos(ola)*dcos(dla)*(dsin((olo-dlo)/2.0d0)**2)
	dist = 2.0d0*dasin(sqrt(dist))
	dist = earthrad*dist

	return
	end