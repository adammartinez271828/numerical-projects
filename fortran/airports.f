	program airports

c  Adam Martinez
c  19 JAN 2010
c  3450:635 Optimization
c  Project #1, Great Cirlce Distance

c  Recieve input

	implicit none
	integer dolat,molat,solat,dolon,molon,solon
	integer ddlat,mdlat,sdlat,ddlon,mdlon,sdlon
	double precision origlat,origlon,destlat,destlon,dlat,dlon

c  Angular distance and great circle distance

	double precision distance,dist,earthrad,pi

c  Distance in kilometers

	earthrad = 6371.0d0
	pi = 3.141592653589793238462643d0

	write(6,*) 'Enter the lat/lon of the origin. (d m s d m s)'
	read(*,*) dolat,molat,solat,dolon,molon,solon
	write(6,*) 'Enter the lat/lon of the destination. (d m s d m s)'
	read(*,*) ddlat,mdlat,sdlat,ddlon,mdlon,sdlon

	origlat = dfloat(dolat) + dfloat(molat)/60.0d0 + dfloat(solat)/360.0d0
	origlon = dfloat(dolon) + dfloat(molon)/60.0d0 + dfloat(solon)/360.0d0
	destlat = dfloat(ddlat) + dfloat(mdlat)/60.0d0 + dfloat(sdlat)/360.0d0
	destlon = dfloat(ddlon) + dfloat(mdlon)/60.0d0 + dfloat(sdlon)/360.0d0

	write(6,*) 'Origin latitude      :',origlat
	write(6,*) 'Origin longitude     :',origlon
	write(6,*) 'Destination latitude :',destlat
	write(6,*) 'Destination longitude:',destlon

	origlat = origlat*2*pi/360
	origlon = origlon*2*pi/360
	destlat = destlat*2*pi/360
	destlon = destlon*2*pi/360

	distance = dist(origlat,origlon,destlat,destlon)

	write(6,*) 'The great circle distance between the two cities is:'
	write(6,*) distance,'km'

	stop
	end

c  Distance finding function
c  Uses the Haversine formula for great circle distance
c  Input is lat/lon of origin and lat/lon of destination in radians
c  Output is in kilometers
c  Inaccurate for antinodal points

	function dist(olat,olon,dlat,dlon)
	implicit none
	double precision dist,olat,olon,dlat,dlon
	double precision pi,earthrad

	earthrad = 6371.0d0
	pi = 3.141592653589793238462643d0

	dist = dsin((olat-dlat)/2.0d0)**2
	dist = dist + dcos(olat)*dcos(dlat)*(dsin((olon-dlon)/2.0d0)**2)
	dist = 2.0d0*dasin(sqrt(dist))
	dist = earthrad*dist

	return
	end