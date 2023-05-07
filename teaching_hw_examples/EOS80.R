###############################################################################
#
# Inputs
#     temp: vector of water temperature values (units: degrees Celcius)
#     sal:  vector of salinity values (units: practical salinity units)
#     pres: vector of depth (units: decibars)
#
# Output
#     vector of water density (units: kg/m^3)
#
# Note that length(temp) == length(sal) == length(pres)
#
###############################################################################

water.density <- function(temp, sal, pres){
	
	st0 = 999.842594 + 6.793952E-2 * temp - 9.095290E-3 * temp^2 + 
	      1.001685E-4 * temp^3 - 1.120083E-6 * temp^4 + 6.536332E-9 * 
	      temp^5 + 8.24493E-1 * sal - 4.0899E-3 * temp * sal + 7.6438E-5 * 
	      temp^2 * sal - 8.2467E-7 * temp^3 * sal + 5.3875E-9 * temp^4 * sal -
	      5.72466E-3 * sal^(3/2) + 1.0227E-4 * temp * sal^(3/2) - 1.6546E-6 *
	      temp^2 * sal^(3/2) + 4.8314E-4 * sal^2

	K = 19652.21 + 148.4206 * temp - 2.327105 * temp^2 + 1.360477E-2 * temp^3 -
	    5.155288E-5 * temp^4 + 3.239908 * pres + 1.43713E-3 * temp * pres + 
	    1.16092E-4 * temp^2 * pres - 5.77905E-7 * temp^3 * pres + 8.50935E-5 *
	    pres^2 - 6.12293E-6 * temp * pres^2 + 5.2787E-8 * temp^2 * pres^2 + 
	    54.6746 * sal - 0.603459 * temp * sal + 1.09987E-2 * temp^2 * sal -
	    6.167E-5 * temp^3 * sal + 7.944E-2 * sal^(3/2) + 1.6483E-2 * temp *
	    sal^(3/2) - 5.3009E-4 * temp^2 * sal^(3/2) + 2.2838E-3 * pres * sal -
	    1.0981E-5 * temp * pres * sal - 1.6078E-6 * temp^2 * pres * sal + 
	    1.91075E-4 * pres * sal^(3/2) - 9.9348E-7 * pres^2 * sal + 2.0816E-8 *
	    temp * pres^2 * sal + 9.1697E-10 * temp^2 * pres^2 * sal
	
	stp = st0 / (1 - (pres / K))
	
	return(stp)
	}
