#include "gpsConversion.h"

DMSData GpsDDToDMS(double latDD, double lonDD){
	//Degrees = Whole part of number
	float latDeg = (int8_t )latDD;
	float lonDeg = (int8_t)lonDD;

	//Minutes - multiply remaining value by 60, take whole part
	float latMin = (int8_t)((latDD - latDeg) * 60);
	float lonMin = (int8_t)((lonDD - lonDeg) * 60);

	//Seconds - multiply new remaining dec by 60 (how much precession do we want)?
	float latSec = (((latDD - latDeg) * 60) - latMin) * 60;
	float lonSec = (((lonDD - lonDeg) * 60) - lonMin) * 60;

	//TODO: Figure out how/when to store this
	//DMS format latDeg° latMin' latSec"

	struct DMSData retData;
	retData.latDegree = latDeg;
	retData.lonDegree = lonDeg;
	retData.latMinute = latMin;
	retData.lonMinute = lonMin;
	retData.latSecond = latSec;
	retData.lonSecond = lonSec;

	return retData;
}

//Convert from gps DD to UTM
//TODO: For now, all is float cause I don't know what I'm doing and converting this code from JS, so come back and make this better
//Original JS thanks to: www,movable-type.co.uk/scripts/latlong-utm-mgrs.html

//Values look correct after hand-computing everything.
UTMData GpsDDtoUTM(double latDD, double lonDD){
	//Assume lat/long are valid. If they're not we should already know because the accuracy value will be cracked

	//Establish easting and northings to ref
	double falseEasting = 500e3, falseNorthing = 10000e3;


	//Figure out zone & lat band
	uint8_t zone = floor((lonDD+180)/6) + 1; //Longitudinal zone;
	double lambda0 = ((zone-1)*6 - 180 + 3)*DEG2RAD; //Longitude of central meridian

	//Handle Norway/Svalbard exceptions
	char mgrsLatBand[21] = "CDEFGHJKLMNPQRSTUVWXX"; //X is repeated for 80-84degN
	char latBand = mgrsLatBand[(int)floor(latDD/8+10)];

	//adjust zone & central meridian for Norway
	if(zone == 31 && latBand=='V' && lonDD>=3){
		zone++;
		lambda0 += 6*DEG2RAD;
	}
	//adjust zone & central meridian for Svalbard
	if(zone==32 && latBand=='X' && lonDD < 9){zone--; lambda0 -= 6*DEG2RAD;}
	if(zone==32 && latBand=='X' && lonDD >= 9){zone++; lambda0 += 6*DEG2RAD;}

	if(zone==34 && latBand=='X' && lonDD < 21){zone--; lambda0 -= 6*DEG2RAD;}
	if(zone==34 && latBand=='X' && lonDD >= 21){zone++; lambda0 += 6*DEG2RAD;}

	if(zone==36 && latBand=='X' && lonDD < 33){zone--; lambda0 -= 6*DEG2RAD;}
	if(zone==36 && latBand=='X' && lonDD >= 33){zone++; lambda0 += 6*DEG2RAD;}

	double phi = latDD*DEG2RAD;
	double lambda = lonDD*DEG2RAD - lambda0;

	//USE WGS84 ellipsoid
	// unused - everything is hand done, but for ref --
	//double a = 6378137;
	//double f = 1/298.257223563;

	double k0 = 0.9996; //UTM scale on central meridian

	//Compute - easting, northing: Karney Eq 7-14, 29, 35:

	//hand compute e [f*(2-f)] -- eccentricity
	double e = 0.08181919084;

	//n is used in computation of alpha and A, hand compute those, so not needed, keeping for future debugging
	//double n = f / (2-f); //3-d flattening

	//So we don't have to re-compute
	double cosLambda = cos(lambda);
	double sinLambda = sin(lambda);

	double tau = tan(phi);
	double sigma = sinh(e*atanh(e * tau/sqrt(1+tau*tau)));

	double tauPrime = tau*sqrt(1+sigma*sigma) - sigma*sqrt(1 + tau*tau);

	double zetaPrime = atan2(tauPrime, cosLambda);
	double etaPrime = asinh(sinLambda / sqrt(tauPrime*tauPrime + cosLambda*cosLambda));

	//Constant value - hand compute - leave eq as reference - note, n4 is n^4 and so on
//	double A = a/(1+n)*(1 + 1/4*n2 + 1/64*n4 + 1/256*n6); //2piA is the circumerence of a meridian
	double A = 6367449.146;

	 //note alpha is a one-based array (6th order Kruger expressions)
//	double alpha[7] = {0,
//				   1/2*n - 2/3*n2 + 	5/16*n3 + 	41/180*n4 - 	127/288*n5 +      7891/37800*n6,
//				   	   	 13/48*n2 -	     3/5*n3 + 557/1440*n4 +	    281/630*n5 - 1983433/1935360*n6,
//						   	   	   	  61/240*n3 -	103/140*n4 + 15061/26880*n5 +   167603/181440*n6,
//											  49561/161280*n4 -     179/168*n5 + 6601661/7257600*n6,
//											  	  	  	  	    34729/80640*n5 - 3418889/1995840*n6,
//															  	  	  	  	 212378941/319334400*n6};

	//alpha array hand compute, leave the above comment as reference
	double alpha[7] = {
			0,
			8.37731821e-4,
			7.60852777e-7,
			1.1976455e-9,
			2.42917061e-12,
			5.71175768e-15,
			1.49111773e-17
	};

	double zeta = zetaPrime;
	for(int j = 1; j<=6;j++){
		zeta+= alpha[j] * sin(2*j*zetaPrime) * cosh(2*j*etaPrime);
	}

	double eta = etaPrime;
	for(int j = 1; j <=6; j++){
		eta += alpha[j] * cos(2*j*zetaPrime) * sinh(2*j*etaPrime);
	}

	double x = k0 * A * eta;
	double y = k0 * A * zeta;

//Compute Convergance & Scale here, we're not going to cause we don't care...

	//Shift x/y to false origins
	x = x+falseEasting; //make x relative to false easting
	if(y < 0){
		y = y + falseNorthing; //make y in southern hemisphere relative to false northing
	}

	//round to reasonable precision.
	//The joys of C; //round to 9 decimal places
	x = round(x*1000000000)/1000000000;
	y = round(y*1000000000)/1000000000;

	char h = latDD >= 0 ? 'N' : 'S';

	//return zone, latBand, h, x, y
	struct UTMData retData;
	retData.easting = x;
	retData.northing = y;
	retData.zone = zone;
	retData.latBand = latBand;
	retData.hemisphere = h;

	return retData;
}
