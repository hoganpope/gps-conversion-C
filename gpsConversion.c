#ifndef SRC_GPS_H_
#define SRC_GPS_H_

extern "C" {


	//----------------------CONSTANTS--------------------
	static const double PI = 3.14159265358979323e0;  /* PI */
	static const double DEG2RAD = PI / 180;
	static const double RAD2DEG = 180 / PI;

	//---------------------Structs----------------------
	typedef struct UTMData {
		double easting;
		double northing;
		uint8_t zone;
		char latBand;
		char hemisphere;
	} UTMData;

	typedef struct DMSUnit {
		double degree;
		double minute;
		double second;
	} DMSUnit;


	typedef struct DMSCoord {
		DMSUnit latitude;
		DMSUnit longitude;
	} DMSCoord;

	typedef struct DDCoord {
		double latitude;
		double longitude;
	};


	//-------------------Prototypes-----------------------

	//Convert the Gps coords to minutes
	DMSCoord GpsDDToDMS(double latDD, double lonDD);

	//Convert DD to UTM
	UTMData GpsDDtoUTM(double latDD, double lonDD);

	//Convert UTM to DD
	DDCoord GpsUTMtoDD(UTMData);

	//------------------Definitions-----------------------

	DMSCoord GpsDDToDMS(double latDD, double lonDD) {
		//Degrees = Whole part of number
		double latDeg = (int8_t)latDD;
		double lonDeg = (int8_t)lonDD;

		//Minutes - multiply remaining value by 60, take whole part
		double latMin = (int8_t)((latDD - latDeg) * 60);
		double lonMin = (int8_t)((lonDD - lonDeg) * 60);

		//Seconds - multiply new remaining dec by 60 (how much precession do we want)?
		double latSec = (((latDD - latDeg) * 60) - latMin) * 60;
		double lonSec = (((lonDD - lonDeg) * 60) - lonMin) * 60;

		//TODO: Figure out how/when to store this
		//DMS format latDeg° latMin' latSec"

		DMSUnit lat;
		lat.degree = latDeg;
		lat.minute = latMin;
		lat.second = latSec;

		DMSUnit lon;
		lon.degree = lonDeg;
		lon.minute = lonMin;
		lon.second = lonSec;

		DMSCoord retData;
		retData.latitude = lat;
		retData.longitude = lon;

		return retData;
	}

	//Convert from gps DD to UTM
	//TODO: For now, all is float cause I don't know what I'm doing and converting this code from JS, so come back and make this better
	//Original JS thanks to: www,movable-type.co.uk/scripts/latlong-utm-mgrs.html
	UTMData GpsDDtoUTM(double latDD, double lonDD) {
		//Assume lat/long are valid. If they're not we should already know because the accuracy value will be cracked

		//Establish easting and northings to ref
		double falseEasting = 500e3, falseNorthing = 10000e3;


		//Figure out zone & lat band
		uint8_t zone = floor((lonDD + 180) / 6) + 1; //Longitudinal zone;
		double lambda0 = ((zone - 1) * 6 - 180 + 3) * DEG2RAD; //Longitude of central meridian

		//Handle Norway/Svalbard exceptions
		char mgrsLatBand[] = "CDEFGHJKLMNPQRSTUVWXX"; //X is repeated for 80-84degN
		char latBand = mgrsLatBand[(int)floor(latDD / 8 + 10)];

		//adjust zone & central meridian for Norway
		if (zone == 31 && latBand == 'V' && lonDD >= 3) {
			zone++;
			lambda0 += 6 * DEG2RAD;
		}
		//adjust zone & central meridian for Svalbard
		if (zone == 32 && latBand == 'X' && lonDD < 9) { zone--; lambda0 -= 6 * DEG2RAD; }
		if (zone == 32 && latBand == 'X' && lonDD >= 9) { zone++; lambda0 += 6 * DEG2RAD; }

		if (zone == 34 && latBand == 'X' && lonDD < 21) { zone--; lambda0 -= 6 * DEG2RAD; }
		if (zone == 34 && latBand == 'X' && lonDD >= 21) { zone++; lambda0 += 6 * DEG2RAD; }

		if (zone == 36 && latBand == 'X' && lonDD < 33) { zone--; lambda0 -= 6 * DEG2RAD; }
		if (zone == 36 && latBand == 'X' && lonDD >= 33) { zone++; lambda0 += 6 * DEG2RAD; }

		double phi = latDD * DEG2RAD;
		double lambda = lonDD * DEG2RAD - lambda0;

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
		double sigma = sinh(e * atanh(e * tau / sqrt(1 + tau * tau)));

		double tauPrime = tau * sqrt(1 + sigma * sigma) - sigma * sqrt(1 + tau * tau);

		double zetaPrime = atan2(tauPrime, cosLambda);
		double etaPrime = asinh(sinLambda / sqrt(tauPrime * tauPrime + cosLambda * cosLambda));

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
		for (int j = 1; j <= 6; j++) {
			zeta += alpha[j] * sin(2 * j * zetaPrime) * cosh(2 * j * etaPrime);
		}

		double eta = etaPrime;
		for (int j = 1; j <= 6; j++) {
			eta += alpha[j] * cos(2 * j * zetaPrime) * sinh(2 * j * etaPrime);
		}

		double x = k0 * A * eta;
		double y = k0 * A * zeta;

		//Compute Convergance & Scale here, we're not going to cause we don't care...

			//Shift x/y to false origins
		x = x + falseEasting; //make x relative to false easting
		if (y < 0) {
			y = y + falseNorthing; //make y in southern hemisphere relative to false northing
		}

		char h = latDD > 0 ? 'N' : 'S';

		//return zone, latBand, h, x, y
		struct UTMData retData;
		retData.easting = x;
		retData.northing = y;
		retData.zone = zone;
		retData.latBand = latBand;
		retData.hemisphere = h;

		return retData;
	}

	DDCoord GpsUTMtoDD(UTMData utmCoord) {
		double falseEasting = 500e3, falseNorthing = 10000e3;

		//From WGS84 ellipsoid
		double a = 6378137, f = 1 / 293.257223563;

		//UTM scale on central meridian;
		double k0 = 0.9996;

		double x = utmCoord.easting - falseEasting;
		double y = utmCoord.hemisphere == 'S' ? utmCoord.northing - falseNorthing : utmCoord.northing;

		//From Karney 2011 Eq 15-22, 36
		double e = 0.08181919084; //SQRT(f*(2-f))
		double n = 0.00170789978; //f / (2-f)

		double A = 6367449.146;
		double eta = x / (6364902.166);
		double xi = y / (6364902.166);

		//6th order Kruger expressions
		double beta[7] = {
			0,
			8.52007199e-4,
			6.10987467e-8,
			1.76063649e-10,
			2.31646621e-13,
			4.12257831e-16,
			8.02400469e-19
		};

		double xi_prime = xi;
		for (int j = 1; j <= 6; j++)
			xi_prime -= beta[j] * sin(2 * j * xi) * cosh(2 * j * eta);


		double eta_prime = eta;
		for (int j = 1; j <= 6; j++)
			eta_prime -= beta[j] * cos(2 * j * xi) * sinh(2 * j * eta);


		double sinh_eta_prime = sinh(eta_prime);
		double sin_xi_prime = sin(xi_prime);
		double cos_xi_prime = cos(xi_prime);

		double tau_prime = sin_xi_prime / sqrt(sinh_eta_prime * sinh_eta_prime + cos_xi_prime * cos_xi_prime);

		double delta_tau_i = 0;
		double tau_i = tau_prime;
		
		do {
			double sigma_i = sinh(e * atanh(e * tau_i / sqrt(1 + tau_i * tau_i)));
			double tau_i_prime = tau_i * sqrt(1 + sigma_i * sigma_i) - sigma_i * sqrt(1 + tau_i * tau_i);
			delta_tau_i = (tau_prime - tau_i_prime) / sqrt(1 + tau_i_prime * tau_i_prime) *
				(1 + (1 - e * e) * tau_i * tau_i) / ((1 - e * e) * sqrt(1 + tau_i * tau_i));
			tau_i += delta_tau_i;
		} while (fabs(delta_tau_i) > 1e-12);

		double tau = tau_i;

		double phi = atan(tau);

		double lambda = atan2(sinh_eta_prime, cos_xi_prime);

		//Convergence
		double p = 1;
		for (int j = 1; j <= 6; j++)
			p -= 2 * j * beta[j] * cos(2 * j * xi) * cosh(2 * j * eta);
		double q = 0;
		for (int j = 1; j <= 6; j++)
			q += 2 * j * beta[j] * sin(2 * j * xi) * sinh(2 * j * eta);


		double gamma_prime = atan(tan(xi_prime) * tanh(eta_prime));
		double gamma_double_prime = atan2(q, p);
		double gamma = gamma_prime + gamma_double_prime;

		double sin_phi = sin(phi);
		double k_prime = sqrt(1 - e * e * sin_phi * sin_phi) * sqrt(1 + tau * tau) *
			sqrt(sinh_eta_prime * sinh_eta_prime + cos_xi_prime * cos_xi_prime);
		double k_double_prime = A / a / sqrt(p * p + q * q);
		double k = k0 * k_prime * k_double_prime;

		double lambda0 = ((utmCoord.zone - 1) * 6 - 180 + 3) * DEG2RAD;
		lambda += lambda0;

		double lat = phi * RAD2DEG;
		double lon = lambda * RAD2DEG;
		double convergence = gamma * RAD2DEG;
		double scale = k;

		//round...? Nah.

		DDCoord retCoord;
		retCoord.latitude = lat;
		retCoord.longitude = lon;

		//DC about convergence or scale.
		//If somebody reading this is...
		//Convergence == Y*RAD2DEG
		//Scale  = k

		return retCoord;

	}
}
#endif /* SRC_GPS_H_ */