#ifndef SRC_GPS_H_
#define SRC_GPS_H_


#define PI    3.14159265358979323e0  /* PI */
#define DEG2RAD	PI/180
#define RAD2DEG 180/PI

typedef struct UTMData{
	double easting;
	double northing;
	uint8_t zone;
	char latBand;
	char hemisphere;
} UTMData;

typedef struct DMSData{
	float latDegree;
	float latMinute;
	float latSecond;
	float lonDegree;
	float lonMinute;
	float lonSecond;
} DMSData;

//Convert the Gps coords to minutes
DMSData GpsDDToDMS(double latDD, double lonDD);

//Convert DD to UTM
UTMData GpsDDtoUTM(double latDD, double lonDD);

#endif /* SRC_GPS_H_ */
