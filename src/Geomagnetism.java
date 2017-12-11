/*	License Statement from the NOAA
* The WMM source code is in the public domain and not licensed or
* under copyright. The information and software may be used freely
* by the public. As required by 17 U.S.C. 403, third parties producing
* copyrighted works consisting predominantly of the material produced
* by U.S. government agencies must provide notice with such work(s)
* identifying the U.S. Government material incorporated and stating
* that such material is not subject to copyright protection.*/

/** <p>Class to calculate magnetic declination, magnetic field strength,
* inclination etc. for any point on the earth.</p>
* <p>Adapted from the geomagc software and World Magnetic Model of the NOAA
* Satellite and Information Service, National Geophysical Data Center</p>
* http://www.ngdc.noaa.gov/geomag/WMM/DoDWMM.shtml
* <p>© Deep Pradhan, 2017</p>*/
class Geomagnetism {

	/** Initialise the instance without calculations*/
	Geomagnetism() {
		// Initialize constants
		maxord = MAX_DEG;
		sp[0] = 0;
		cp[0] = snorm[0] = pp[0] = 1;
		dp[0][0] = 0;

		c[0][0] = 0;
		cd[0][0] = 0;

		epoch = Double.parseDouble(WMM_COF[0].trim().split("[\\s]+")[0]);

		String[] tokens;

		double gnm, hnm, dgnm, dhnm;
		for (int i = 1, m, n; i < WMM_COF.length; i++) {
			tokens = WMM_COF[i].trim().split("[\\s]+");
			n = Integer.parseInt(tokens[0]);
			m = Integer.parseInt(tokens[1]);
			gnm = Double.parseDouble(tokens[2]);
			hnm = Double.parseDouble(tokens[3]);
			dgnm = Double.parseDouble(tokens[4]);
			dhnm = Double.parseDouble(tokens[5]);
			if (m <= n) {
				c[m][n] = gnm;
				cd[m][n] = dgnm;
				if (m != 0) {
					c[n][m - 1] = hnm;
					cd[n][m - 1] = dhnm;
				}
			}			
		}
		// Convert schmidt normalized gauss coefficients to unnormalized
		snorm[0] = 1;
		double flnmj;
		for (int j, n = 1; n <= maxord; n++) {
			snorm[n] = snorm[n - 1] * (2 * n - 1) / n;
			j = 2;
			for (int m = 0, d1 = 1, d2 = (n - m + d1) / d1; d2 > 0; d2--, m += d1) {
				k[m][n] = (double) (((n - 1) * (n - 1)) - (m * m)) / (double) ((2 * n - 1) * (2 * n - 3));
				if (m > 0) {
					flnmj = ((n - m + 1) * j) / (double) (n + m);
					snorm[n + m * 13] = snorm[n + (m -1) * 13] * Math.sqrt(flnmj);
					j = 1;
					c[n][m - 1] = snorm[n + m * 13] * c[n][m - 1];
					cd[n][m - 1] = snorm[n + m * 13] * cd[n][m - 1];
				}
				c[m][n] = snorm[n + m * 13] * c[m][n];
				cd[m][n] = snorm[n + m * 13] * cd[m][n];
			}
			fn[n] = (n + 1);
			fm[n] = n;
		}
		k[1][1] = 0;
		fm[0] = 0;
		otime = oalt = olat = olon = -1000;
	}

	/** Initialise the instance and calculate for given location and time
	*	@param longitude	Longitude in decimal degrees
	*	@param latitude		Latitude in decimal degrees
	*	@param altitude		Altitude in metres
	*	@param yearFraction	Date of the calculation in decimal years*/
	Geomagnetism(double longitude, double latitude, double altitude, double yearFraction) {
		this();
		calculate(longitude, latitude, altitude, yearFraction);
	}

	/** Calculate for given location and time
	*	@param longitude	Longitude in decimal degrees
	*	@param latitude		Latitude in decimal degrees
	*	@param altitude		Altitude in metres
	*	@param yearFraction	Date of the calculation in decimal years*/
	void calculate(double longitude, double latitude, double altitude, double yearFraction) {
		double rlon = Math.toRadians(longitude),
				rlat = Math.toRadians(latitude),
				altitudeKm = altitude / 1000,
				dt = yearFraction - epoch,
				srlon = Math.sin(rlon),
				srlat = Math.sin(rlat),
				crlon = Math.cos(rlon),
				crlat = Math.cos(rlat),
				srlat2 = srlat * srlat,
				crlat2 = crlat * crlat,
				a2 = WGS84_A * WGS84_A,
				b2 = WGS84_B * WGS84_B,
				c2 = a2 - b2,
				a4 = a2 * a2,
				b4 = b2 * b2,
				c4 = a4 - b4;

		sp[1] = srlon;
		cp[1] = crlon;

		// Convert from geodetic coords. to spherical coords.
		if (altitudeKm != oalt || latitude != olat) {
			double q = Math.sqrt(a2 - c2 * srlat2),
					q1 = altitudeKm * q,
					q2 = ((q1 + a2) / (q1 + b2)) * ((q1 + a2) / (q1 + b2)),
					r2 = ((altitudeKm * altitudeKm) + 2 * q1 + (a4 - c4 * srlat2) / (q * q));
			ct = srlat / Math.sqrt(q2 * crlat2 + srlat2);
			st = Math.sqrt(1 - (ct * ct));
			r = Math.sqrt(r2);
			d = Math.sqrt(a2 * crlat2 + b2 * srlat2);
			ca = (altitudeKm + d) / r;
			sa = c2 * crlat * srlat / (r * d);
		}
		if (longitude != olon) {
			for (int m = 2; m <= maxord; m++) {
				sp[m] = sp[1] * cp[m - 1] + cp[1] * sp[m - 1];
				cp[m] = cp[1] * cp[m - 1] - sp[1] * sp[m - 1];
			}
		}
		double aor = IAU66_RADIUS / r,
				ar = aor * aor,
				br = 0, bt = 0, bp = 0, bpp = 0,
				par, parp, temp1, temp2;

		for (int n = 1; n <= maxord; n++) {
			ar = ar * aor;
			for (int m = 0, d3 = 1, d4 = (n + m + d3) / d3; d4 > 0; d4--, m += d3) {

				// Compute unnormalized associated legendre polynomials and derivatives via recursion relations
				if (altitudeKm != oalt || latitude != olat) {
					if (n == m) {
						snorm[n + m * 13] = st * snorm[n - 1 + (m - 1) * 13];				
						dp[m][n] = st * dp[m - 1][n - 1]+ ct * snorm[n - 1 + (m - 1) * 13];
					}
					if (n == 1 && m == 0) {
						snorm[n + m * 13] = ct * snorm[n - 1 + m * 13];
						dp[m][n] = ct * dp[m][n - 1] - st * snorm[n - 1 + m * 13];
					}
					if (n > 1 && n != m) {
						if (m > n - 2) 
							snorm[n - 2 + m * 13] = 0;
						if (m > n - 2)
							dp[m][n - 2] = 0;
						snorm[n + m * 13] = ct * snorm[n - 1 + m * 13] - k[m][n] * snorm[n - 2 + m * 13];
						dp[m][n] = ct * dp[m][n - 1] - st * snorm[n - 1 + m * 13] - k[m][n] * dp[m][n - 2];
					}
				}

				// Time adjust the gauss coefficients
				if (yearFraction != otime) {
					tc[m][n] = c[m][n] + dt * cd[m][n];

					if (m != 0)
						tc[n][m - 1] = c[n][m - 1]+ dt * cd[n][m - 1];
				}

				//Accumulate terms of the spherical harmonic expansions
				par = ar * snorm[ n + m * 13];
				if (m == 0) {
					temp1 = tc[m][n] * cp[m];
					temp2 = tc[m][n] * sp[m];
				} else {
					temp1 = tc[m][n] * cp[m] + tc[n][m - 1] * sp[m];
					temp2 = tc[m][n] * sp[m] - tc[n][m - 1] * cp[m];
				}

				bt = bt - ar * temp1 * dp[m][n];
				bp += (fm[m] * temp2 * par);
				br += (fn[n] * temp1 * par);

				// Special case: north/south geographic poles
				if (st == 0 && m == 1) {
					if (n == 1)
						pp[n] = pp[n - 1];
					else 
						pp[n] = ct * pp[n - 1] - k[m][n] * pp[n - 2];
					parp = ar * pp[n];
					bpp += (fm[m] * temp2 * parp);
				}
			}
		}

		if (st == 0)
			bp = bpp;
		else 
			bp /= st;

		// Rotate magnetic vector components from spherical to geodetic coordinates
		// bx must be the east-west field component
		// by must be the north-south field component
		// bz must be the vertical field component.
		bx = -bt * ca - br * sa;
		by = bp;
		bz = bt * sa - br * ca;

		// Compute declination (dec), inclination (dip) and total intensity (ti)
		bh = Math.sqrt((bx * bx)+(by * by));
		intensity = Math.sqrt((bh * bh)+(bz * bz));
		//	Calculate the declination.
		declination = Math.toDegrees(Math.atan2(by, bx));
		inclination = Math.toDegrees(Math.atan2(bz, bh));

		otime = yearFraction;
		oalt = altitudeKm;
		olat = latitude;
		olon = longitude;
	}

	/** @return The declination in degrees*/
	double getDeclination() {
		return declination;
	}

	/** @return The magnetic field intensity/strength in nano Tesla*/
	double getIntensity() {
		return intensity;
	}

	/** @return The horizontal magnetic field intensity/strength in nano Tesla*/
	double getHorizontalIntensity() {
		return bh;
	}

	/** @return The vertical magnetic field intensity/strength in nano Tesla*/
	double getVerticalIntensity() {
		return bz;
	}

	/** @return The northerly component of the magnetic field strength in nano Tesla*/
	double getNorthIntensity() {
		return bx;
	}

	/** @return The easterly component of the magnetic field strength in nano Tesla*/
	double getEastIntensity() {
		return by;
	}

	/** @return The magnetic field dip angle, in degrees*/
	double getDipAngle() {
		return inclination;
	}

	/**	The input string array which contains each line of input for the wmm.cof input file.
	*	The columns in this file are as follows:	n,	m,	gnm,	hnm,	dgnm,	dhnm*/
	private final static String [] WMM_COF = {"    2015.0            WMM-2015        12/15/2014",
			"  1  0  -29438.5       0.0       10.7        0.0",
			"  1  1   -1501.1    4796.2       17.9      -26.8",
			"  2  0   -2445.3       0.0       -8.6        0.0",
			"  2  1    3012.5   -2845.6       -3.3      -27.1",
			"  2  2    1676.6    -642.0        2.4      -13.3",
			"  3  0    1351.1       0.0        3.1        0.0",
			"  3  1   -2352.3    -115.3       -6.2        8.4",
			"  3  2    1225.6     245.0       -0.4       -0.4",
			"  3  3     581.9    -538.3      -10.4        2.3",
			"  4  0     907.2       0.0       -0.4        0.0",
			"  4  1     813.7     283.4        0.8       -0.6",
			"  4  2     120.3    -188.6       -9.2        5.3",
			"  4  3    -335.0     180.9        4.0        3.0",
			"  4  4      70.3    -329.5       -4.2       -5.3",
			"  5  0    -232.6       0.0       -0.2        0.0",
			"  5  1     360.1      47.4        0.1        0.4",
			"  5  2     192.4     196.9       -1.4        1.6",
			"  5  3    -141.0    -119.4        0.0       -1.1",
			"  5  4    -157.4      16.1        1.3        3.3",
			"  5  5       4.3     100.1        3.8        0.1",
			"  6  0      69.5       0.0       -0.5        0.0",
			"  6  1      67.4     -20.7       -0.2        0.0",
			"  6  2      72.8      33.2       -0.6       -2.2",
			"  6  3    -129.8      58.8        2.4       -0.7",
			"  6  4     -29.0     -66.5       -1.1        0.1",
			"  6  5      13.2       7.3        0.3        1.0",
			"  6  6     -70.9      62.5        1.5        1.3",
			"  7  0      81.6       0.0        0.2        0.0",
			"  7  1     -76.1     -54.1       -0.2        0.7",
			"  7  2      -6.8     -19.4       -0.4        0.5",
			"  7  3      51.9       5.6        1.3       -0.2",
			"  7  4      15.0      24.4        0.2       -0.1",
			"  7  5       9.3       3.3       -0.4       -0.7",
			"  7  6      -2.8     -27.5       -0.9        0.1",
			"  7  7       6.7      -2.3        0.3        0.1",
			"  8  0      24.0       0.0        0.0        0.0",
			"  8  1       8.6      10.2        0.1       -0.3",
			"  8  2     -16.9     -18.1       -0.5        0.3",
			"  8  3      -3.2      13.2        0.5        0.3",
			"  8  4     -20.6     -14.6       -0.2        0.6",
			"  8  5      13.3      16.2        0.4       -0.1",
			"  8  6      11.7       5.7        0.2       -0.2",
			"  8  7     -16.0      -9.1       -0.4        0.3",
			"  8  8      -2.0       2.2        0.3        0.0",
			"  9  0       5.4       0.0        0.0        0.0",
			"  9  1       8.8     -21.6       -0.1       -0.2",
			"  9  2       3.1      10.8       -0.1       -0.1",
			"  9  3      -3.1      11.7        0.4       -0.2",
			"  9  4       0.6      -6.8       -0.5        0.1",
			"  9  5     -13.3      -6.9       -0.2        0.1",
			"  9  6      -0.1       7.8        0.1        0.0",
			"  9  7       8.7       1.0        0.0       -0.2",
			"  9  8      -9.1      -3.9       -0.2        0.4",
			"  9  9     -10.5       8.5       -0.1        0.3",
			" 10  0      -1.9       0.0        0.0        0.0",
			" 10  1      -6.5       3.3        0.0        0.1",
			" 10  2       0.2      -0.3       -0.1       -0.1",
			" 10  3       0.6       4.6        0.3        0.0",
			" 10  4      -0.6       4.4       -0.1        0.0",
			" 10  5       1.7      -7.9       -0.1       -0.2",
			" 10  6      -0.7      -0.6       -0.1        0.1",
			" 10  7       2.1      -4.1        0.0       -0.1",
			" 10  8       2.3      -2.8       -0.2       -0.2",
			" 10  9      -1.8      -1.1       -0.1        0.1",
			" 10 10      -3.6      -8.7       -0.2       -0.1",
			" 11  0       3.1       0.0        0.0        0.0",
			" 11  1      -1.5      -0.1        0.0        0.0",
			" 11  2      -2.3       2.1       -0.1        0.1",
			" 11  3       2.1      -0.7        0.1        0.0",
			" 11  4      -0.9      -1.1        0.0        0.1",
			" 11  5       0.6       0.7        0.0        0.0",
			" 11  6      -0.7      -0.2        0.0        0.0",
			" 11  7       0.2      -2.1        0.0        0.1",
			" 11  8       1.7      -1.5        0.0        0.0",
			" 11  9      -0.2      -2.5        0.0       -0.1",
			" 11 10       0.4      -2.0       -0.1        0.0",
			" 11 11       3.5      -2.3       -0.1       -0.1",
			" 12  0      -2.0       0.0        0.1        0.0",
			" 12  1      -0.3      -1.0        0.0        0.0",
			" 12  2       0.4       0.5        0.0        0.0",
			" 12  3       1.3       1.8        0.1       -0.1",
			" 12  4      -0.9      -2.2       -0.1        0.0",
			" 12  5       0.9       0.3        0.0        0.0",
			" 12  6       0.1       0.7        0.1        0.0",
			" 12  7       0.5      -0.1        0.0        0.0",
			" 12  8      -0.4       0.3        0.0        0.0",
			" 12  9      -0.4       0.2        0.0        0.0",
			" 12 10       0.2      -0.9        0.0        0.0",
			" 12 11      -0.9      -0.2        0.0        0.0",
			" 12 12       0.0       0.7        0.0        0.0",};

	/** Mean radius of IAU-66 ellipsoid, in km.*/
	private static final double IAU66_RADIUS = 6371.2;

	/** Semi-major axis of WGS-84 ellipsoid, in km.*/
	private static final double WGS84_A = 6378.137;

	/** Semi-minor axis of WGS-84 ellipsoid, in km.*/
	private static final double WGS84_B = 6356.7523142;

	/** The maximum number of degrees of the spherical harmonic model.*/
	private static final int MAX_DEG = 12;

	/** Geomagnetic declination in deg.
	*	East is positive, West is negative.
	*	(The negative of variation.)*/
	private double declination = 0;

	/** Geomagnetic inclination in deg.
	*	Down is positive, up is negative.*/
	private double inclination = 0;

	/** Geomagnetic total intensity, in nano Teslas.*/
	private double intensity = 0;

	/** The maximum order of spherical harmonic model.*/
	private int maxord;

	/** The Gauss coefficients of main geomagnetic model (nt).*/
	private double c[][] = new double[13][13];

	/** The Gauss coefficients of secular geomagnetic model (nt/yr).*/
	private double cd[][] = new double[13][13];

	/** The time adjusted geomagnetic gauss coefficients (nt).*/
	private double tc[][] = new double[13][13];

	/** The theta derivative of p(n,m) (unnormalized).*/
	private double dp[][] = new double[13][13];

	/** The Schmidt normalization factors.*/
	private double snorm[] = new double[169];

	/** The sine of (m*spherical coord. longitude).*/
	private double sp[] = new double[13];

	/** The cosine of (m*spherical coord. longitude).*/
	private double cp[] = new double[13];
	private double fn[] = new double[13];
	private double fm[] = new double[13];

	/** The associated Legendre polynomials for m = 1 (unnormalized).*/
	private double pp[] = new double[13];

	private double k[][] = new double[13][13];

	/** The variables otime (old time), oalt (old altitude),
	*	olat (old latitude), olon (old longitude), are used to
	*	store the values used from the previous calculation to
	*	save on calculation time if some inputs don't change.*/
	private double otime, oalt, olat, olon;

	/** The date in years, for the start of the valid time of the fit coefficients*/
	private double epoch;

	/** bx is the north south field intensity
	*	by is the east west field intensity
	*	bz is the vertical field intensity positive downward
	*	bh is the horizontal field intensity*/
	private double bx, by, bz, bh;

	private double r, d, ca, sa, ct, st;
}
