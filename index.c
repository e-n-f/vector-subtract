#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int get_bbox_zoom(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2) {
	int z;
	for (z = 0; z < 28; z++) {
		int mask = 1 << (31 - z);

		if (((x1 & mask) != (x2 & mask)) ||
		    ((y1 & mask) != (y2 & mask))) {
			return z;
		}
	}

	return 28;
}

/*
 *  5 bits for zoom             (<< 59)
 * 56 bits for interspersed yx  (<< 3)
 *  3 bits for tags             (<< 0)
 */
long long encode_bbox(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2, int tags) {
	int z = get_bbox_zoom(x1, y1, x2, y2);
	long long out = ((long long) z) << 59;

	int i;
	for (i = 0; i < 28; i++) {
		long long v = ((y1 >> (31 - i)) & 1) << 1;
		v |= (x1 >> (31 - i)) & 1;
		v = v << (57 - 2 * i);

		out |= v;
	}

	return out;
}

void decode_bbox(long long code, int *z, unsigned int *wx, unsigned int *wy) {
	*z = code >> 59;
	*wx = 0;
	*wy = 0;

	int i;
	for (i = 0; i < 28; i++) {
		long long v = code >> (57 - 2 * i);

		*wy |= ((v & 2) >> 1) << (31 - i);
		*wx |= (v & 1) << (31 - i);
	}
}

// http://wiki.openstreetmap.org/wiki/Slippy_map_tilenames
void latlon2tile(double lat, double lon, int zoom, unsigned int *x, unsigned int *y) {
	double lat_rad = lat * M_PI / 180;
	unsigned long long n = 1LL << zoom;

	*x = n * ((lon + 180) / 360);
	*y = n * (1 - (log(tan(lat_rad) + 1/cos(lat_rad)) / M_PI)) / 2;
}

// http://wiki.openstreetmap.org/wiki/Slippy_map_tilenames
void tile2latlon(unsigned int x, unsigned int y, int zoom, double *lat, double *lon) {
	unsigned long long n = 1LL << zoom;
	*lon = 360.0 * x / n - 180.0;
	double lat_rad = atan(sinh(M_PI * (1 - 2.0 * y / n)));
	*lat = lat_rad * 180 / M_PI;
}

int main(int argc, char **argv) {
	char s[2000];

	while (fgets(s, 2000, stdin)) {
		double lat1, lon1, lat2, lon2;

		if (sscanf(s, "%lf,%lf %lf,%lf", &lat1, &lon1, &lat2, &lon2) == 4) {
			unsigned int x1, y1, x2, y2;
			latlon2tile(lat1, lon1, 32, &x1, &y1);
			latlon2tile(lat2, lon2, 32, &x2, &y2);
			int z = get_bbox_zoom(x1, y1, x2, y2);
			long long enc = encode_bbox(x1, y1, x2, y2, 0);

			printf("%f,%f %f,%f %x %x %x %x: %d %llx ", lat1, lon1, lat2, lon2,
				x1, y1, x2, y2, z, enc);

			decode_bbox(enc, &z, &x1, &y1);
			tile2latlon(x1, y1, 32, &lat1, &lon1);

			printf("%x %x %f %f\n", x1, y1, lat1, lon1);
		}
	}
}
