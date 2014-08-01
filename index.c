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
unsigned long long encode_bbox(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2, int tags) {
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

void encode_tile(int zz, int z, unsigned int x, unsigned int y, unsigned long long *start, unsigned long long *end) {
	long long out = ((long long) zz) << 59;

	x <<= (32 - z);
	y <<= (32 - z);

	int i;
	for (i = 0; i < 28; i++) {
		long long v = ((y >> (31 - i)) & 1) << 1;
		v |= (x >> (31 - i)) & 1;
		v = v << (57 - 2 * i);

		out |= v;
	}

	*start = out;
	*end = out | (((unsigned long long) -1LL) >> (2 * z + 5));
}

void decode_bbox(unsigned long long code, int *z, unsigned int *wx, unsigned int *wy) {
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

struct point {
	unsigned long long index;
	double minlat, minlon;
	double maxlat, maxlon;
};

int pointcmp(const void *v1, const void *v2) {
	const unsigned long long *p1 = v1;
	const unsigned long long *p2 = v2;

	if (*p1 < *p2) {
		return -1;
	} else if (*p1 > *p2) {
		return 1;
	} else {
		return 0;
	}
}

// http://www.tbray.org/ongoing/When/200x/2003/03/22/Binary
void *search(const void *key, const void *base, size_t nel, size_t width,
		int (*cmp)(const void *, const void *)) {

	long long high = nel, low = -1, probe;
	while (high - low > 1) {
		probe = (low + high) >> 1;
		int c = cmp(((char *) base) + probe * width, key);
		if (c > 0) {
			high = probe;
		} else {
			low = probe;
		}
	}

	if (low < 0) {
		low = 0;
	}

	return ((char *) base) + low * width;
}

int main(int argc, char **argv) {
	char s[2000];

	struct point *points = NULL;
	int npalloc = 0;
	int npoints = 0;

	while (fgets(s, 2000, stdin)) {
		double minlat, minlon, maxlat, maxlon;

		if (sscanf(s, "%lf,%lf %lf,%lf", &minlat, &minlon, &maxlat, &maxlon) == 4) {
			unsigned int x1, y1, x2, y2;

			// act like these are bboxes
			if (minlat > maxlat) {
				double swap = minlat;
				minlat = maxlat;
				maxlat = swap;
			}
			if (minlon > maxlon) {
				double swap = minlon;
				minlon = maxlon;
				maxlon = swap;
			}

			latlon2tile(minlat, minlon, 32, &x1, &y1);
			latlon2tile(maxlat, maxlon, 32, &x2, &y2);
			//int z = get_bbox_zoom(x1, y1, x2, y2);
			unsigned long long enc = encode_bbox(x1, y1, x2, y2, 0);

			if (npoints + 1 >= npalloc) {
				npalloc = (npalloc + 1000) * 3 / 2;
				points = realloc(points, npalloc * sizeof(points[0]));
			}

			points[npoints].index = enc;
			points[npoints].minlat = minlat;
			points[npoints].minlon = minlon;
			points[npoints].maxlat = maxlat;
			points[npoints].maxlon = maxlon;
			npoints++;
		}
	}

	qsort(points, npoints, sizeof(points[0]), pointcmp);

	int i;
	for (i = 0; i < npoints; i++) {
		unsigned wx, wy;
		int z;
		double minlat, minlon;

		decode_bbox(points[i].index, &z, &wx, &wy);
		tile2latlon(wx, wy, 32, &minlat, &minlon);

		// printf("%llx %d %f,%f    %f,%f %f,%f\n", points[i].index, z, minlat, minlon, points[i].minlat, points[i].minlon, points[i].maxlat, points[i].maxlon);

		unsigned x = wx >> (32 - z);
		unsigned y = wy >> (32 - z);

		// printf("\t%d/%u/%u\n", z, x, y);

		int possible = 0;
		int matchedself = 0;

		int zz;
		for (zz = 0; zz <= 28; zz++) {
			unsigned long long start, end;

			if (zz < z) {
				encode_tile(zz, zz, x >> (z - zz), y >> (z - zz), &start, &end);
			} else {
				encode_tile(zz, z, x, y, &start, &end);
			}

			// printf("\t%016llx  %d\n", start, zz);

			struct point *pstart = search(&start, points, npoints, sizeof(points[0]), pointcmp);
			struct point *pend = search(&end, points, npoints, sizeof(points[0]), pointcmp);

			if (pend >= points + npoints) {
				pend = points + npoints - 1;
			}

			if (pointcmp(pstart, &start) < 0) {
				pstart++;
			}
			if (pointcmp(pend, &end) > 0) {
				pend--;
			}

			struct point *j;
			for (j = pstart; j <= pend; j++) {
				int dz;
				unsigned int dwx, dwy;

				decode_bbox(j->index, &dz, &dwx, &dwy);

#if 0
				// reject by bbox
				if (j->minlat > points[i].maxlat ||
				    j->minlon > points[i].maxlon ||
				    points[i].minlat > j->maxlat ||
				    points[i].minlon > j->maxlon) {
					continue;
				}
#endif

				// printf("\t%llx  %d  %f,%f %f,%f\n", j->index, dz, j->minlat, j->minlon, j->maxlat, j->maxlon);
				possible++;

				if (j == points + i) {
					matchedself = 1;
				}
			}

			// printf("\t%016llx  %d\n", end, zz);
		}

		printf("%d\n", possible);

		if (!matchedself) {
			fprintf(stderr, "did not match self: %llx %d %f,%f    %f,%f %f,%f\n", points[i].index, z, minlat, minlon, points[i].minlat, points[i].minlon, points[i].maxlat, points[i].maxlon);
		}
	}
}
