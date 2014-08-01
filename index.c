#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_ZOOM 28
#define ZOOM_BITS 5

int get_bbox_zoom(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2) {
	int z;
	for (z = 0; z < MAX_ZOOM; z++) {
		int mask = 1 << (32 - (z + 1));

		if (((x1 & mask) != (x2 & mask)) ||
		    ((y1 & mask) != (y2 & mask))) {
			return z;
		}
	}

	return MAX_ZOOM;
}

void get_bbox_tile(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2, int *z, unsigned int *x, unsigned int *y) {
	*z = get_bbox_zoom(x1, y1, x2, y2);
	*x = x1 >> (32 - *z);
	*y = y1 >> (32 - *z);
}

/*
 *  5 bits for zoom             (<< 59)
 * 56 bits for interspersed yx  (<< 3)
 *  3 bits for tags             (<< 0)
 */
unsigned long long encode_bbox(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2, int tags) {
	int z = get_bbox_zoom(x1, y1, x2, y2);
	long long out = ((long long) z) << (64 - ZOOM_BITS);

	int i;
	for (i = 0; i < MAX_ZOOM; i++) {
		long long v = ((y1 >> (32 - (i + 1))) & 1) << 1;
		v |= (x1 >> (32 - (i + 1))) & 1;
		v = v << (64 - ZOOM_BITS - 2 * (i + 1));

		out |= v;
	}

	return out;
}

void encode_tile(int zz, int z, unsigned int x, unsigned int y, unsigned long long *start, unsigned long long *end) {
	long long out = ((long long) zz) << (64 - ZOOM_BITS);

	x <<= (32 - z);
	y <<= (32 - z);

	int i;
	for (i = 0; i < MAX_ZOOM; i++) {
		long long v = ((y >> (32 - (i + 1))) & 1) << 1;
		v |= (x >> (32 - (i + 1))) & 1;
		v = v << (64 - ZOOM_BITS - 2 * (i + 1));

		out |= v;
	}

	*start = out;
	*end = out | (((unsigned long long) -1LL) >> (2 * z + ZOOM_BITS));
}

void decode_bbox(unsigned long long code, int *z, unsigned int *wx, unsigned int *wy) {
	*z = code >> (64 - ZOOM_BITS);
	*wx = 0;
	*wy = 0;

	int i;
	for (i = 0; i < MAX_ZOOM; i++) {
		long long v = code >> (64 - ZOOM_BITS - 2 * (i + 1));

		*wy |= ((v & 2) >> 1) << (32 - (i + 1));
		*wx |= (v & 1) << (32 - (i + 1));
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

struct index {
	struct point *points;
	int npoints;
	int npalloc;
};

void index_add(struct index *i, double minlat, double minlon, double maxlat, double maxlon) {
	unsigned int x1, y1, x2, y2;

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
	unsigned long long enc = encode_bbox(x1, y1, x2, y2, 0);

	if (i->npoints + 1 >= i->npalloc) {
		i->npalloc = (i->npalloc + 1000) * 3 / 2;
		i->points = realloc(i->points, i->npalloc * sizeof(i->points[0]));
	}

	i->points[i->npoints].index = enc;
	i->points[i->npoints].minlat = minlat;
	i->points[i->npoints].minlon = minlon;
	i->points[i->npoints].maxlat = maxlat;
	i->points[i->npoints].maxlon = maxlon;
	i->npoints++;
}

void index_lookup(struct index *ix, double minlat, double minlon, double maxlat, double maxlon, void (*callback)(struct point *, void *), void *data) {
	unsigned x1, y1, x2, y2;
	int z;
	unsigned x, y;

	latlon2tile(minlat, minlon, 32, &x1, &y1);
	latlon2tile(maxlat, maxlon, 32, &x2, &y2);
	get_bbox_tile(x1, y1, x2, y2, &z, &x, &y);

	int zz;
	for (zz = 0; zz <= MAX_ZOOM; zz++) {
		unsigned long long start, end;

		if (zz < z) {
			encode_tile(zz, zz, x >> (z - zz), y >> (z - zz), &start, &end);
		} else {
			encode_tile(zz, z, x, y, &start, &end);
		}

		struct point *pstart = search(&start, ix->points, ix->npoints, sizeof(ix->points[0]), pointcmp);
		struct point *pend = search(&end, ix->points, ix->npoints, sizeof(ix->points[0]), pointcmp);

		if (pend >= ix->points + ix->npoints) {
			pend = ix->points + ix->npoints - 1;
		}
		while (pstart > ix->points && pointcmp(pstart - 1, &start) == 0) {
			pstart--;
		}
		if (pointcmp(pstart, &start) < 0) {
			pstart++;
		}
		if (pointcmp(pend, &end) > 0) {
			pend--;
		}

		struct point *j;
		for (j = pstart; j <= pend; j++) {
			// reject by bbox
			if (j->minlat > maxlat ||
			    j->minlon > maxlon ||
			    minlat > j->maxlat ||
			    minlon > j->maxlon) {
				continue;
			}

			callback(j, data);
		}

		// printf("\t%016llx  %d\n", end, zz);
	}
}

void callback(struct point *p, void *v) {
	int *i = v;
	(*i)++;
}

int main(int argc, char **argv) {
	char s[2000];

	struct index ix;
	ix.points = NULL;
	ix.npoints = 0;
	ix.npalloc = 0;

	while (fgets(s, 2000, stdin)) {
		double minlat, minlon, maxlat, maxlon;

		if (sscanf(s, "%lf,%lf %lf,%lf", &minlat, &minlon, &maxlat, &maxlon) == 4) {
			index_add(&ix, minlat, minlon, maxlat, maxlon);
		}
	}

	qsort(ix.points, ix.npoints, sizeof(ix.points[0]), pointcmp);

	int i;
	for (i = 0; i < ix.npoints; i++) {
		int possible = 0;
		int matchedself = 0;

		index_lookup(&ix, ix.points[i].minlat, ix.points[i].minlon, ix.points[i].maxlat, ix.points[i].maxlon, callback, &possible);

		printf("%d\n", possible);

#if 0
		if (!matchedself) {
			fprintf(stderr, "did not match self: %llx %d    %f,%f %f,%f\n", ix.points[i].index, z, ix.points[i].minlat, ix.points[i].minlon, ix.points[i].maxlat, ix.points[i].maxlon);
		}
#endif
	}
}
