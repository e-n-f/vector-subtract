#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MAX_ZOOM 28
#define ZOOM_BITS 5

#define FOOT .00000274
#define BUFFER (100 * FOOT)

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
void latlon2tile(float lat, float lon, int zoom, unsigned int *x, unsigned int *y) {
	float lat_rad = lat * M_PI / 180;
	unsigned long long n = 1LL << zoom;

	*x = n * ((lon + 180) / 360);
	*y = n * (1 - (log(tan(lat_rad) + 1/cos(lat_rad)) / M_PI)) / 2;
}

// http://wiki.openstreetmap.org/wiki/Slippy_map_tilenames
void tile2latlon(unsigned int x, unsigned int y, int zoom, float *lat, float *lon) {
	unsigned long long n = 1LL << zoom;
	*lon = 360.0 * x / n - 180.0;
	float lat_rad = atan(sinh(M_PI * (1 - 2.0 * y / n)));
	*lat = lat_rad * 180 / M_PI;
}

struct point {
	unsigned long long index;
	float minlat, minlon;
	float maxlat, maxlon;

	int n;
	float *lats, *lons;
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

struct index *index_init() {
	struct index *ix = malloc(sizeof(struct index));

	ix->points = NULL;
	ix->npoints = 0;
	ix->npalloc = 0;

	return ix;
}

void index_destroy(struct index *ix) {
	free(ix->points);
	free(ix);
}

void index_add(struct index *i, float minlat, float minlon, float maxlat, float maxlon, int n, float *lats, float *lons) {
	unsigned int x1, y1, x2, y2;

	if (minlat > maxlat) {
		float swap = minlat;
		minlat = maxlat;
		maxlat = swap;
	}
	if (minlon > maxlon) {
		float swap = minlon;
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

	i->points[i->npoints].n = n;
	i->points[i->npoints].lats = malloc(n * sizeof(float));
	i->points[i->npoints].lons = malloc(n * sizeof(float));
	memcpy(i->points[i->npoints].lats, lats, n * sizeof(float));
	memcpy(i->points[i->npoints].lons, lons, n * sizeof(float));

	i->npoints++;
}

void index_sort(struct index *ix) {
	qsort(ix->points, ix->npoints, sizeof(ix->points[0]), pointcmp);
}

void index_lookup(struct index *ix, float minlat, float minlon, float maxlat, float maxlon, int (*callback)(struct point *, void *), void *data) {
	unsigned x1, y1, x2, y2;
	int z;
	unsigned x, y;

	latlon2tile(minlat, minlon, 32, &x1, &y1);
	latlon2tile(maxlat, maxlon, 32, &x2, &y2);
	get_bbox_tile(x1, y1, x2, y2, &z, &x, &y);

	int zz;
	for (zz = MAX_ZOOM; zz >= 0; zz--) {
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

			if (callback(j, data)) {
				return;
			}
		}

		// printf("\t%016llx  %d\n", end, zz);
	}
}

struct seg {
	float lat1;
	float lon1;

	float lat2;
	float lon2;

	struct seg *next;
};

// http://www.ecse.rpi.edu/~wrf/Research/Short_Notes/pnpoly.html
int pnpoly(int nvert, float *vertx, float *verty, float testx, float testy)
{
  int i, j, c = 0;
  for (i = 0, j = nvert-1; i < nvert; j = i++) {
    if ( ((verty[i]>testy) != (verty[j]>testy)) &&
	 (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
       c = !c;
  }
  return c;
}

// http://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect
int intersect(float p0_x, float p0_y, float p1_x, float p1_y, 
    float p2_x, float p2_y, float p3_x, float p3_y, float *i_x, float *i_y)
{
    float s1_x, s1_y, s2_x, s2_y;
    s1_x = p1_x - p0_x;     s1_y = p1_y - p0_y;
    s2_x = p3_x - p2_x;     s2_y = p3_y - p2_y;

    float s, t;
    s = (-s1_y * (p0_x - p2_x) + s1_x * (p0_y - p2_y)) / (-s2_x * s1_y + s1_x * s2_y);
    t = ( s2_x * (p0_y - p2_y) - s2_y * (p0_x - p2_x)) / (-s2_x * s1_y + s1_x * s2_y);

    if (s >= 0 && s <= 1 && t >= 0 && t <= 1)
    {
        // Collision detected
        if (i_x != NULL)
            *i_x = p0_x + (t * s1_x);
        if (i_y != NULL)
            *i_y = p0_y + (t * s1_y);
        return 1;
    }

    return 0; // No collision
}

int callback(struct point *p, void *v) {
	struct seg **s = v;

	while (*s != NULL) {
		int p1 = pnpoly(p->n, p->lats, p->lons, (*s)->lat1, (*s)->lon1);
		int p2 = pnpoly(p->n, p->lats, p->lons, (*s)->lat2, (*s)->lon2);

		if (p1 && p2) {
			struct seg *cur = *s;
			struct seg *next = (*s)->next;
			*s = next;
			free(cur);
			continue;
		}

		int intersects[p->n];
		float intersect_lat[p->n];
		float intersect_lon[p->n];
		int nintersect = 0;

		int i;
		for (i = 0; i < p->n; i++) {
			intersects[i] = intersect(p->lats[i], p->lons[i], p->lats[(i + 1) % p->n], p->lons[(i + 1) % p->n], (*s)->lat1, (*s)->lon1, (*s)->lat2, (*s)->lon2, &intersect_lat[nintersect], &intersect_lon[nintersect]);

			if (intersects[i]) {
				nintersect++;
			}
		}

		if (p1 + p2 == 0 && (nintersect != 2 && nintersect != 0)) {
			fprintf(stderr, "0 within should intersect 0 or 2, not %d\n", nintersect);
			break;
		}
		if (p1 + p2 == 1 && nintersect != 1) {
			fprintf(stderr, "1 within should intersect 1, not %d\n", nintersect);
			break;
		}

		if (p1 || p2) {
			if (p1) {
				(*s)->lat1 = intersect_lat[0];
				(*s)->lon1 = intersect_lon[0];
			} else {
				(*s)->lat2 = intersect_lat[0];
				(*s)->lon2 = intersect_lon[0];
			}
		} else if (nintersect == 2) {
			float rat = cos((*s)->lat1 * M_PI / 180);
			float latd1 = (*s)->lat1 - intersect_lat[0];
			float lond1 = ((*s)->lon1 - intersect_lon[0]) * rat;
			float d1 = sqrt(latd1 * latd1 + lond1 * lond1);

			float latd2 = (*s)->lat1 - intersect_lat[1];
			float lond2 = ((*s)->lon1 - intersect_lon[1]) * rat;
			float d2 = sqrt(latd2 * latd2 + lond2 * lond2);

			struct seg *n = malloc(sizeof(struct seg));
			n->next = (*s)->next;
			n->lat2 = (*s)->lat2;
			n->lon2 = (*s)->lon2;
			(*s)->next = n;

			if (d1 < d2) {
				(*s)->lat2 = intersect_lat[0];
				(*s)->lon2 = intersect_lon[0];
				n->lat1 = intersect_lat[1];
				n->lon1 = intersect_lon[1];
			} else {
				(*s)->lat2 = intersect_lat[1];
				(*s)->lon2 = intersect_lon[1];
				n->lat1 = intersect_lat[0];
				n->lon1 = intersect_lon[0];
			}

			// don't recheck the split second half, which will
			// often mismatch by having one end right on the line
			s = &((*s)->next);
		}

		s = &((*s)->next);
	}

	s = v;
	if (*s == NULL) {
		return 1;
	} else {
		return 0;
	}
}

int main(int argc, char **argv) {
	char s[2000];
	long long seq = 0;

	struct index *ix = index_init();

	while (fgets(s, 2000, stdin)) {
		float lat1, lon1, lat2, lon2;

		if (seq++ % 100000 == 0) {
			fprintf(stderr, "%.1f million\r", seq / 1000000.0);
		}

		if (strcmp(s, "--\n") == 0) {
			break;
		}

		if (sscanf(s, "%f,%f %f,%f", &lat1, &lon1, &lat2, &lon2) == 4) {
			float rat = cos(lat1 * M_PI / 180);
			float ang = atan2(lat2 - lat1, (lon2 - lon1) * rat);

			float lats[] = {
				lat2 + BUFFER * sin(ang + M_PI / 4),
				lat2 + BUFFER * sin(ang + M_PI * 7 / 4),
				lat1 + BUFFER * sin(ang + M_PI * 5 / 4),
				lat1 + BUFFER * sin(ang + M_PI * 3 / 4),
			};

			float lons[] = {
				lon2 + BUFFER * cos(ang + M_PI / 4) / rat,
				lon2 + BUFFER * cos(ang + M_PI * 7 / 4) / rat,
				lon1 + BUFFER * cos(ang + M_PI * 5 / 4) / rat,
				lon1 + BUFFER * cos(ang + M_PI * 3 / 4) / rat,
			};

			float minlat = 360, minlon = 360, maxlat = -360, maxlon = -360;

			int i;
			for (i = 0; i < sizeof(lats) / sizeof(lats[0]); i++) {
				if (lats[i] < minlat) {
					minlat = lats[i];
				}
				if (lons[i] < minlon) {
					minlon = lons[i];
				}
				if (lats[i] > maxlat) {
					maxlat = lats[i];
				}
				if (lons[i] > maxlon) {
					maxlon = lons[i];
				}
			}

			index_add(ix, minlat, minlon, maxlat, maxlon, sizeof(lats) / sizeof(lats[0]), lats, lons);
		}
	}

	index_sort(ix);
	seq = 0;

	while (fgets(s, 2000, stdin)) {
		float lat1, lon1, lat2, lon2;

		if (seq++ % 10000 == 0) {
			fprintf(stderr, "checked %.3f million\r", seq / 1000000.0);
		}

		if (sscanf(s, "%f,%f %f,%f", &lat1, &lon1, &lat2, &lon2) == 4) {
			float minlat, minlon, maxlat, maxlon;

			if (lat1 < lat2) {
				minlat = lat1;
			} else {
				minlat = lat2;
			}

			if (lon1 < lon2) {
				minlon = lon1;
			} else {
				minlon = lon2;
			}

			if (lat1 > lat2) {
				maxlat = lat1;
			} else {
				maxlat = lat2;
			}

			if (lon1 > lon2) {
				maxlon = lon1;
			} else {
				maxlon = lon2;
			}

			struct seg *s = malloc(sizeof(struct seg));
			s->lat1 = lat1;
			s->lon1 = lon1;
			s->lat2 = lat2;
			s->lon2 = lon2;
			s->next = NULL;

			index_lookup(ix, minlat, minlon, maxlat, maxlon, callback, &s);

			struct seg *next;
			for (; s != NULL; s = next) {
				next = s->next;
				printf("%f,%f %f,%f\n", s->lat1, s->lon1, s->lat2, s->lon2);
				free(s);
			}
		}
	}

	index_destroy(ix);
	return 0;
}
