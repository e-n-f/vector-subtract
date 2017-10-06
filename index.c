#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "jsonpull/jsonpull.h"

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

	if (*z == 0) {
		*x = *y = 0;
	} else {
		*x = x1 >> (32 - *z);
		*y = y1 >> (32 - *z);
	}
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
	*y = n * (1 - (log(tan(lat_rad) + 1 / cos(lat_rad)) / M_PI)) / 2;
}

// http://wiki.openstreetmap.org/wiki/Slippy_map_tilenames
void tile2latlon(unsigned int x, unsigned int y, int zoom, float *lat, float *lon) {
	unsigned long long n = 1LL << zoom;
	*lon = 360.0 * x / n - 180.0;
	float lat_rad = atan(sinh(M_PI * (1 - 2.0 * y / n)));
	*lat = lat_rad * 180 / M_PI;
}

struct pointdata {
	float minlat, minlon;
	float maxlat, maxlon;

	int n;
	float latlon[0];
};

struct point {
	unsigned long long index;
	struct pointdata *pointdata;
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

	i->points[i->npoints].pointdata = malloc(sizeof(struct pointdata) + 2 * n * sizeof(float));
	i->points[i->npoints].pointdata->minlat = minlat;
	i->points[i->npoints].pointdata->minlon = minlon;
	i->points[i->npoints].pointdata->maxlat = maxlat;
	i->points[i->npoints].pointdata->maxlon = maxlon;
	i->points[i->npoints].pointdata->n = n;

	memcpy(i->points[i->npoints].pointdata->latlon, lats, n * sizeof(float));
	memcpy(i->points[i->npoints].pointdata->latlon + n, lons, n * sizeof(float));

	i->npoints++;
}

void index_sort(struct index *ix) {
	qsort(ix->points, ix->npoints, sizeof(ix->points[0]), pointcmp);
}

int range_lookup(struct index *ix, float minlat, float minlon, float maxlat, float maxlon, int (*callback)(struct point *, void *), void *data, long long start, long long end) {
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
		if (j->pointdata->minlat > maxlat ||
		    j->pointdata->minlon > maxlon ||
		    minlat > j->pointdata->maxlat ||
		    minlon > j->pointdata->maxlon) {
			continue;
		}

		if (callback(j, data)) {
			return 1;
		}
	}

	return 0;
}

void index_lookup(struct index *ix, float minlat, float minlon, float maxlat, float maxlon, int (*callback)(struct point *, void *), void *data) {
	unsigned x1, y1, x2, y2;
	int z;
	unsigned x, y;

	latlon2tile(minlat, minlon, 32, &x1, &y1);
	latlon2tile(maxlat, maxlon, 32, &x2, &y2);
	get_bbox_tile(x1, y1, x2, y2, &z, &x, &y);

	int z1 = z;
	int z2 = z + 1;

	for (; z1 >= 0 || z2 <= MAX_ZOOM; z1--, z2++) {
		unsigned long long start, end;

		if (z1 >= 0) {
			encode_tile(z1, z1, x >> (z - z1), y >> (z - z1), &start, &end);
			if (range_lookup(ix, minlat, minlon, maxlat, maxlon, callback, data, start, end)) {
				return;
			}
		}
		if (z2 <= MAX_ZOOM) {
			encode_tile(z2, z, x, y, &start, &end);
			if (range_lookup(ix, minlat, minlon, maxlat, maxlon, callback, data, start, end)) {
				return;
			}
		}
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
int pnpoly(int nvert, float *vertx, float *verty, float testx, float testy) {
	int i, j, c = 0;
	for (i = 0, j = nvert - 1; i < nvert; j = i++) {
		if (((verty[i] > testy) != (verty[j] > testy)) &&
		    (testx < (vertx[j] - vertx[i]) * (testy - verty[i]) / (verty[j] - verty[i]) + vertx[i]))
			c = !c;
	}
	return c;
}

// http://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect
int intersect(float p0_x, float p0_y, float p1_x, float p1_y,
	      float p2_x, float p2_y, float p3_x, float p3_y, float *i_x, float *i_y) {
	float s1_x, s1_y, s2_x, s2_y;
	s1_x = p1_x - p0_x;
	s1_y = p1_y - p0_y;
	s2_x = p3_x - p2_x;
	s2_y = p3_y - p2_y;

	float s, t;
	s = (-s1_y * (p0_x - p2_x) + s1_x * (p0_y - p2_y)) / (-s2_x * s1_y + s1_x * s2_y);
	t = (s2_x * (p0_y - p2_y) - s2_y * (p0_x - p2_x)) / (-s2_x * s1_y + s1_x * s2_y);

	if (s >= 0 && s <= 1 && t >= 0 && t <= 1) {
		// Collision detected
		if (i_x != NULL)
			*i_x = p0_x + (t * s1_x);
		if (i_y != NULL)
			*i_y = p0_y + (t * s1_y);
		return 1;
	}

	return 0;  // No collision
}

int callback(struct point *p, void *v) {
	struct seg **s = v;

	while (*s != NULL) {
		int p1 = pnpoly(p->pointdata->n, p->pointdata->latlon, p->pointdata->latlon + p->pointdata->n, (*s)->lat1, (*s)->lon1);
		int p2 = pnpoly(p->pointdata->n, p->pointdata->latlon, p->pointdata->latlon + p->pointdata->n, (*s)->lat2, (*s)->lon2);

		if (p1 && p2) {
			struct seg *cur = *s;
			struct seg *next = (*s)->next;
			*s = next;
			free(cur);
			continue;
		}

		int intersects[p->pointdata->n];
		float intersect_lat[p->pointdata->n];
		float intersect_lon[p->pointdata->n];
		int nintersect = 0;

		int n = p->pointdata->n;
		float *lats = p->pointdata->latlon;
		float *lons = p->pointdata->latlon + n;

		int i;
		for (i = 0; i < n; i++) {
			intersects[i] = intersect(lats[i], lons[i], lats[(i + 1) % n], lons[(i + 1) % n], (*s)->lat1, (*s)->lon1, (*s)->lat2, (*s)->lon2, &intersect_lat[nintersect], &intersect_lon[nintersect]);

			if (intersects[i]) {
				nintersect++;
			}
		}

		if (p1 + p2 == 0 && (nintersect != 2 && nintersect != 0)) {
			// fprintf(stderr, "0 within should intersect 0 or 2, not %d\n", nintersect);
			s = &((*s)->next);
			continue;
		}
		if (p1 + p2 == 1 && nintersect != 1) {
			// fprintf(stderr, "1 within should intersect 1, not %d\n", nintersect);
			s = &((*s)->next);
			continue;
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

void index_add1(struct index *ix, double lat1, double lon1, double lat2, double lon2) {
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

void compare(struct index *ix, double lat1, double lon1, double lat2, double lon2, const char *props) {
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

	struct seg *seg = malloc(sizeof(struct seg));
	seg->lat1 = lat1;
	seg->lon1 = lon1;
	seg->lat2 = lat2;
	seg->lon2 = lon2;
	seg->next = NULL;

	index_lookup(ix, minlat, minlon, maxlat, maxlon, callback, &seg);

	struct seg *next;
	for (; seg != NULL; seg = next) {
		next = seg->next;
		printf("{\"type\":\"Feature\",\"properties\":%s,\"geometry\":{\"type\":\"LineString\",\"coordinates\":[[%f,%f],[%f,%f]]}}\n",
		       props, seg->lon1, seg->lat1, seg->lon2, seg->lat2);
		free(seg);
	}
}

int main(int argc, char **argv) {
	long long seq = 0;

	struct index *ix = index_init();

	json_pull *jp = json_begin_file(stdin);

	json_object *j;
	while ((j = json_read(jp)) != NULL) {
		if (j->type == JSON_STRING && strcmp(j->string, "end part 1") == 0) {
			break;
		}

		if (j->type == JSON_HASH) {
			json_object *type = json_hash_get(j, "type");

			if (type != NULL && type->type == JSON_STRING && strcmp(type->string, "Feature") == 0) {
				json_object *geometry = json_hash_get(j, "geometry");
				int has_coords = 0;

				if (geometry != NULL && geometry->type == JSON_HASH) {
					json_object *geom_type = json_hash_get(geometry, "type");

					if (geom_type != NULL && geom_type->type == JSON_STRING && strcmp(geom_type->string, "LineString") == 0) {
						json_object *coordinates = json_hash_get(geometry, "coordinates");

						if (coordinates != NULL && coordinates->type == JSON_ARRAY) {
							for (size_t i = 0; i + 1 < coordinates->length; i++) {
								if (coordinates->array[i]->type == JSON_ARRAY &&
								    coordinates->array[i + 1]->type == JSON_ARRAY &&
								    coordinates->array[i]->length >= 2 &&
								    coordinates->array[i + 1]->length >= 2 &&
								    coordinates->array[i]->array[0]->type == JSON_NUMBER &&
								    coordinates->array[i]->array[1]->type == JSON_NUMBER &&
								    coordinates->array[i + 1]->array[0]->type == JSON_NUMBER &&
								    coordinates->array[i + 1]->array[1]->type == JSON_NUMBER) {
									index_add1(ix,
										   coordinates->array[i]->array[1]->number,
										   coordinates->array[i]->array[0]->number,
										   coordinates->array[i + 1]->array[1]->number,
										   coordinates->array[i + 1]->array[0]->number);

									has_coords = 1;
								}
							}
						}
					}
				}

				if (has_coords) {
					if (seq++ % (1000000 / 100) == 0) {
						fprintf(stderr, "%.2f million\r", seq / 1000000.0);
					}
				}

				json_free(j);
			} else if (type != NULL && type->type == JSON_STRING && strcmp(type->string, "FeatureCollection") == 0) {
				json_free(j);
			}
		}
	}

	index_sort(ix);
	seq = 0;

	while ((j = json_read(jp)) != NULL) {
		if (j->type == JSON_STRING && strcmp(j->string, "end part 1") == 0) {
			break;
		}

		if (j->type == JSON_HASH) {
			json_object *type = json_hash_get(j, "type");

			if (type != NULL && type->type == JSON_STRING && strcmp(type->string, "Feature") == 0) {
				json_object *properties = json_hash_get(j, "properties");
				char *props = NULL;
				if (properties != NULL) {
					props = json_stringify(properties);
				}

				json_object *geometry = json_hash_get(j, "geometry");

				if (geometry != NULL && geometry->type == JSON_HASH) {
					json_object *geom_type = json_hash_get(geometry, "type");

					if (geom_type != NULL && geom_type->type == JSON_STRING && strcmp(geom_type->string, "LineString") == 0) {
						json_object *coordinates = json_hash_get(geometry, "coordinates");

						if (coordinates != NULL && coordinates->type == JSON_ARRAY) {
							for (size_t i = 0; i + 1 < coordinates->length; i++) {
								if (coordinates->array[i]->type == JSON_ARRAY &&
								    coordinates->array[i + 1]->type == JSON_ARRAY &&
								    coordinates->array[i]->length >= 2 &&
								    coordinates->array[i + 1]->length >= 2 &&
								    coordinates->array[i]->array[0]->type == JSON_NUMBER &&
								    coordinates->array[i]->array[1]->type == JSON_NUMBER &&
								    coordinates->array[i + 1]->array[0]->type == JSON_NUMBER &&
								    coordinates->array[i + 1]->array[1]->type == JSON_NUMBER) {
									compare(ix,
										coordinates->array[i]->array[1]->number,
										coordinates->array[i]->array[0]->number,
										coordinates->array[i + 1]->array[1]->number,
										coordinates->array[i + 1]->array[0]->number, props ? props : "{}");

									if (seq++ % 10000 == 0) {
										fprintf(stderr, "checked %.3f million\r", seq / 1000000.0);
									}
								}
							}
						}
					}
				}

				if (props != NULL) {
					free(props);
				}

				json_free(j);
			} else if (type != NULL && type->type == JSON_STRING && strcmp(type->string, "FeatureCollection") == 0) {
				json_free(j);
			}
		}
	}

	json_end(jp);

	index_destroy(ix);
	return 0;
}
