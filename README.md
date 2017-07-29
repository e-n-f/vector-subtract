vector-subtract
===============

To find TIGER streets that are not in OSM:

```
(
	cat osm.geojson
	echo '"end part 1"'
	cat tiger.geojson
) | ./index > unmatched.geojson
```
