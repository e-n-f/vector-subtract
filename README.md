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

More detailed example, subtracting the Los Angeles OSM metro extract from the TIGER roads for Los Angeles County:

```
➤ unzip ../Downloads/tl_2017_06037_roads.zip
Archive:  ../Downloads/tl_2017_06037_roads.zip
 extracting: tl_2017_06037_roads.cpg
  inflating: tl_2017_06037_roads.dbf
  inflating: tl_2017_06037_roads.prj
  inflating: tl_2017_06037_roads.shp
  inflating: tl_2017_06037_roads.shp.ea.iso.xml
  inflating: tl_2017_06037_roads.shp.iso.xml
  inflating: tl_2017_06037_roads.shp.xml
  inflating: tl_2017_06037_roads.shx
➤ ogr2ogr -f GeoJSON tl_2017_06037_roads.json tl_2017_06037_roads.shp
➤ ( unzip -p ../Downloads/los-angeles_california.imposm-geojson.zip los-angeles_california_roads.geojson; echo '"end part 1"'; cat tl_2017_06037_roads.json ) | ./index > unmatched.geojson
➤ tippecanoe -zg --coalesce --drop-densest-as-needed --extend-zooms-if-still-dropping -f -o unmatched.mbtiles unmatched.geojson
```

(Leaves in place the roads in the northern part of the county that are not part of the metro extract.)
