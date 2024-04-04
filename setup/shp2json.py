
import fiona
import json
import glob, os

import sys



for file in glob.glob(sys.argv[1]+"/result/*.shp"):

    shp = fiona.open(file,"r")

    geojson = { "type": "FeatureCollection",
                "features": [] }

    for f in shp:
         geojson["features"].append(f)

    outf = file.split(".")[0]+".geojson"

    with open(outf, "w") as gjs:
        json.dump(geojson, gjs)
