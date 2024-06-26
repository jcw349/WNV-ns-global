from __future__ import print_function, division
import sys, re
from collections import defaultdict
import numpy as np

## augur export requires a TSV of lat/longs. This can create them.
## this allows the division lat/long to be dynamically generated from GPS points in the metadata

with open(sys.argv[1], 'r') as f:
  raw = f.read().splitlines() 
  header = raw[0].split("\t")
  metadata = []
  for line in raw[1:]:
    metadata.append({header[i]:data for i, data in enumerate(line.split("\t"))})
    

f = open(sys.argv[2], "w")

states = {
	"USA": (37.0902,	-95.7129),
	"NY": (40.7642499, -73.9545249),
	"NYC": (40.7095340064939, -73.93966854),
	"CA": (36.7014631, -118.7559973),
	"TX": (31.8160381, -99.5120985),
	"AZ": (34.395342, -111.7632754),
	"SD": (44.6471761, -100.3487609),
	"ND": (47.6201461, -100.5407369),
	"CO": (38.7251776, -105.6077166),
	"MD": (39.5162234, -76.9382068),
	"MN": (45.9896587, -94.6113287),
	"WI": (44.4308975, -89.6884636),
	"NJ": (40.0757384, -74.4041621),
	"CT": (41.6500201, -72.7342162),
	"AR": (35.2048883, -92.4479107),
	"MO": (38.7604815, -92.5617874),
	"IL": (40.0796319, -89.4339808),
	"PA": (40.9699889, -77.727883),
	"KY": (37.5726028, -85.155141),
	"MT": (47.3752671, -109.6387578),
	"IA": (41.9216734, -93.3122704),
	"AL": (33.2588817, -86.8295336),
	"VA": (37.1232245, -78.492772),
	"OK": (34.9550817, -97.2684062),
	"LA": (30.8703881, -92.0071259),
	"FL": (27.7567667, -81.4639834),
	"GA": (32.3293809, -83.1137365),
	"NV": (39.5158825, -116.8537226),
	"UT": (39.4225192, -111.7143583),
	"NM": (34.5708167, -105.9930069),
	"MS": (32.9715645, -89.7348496),
	"IN": (40.3270127, -86.1746932),
	"OH": (40.2253569, -82.6881394),
	"MI": (43.6211955, -84.6824345),
	"WY": (42.83, -107.48),
	"TN": (35.84, -86.85),
	"KS": (38.65, -98.38),
	"MA": (42.34, -72.11),
	"NC": (35.58, -80.1),
	"ID": (43.93, -114.24),
	"WA": (47.60, -122.33),
	"DC": (38.9, -77.0),
	"NE": (43.6211955, -84.6824345),
	"OR": (44, -120.5),
	"NH": (44, -71.5),
	"US-VI": (17.789187, -64.7080574),
        "MEX": (23.2127575,-113.2579052),
        "BRA": (-9.6711986,-99.7670981),
        "ARG": (-36.4741671,-84.8171707),
        "ISR": (31.4018434,33.7625888),
        "COL": (5.8422448,-85.0076846),
        "CAN": (54.2213959,-116.7597504),
	"CAN/QC": (45.508889, -73.561667),
	"MEX/BCN": (30, -115.166667),
	"MEX/CHH": (28.816667, -106.433333),
	"MEX/SON": (29.646111, -110.868889),
	"MEX/TAM": (24.283333, -98.566667),
	"VGB": (18.43596, -64.616942),
	"COL/ANT": (6.223333, -75.58),
	"BRA/ES": (-18.715833, -39.858889),
	"ARG/B": (-34.25, -59.466667),
	"ISR/D": (31.0, 35.0)
}

for key, value in states.items():
  f.write("{}\t{}\t{}\t{}\n".format("state", key, value[0], value[1]))

## DIVISIONS

# step 1: quick check to make sure there are no divisions of the same name in different states
div_state = {}
for x in metadata:
  div = x["division"]
  if div != "Unknown":
    if div in div_state:
      if div_state[div] != x["state"]:
        print("PROBLEM! {} - {} & {}".format(x["state"], div, div_state[div]))
    else:
      div_state[div] = x["state"]

# step 2: what are all the divisions and their corresponding GPS co-ords?
divisions_gps = defaultdict(lambda: [])
for x in metadata:
  if x["division"] != "Unknown":
    divisions_gps[x["division"]].append( [ x["latitude"], x["longitude"] ] )

# step 3: average the lat/longs:
divisions = {}
for key, values in divisions_gps.items():
  ll = list(filter(lambda x: "Unknown" not in x and "XXX" not in x, values))
  if len(ll) == 0:
    # no lat longs for these isolates, fall back to state instead
    state = key.split("/")[0]
    divisions[key] = states[state]
  else:
    divisions[key] = ( np.mean([float(x[0]) for x in ll]), np.mean([float(x[1]) for x in ll]) )

for key, value in divisions.items():
  f.write("{}\t{}\t{}\t{}\n".format("division", key, value[0], value[1]))
for key, value in states.items():
  f.write("{}\t{}\t{}\t{}\n".format("division", key, value[0], value[1]))
