from __future__ import print_function, division
import numpy as np
import matplotlib as mpl
from collections import defaultdict
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import sys

# This script attempts to construct a colour palette for the states, locations etc used by WNV nextstrain.
# It is called by Snakemake, just before the final step (export)
# ARG1: METADATA ARG2: OUTPUT

with open(sys.argv[1], 'r') as f:
  raw = f.read().splitlines() 
  header = raw[0].split("\t")
  metadata = []
  for line in raw[1:]:
    metadata.append({header[i]:data for i, data in enumerate(line.split("\t"))})

fh = open(sys.argv[2], "w")

#########################################################################################
# COUNTRIES
countries    = ["United-States", "Mexico", "Israel", "British-Virgin-Islands", "Canada", "Colombia", "Brazil", "Argentina", "US-Virgin-Islands","Unknown"]
country_cols = ["#969696", "#7E00A8", "#0054A8", "#387A47", "#7BDE00", "#E5C800", "#FF7F00", "#A80000", "#4C2F00","#00E5E5"]
fh.write("## COUNTRIES ##\n")
for pair in zip(countries, country_cols):
  fh.write("{}\t{}\t{}\n".format("country", pair[0], pair[1]))

#########################################################################################
# LINEAGES / STRAINS / CLADES (they're not monophyletic)
#wnv_strain       = ["NY99",    "SW03",    "WN02", "pre-NY","Unknown"]
#wnv_strain_cols  = ["#CBB742", "#7EB876", "#4988C5", "#A80000","#00E5E5"]
strains = ['kunjin', 'Kunjin MRM61C', 'NY99', 'NY99a', 'NY99b', 'RabV', 'SW03', 'WN02', 'Unknown']
strain_cols = ['#44AA99', '#44AA99', '#332288', '#117733', '#88CCEE', '#AA4499', '#DDCC77', '#882255', '#A0A0A0']
fh.write("## WNV STRAINS / LINEAGES ##\n")
for pair in zip(strains, strain_cols):
	fh.write("{}\t{}\t{}\n".format("strains", pair[0], pair[1]))

#########################################################################################
# LINEAGES / STRAINS / CLADES (they're not monophyletic)
clades = ['1', '1a', '1b', '2', '3', '4', '7', 'Unknown']
clade_cols = ['#882255', '#AA4499', '#DDCC77', '#332288', '#117733', '#44AA99', '#88CCEE', '#A0A0A0']
#wnv_clades       = ["1a",    "1b",    "2", "3", "4","Unknown"]
#wnv_clade_cols  = ["#969696", "#7E00A8", "#0054A8", "#387A47", "#7BDE00","#00E5E5"]
fh.write("## WNV CLADES / CLADES ##\n")
for pair in zip(clades, clade_cols):
	fh.write("{}\t{}\t{}\n".format("clade_membership", pair[0], pair[1]))
    
#########################################################################################
# STATES
# different colour scales where used to colour states based on five main US geographic regions
states = {
	"west": ["AK", "WA", "ID", "MT", "OR", "NV", "WY", "CA", "UT", "CO", "HI"],
	"southwest": ["AZ", "NM", "OK", "TX"],
	"midwest": ["ND", "MN", "IL", "WI", "MI", "SD", "IA", "IN", "OH", "NE", "MO", "KS"],
	"southeast": ["KY", "WV", "VA", "AR", "TN", "NC", "SC", "LA", "MS", "AL", "GA", "FL"],
	"northeast": ["ME", "VT", "NH", "NY", "MA", "RI", "PA", "NJ", "CT", "MD", "DC", "DE"],
	"caribbean": ["US-VI", "VGB"],
	"canada": ["CAN/QC"],
	"mexico": ["MEX/BCN", "MEX/SON", "MEX/CHH", "MEX/TAM"],
	"southamerica": ["COL/ANT", "BRA/ES", "ARG/B"],
	"israel": ["ISR/D"]
}
# prune out the states that _aren't_ in the metadata (before we create the colour scale)
# states_present = set([x["state"] for x in metadata])
# for key, values in states.items():
#   states[key] = list(filter(lambda x: x in states_present, values))
  # generate color maps
states_cols = {
	"west": ["#590000", "#690909", "#7A1515", "#8A2424", "#9B3636", "#AB4B4B", "#BC6363", "#CD7E7E", "#DD9C9C", "#EEBDBD", "#FFE1E1"],
	"southwest": ["#AD4700", "#C36A11", "#D98D22", "#EFB034"],
	"midwest": ["#001C00", "#022B02", "#063906", "#0C480C", "#145814", "#1D671D", "#297529", "#358535", "#449344", "#54A354", "#67B267", "#7BC17B"],
	"southeast": ["#001833", "#032142", "#092C51", "#113860", "#1A466F", "#25547E", "#32638D", "#41729C", "#5183AB", "#6495BA", "#78A7C9", "#8EBAD9"],
	"northeast": ["#1C001C", "#2E0330", "#3F0B45", "#501559", "#60236E", "#713483", "#824997", "#9460AC", "#A77BC1", "#BD9AD5", "#D5BBEA", "#F0E1FF"],
	"caribbean": ["#4C2F00", "#7F6233"],
	"canada": ["#7FD000"],
	"mexico": ["#4E5900", "#7D882E", "#ACB75C", "#DBE78A"],
	"southamerica": ["#00594E", "#00A08C", "#00E7CA"],
	"israel": ["#E533E5"],
	"Unknown": ["#00E5e5"]
}

fh.write("## STATES ##\n")
for category, names in states.items():
  for pair in zip(names, states_cols[category]):
    fh.write("{}\t{}\t{}\n".format("state", pair[0], pair[1]))

#########################################################################################
# HOSTS
# use the tab20c scale - https://matplotlib.org/examples/color/colormaps_reference.html
# tab20     = [mpl.colors.rgb2hex(mpl.cm.tab20c(i)) for i in range(0,20)]
#host_categories = ["Bird-crow", "Bird-other", "Bird-unknown", "Mosquito-Aedes", "Mosquito-Culex", "Mosquito-Culiseta", "Mosquito-other", "Mosquito-unknown", "Human", "Horse", "Squirrel", "Unknown"]
#host_cols = ["#000000", "#41ab5d", "#addd8e","#969696", "#8c96c6", "#8c6bb1","#88419d", "#810f7c", "#B20023","#D86239", "#FEC34F", "#00E5E5"]
host_categories = ['Birds', 'Birds-Hawks and Eagles', 'Birds-Landfowls', 'Birds-Other', 'Birds-Owls', 'Birds-Parrot', 'Birds-Perching Birds', 'Birds-Pigeons and Doves', 'Birds-Seabirds', 'Birds-Shorebirds', 'Birds-Vultures', 'Birds-Wading Birds', 'Birds-Waterfowls', 'Mammals-Alpaca', 'Mammals-Bat', 'Mammals-Camels', 'Mammals-Dogs', 'Mammals-Donkey', 'Mammals-Hamster', 'Mammals-Horse', 'Mammals-Human', 'Mammals-Rodents', 'Mosquitoes', 'Mosquitoes-Aedes', 'Mosquitoes-Anopheles', 'Mosquitoes-Coquillettidia', 'Mosquitoes-Culex', 'Mosquitoes-Culiseta', 'Mosquitoes-Mansonia', 'Mosquitoes-Ochlerotatus', 'Mosquitoes-Other', 'Mosquitoes-Psorophora', 'Mosquitoes-Unknown', 'Mosquitoes-Uranotaenia', 'Reptiles-Crocodile', 'Ticks', 'Unknown']
host_cols = ['#08306b', '#08519c', '#2171b5', '#4292c6', '#6baed6', '#9ecae1', '#c6dbef', '#08306b', '#08519c', '#08519c', '#2171b5', '#08519c', '#08519c', '#f7fcf5', '#006d2c', '#238b45', '#41ab5d', '#74c476', '#a1d99b', '#c7e9c0', '#00441b', '#a1d99b', '#67000d', '#67000d', '#cb181d', '#ef3b2c', '#fb6a4a', '#ffffcc', '#fcbba1', '#fee0d2', '#fee0d2', '#fee0d2', '#fee0d2', '#fee0d2', '#FFFFFF', '#fff5f0', '#A0A0A0']
fh.write("## HOSTS ##\n")
for pair in zip(host_categories, host_cols):
  fh.write("{}\t{}\t{}\n".format("host_categories", pair[0], pair[1]))

#########################################################################################
# DIVISIONS
divisions = set(x["division"] for x in metadata)
divisions.discard("Unknown")
divisions.discard("division")
divisions = list(divisions)
divisions.sort(key=str.lower)
divisions_cols = [mpl.colors.rgb2hex(mpl.cm.viridis(i)) for i in np.linspace(0,1,len(divisions))]
fh.write("## DIVISIONS ##\n")
for pair in zip(divisions, divisions_cols):
  fh.write("{}\t{}\t{}\n".format("division", pair[0], pair[1]))

#########################################################################################

# LOCATIONS
location = ["NYC", "NYC-PHL", "NY","USA","Global"]
location_cols = ["#CC6677", "#CC6677", "#501559", "#C5E021", "#00E5E5"]
fh.write("## LOCATIONS ##\n")
for pair in zip(location, location_cols):
  fh.write("{}\t{}\t{}\n".format("location", pair[0], pair[1]))

#########################################################################################

#########################################################################################

# SOURCE
source = ["NYC-PHL", "NS_NA_NCBI", "NS_NA_WN4K","WN4K","NCBI_Global"]
source_cols = ["#CC6677", "#C5E021", "#C5E021", "#8382DC", "#00E5E5"]
fh.write("## LOCATIONS ##\n")
for pair in zip(location, location_cols):
  fh.write("{}\t{}\t{}\n".format("DataSource", pair[0], pair[1]))

#########################################################################################

# END
fh.close()
