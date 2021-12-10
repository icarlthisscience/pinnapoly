import csv

import matplotlib.pyplot as plt
import numpy as np

from matplotlib.patches import Polygon

# script could loop over several filenames to aregate data
# possibly, loaded from another txt or csv file  
filename = 'Sky92FNC_13 days_left pinna redo.csv'

with open(filename, newline = '') as f:
    reader = csv.DictReader(f)

    # load length points; first 2 lines in table
    pinna_extrema = []
    for _ in range(2):
        # read next line
        row = next(reader)
        # store X and Y column values
        point = (float(row['X']), float(row['Y']))
        pinna_extrema.append(point)

    # load ruler measurment; 3rd line in table
    for _ in range(1):
        # read next line
        row = next(reader)
        # store Length column value
        ruler_2cm = float(row['Length'])

    # load pinna outline; remaining lines in table
    pinna_points = []
    while True:
        try:
            # read next line
            row = next(reader)
            # store X and Y column values
            pinna_points.append((float(row['X']), float(row['Y'])))
        except StopIteration:
            # stop reading at end of file
            break

# convert data to matrix
length_data = np.array(pinna_extrema)
pinna_data = np.array(pinna_points)

# create shapes to be plotted
pinna_poly = Polygon(pinna_data, fill = False, color = 'm')
pinna_length = Polygon(length_data, fill = False, color = 'y')

# seperate domain and range (to be used to crop plot)
pinna_x = np.array(pinna_data.T)[0]
pinna_y = np.array(pinna_data.T)[1]

# create plot
ax = plt.subplot()
# crop plot with 25 pixel margins
ax.set_xlim((min(pinna_x) - 25, max(pinna_x) + 25))
ax.set_ylim((max(pinna_y) + 25), (min(pinna_y) - 25))
# note y axis inverted (max -> min)

# add shapes to plot
ax.add_patch(pinna_poly)
ax.add_patch(pinna_length)

# render
plt.show()
