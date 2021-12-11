import csv
import math

import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry import Polygon as Area

from matplotlib.patches import Polygon

# script could loop over several filenames to aregate data
# possibly, loaded from another txt or csv file  
filename = 'Sky92FNC_13 days_left pinna redo.csv'

def average(x):
    return sum(x)/len(x)

def direction(point):

    return math.degrees(math.atan2(point[1], point[0]))

def rotate(origin, point, angle):
    
    ox, oy = origin
    px, py = point

    dx = px - ox
    dy = py - oy
    dcos = math.cos(math.radians(angle)) 
    dsin = math.sin(math.radians(angle)) 

    qx = ox + (dcos * dx) - (dsin * dy)
    qy = oy + (dsin * dx) + (dcos * dy)

    return qx, qy

def max_width(data):
    
    y_data = [p[1] for p in data]
    min_y = min(y_data)
    max_y = max(y_data)

    # find maximum
    i = 0
    start = 0
    while True:
        p = data[i][1]
        if (p == max_y): 
            start = i
            break

        i += 1
    
    # collect points until min
    j = i
    left_edge = []
    while data[j][1] > min_y:
        # remove convex segments 
        k = 1
        while i - k > start and data[j - k][1] < data[j][1]:
            k += 1
        while k > 1:
            left_edge.pop()
            k -= 1

        left_edge.append(data[j])

        i += 1
        j = i if i < len(data) else i - len(data)
    
    left_edge.append(data[j])
    checkpoint = i - 1

    # collect points until max
    right = [left_edge[-1]]
    while data[j][1] < max_y:
        # remove convex segments 
        k = 1
        while i - k > checkpoint and data[j - k][1] > data[j][1]:
            k += 1
        while k > 1:
            right.pop()
            k -= 1
        
        right.append(data[j])

        i += 1
        j = i if i < len(data) else i - len(data)

    right.append(data[start])

    left_map = {p[1] : p[0] for p in left_edge}
    right_map = {p[1] : p[0] for p in right}

    # interpolate critical points
    left_int = {}
    y1 = None
    for y2 in left_map:
        if y1 is None:
            y1 = y2
            left_int[y1] = left_map[y1]
            continue
        
        # create left edge node for each right edge node
        for r in right_map:
            if (y1 < r < y2) or (y1 > r > y2):
                x1 = left_map[y1]
                x2 = left_map[y2]

                y = r
                x = x1 + ((x2 - x1) * ((y - y1)/(y2 - y1)))

                left_int[y] = x

        y1 = y2
        left_int[y1] = left_map[y1]
    
    # interpolate critical points
    right_int = {}
    y1 = None
    for y2 in right_map:
        if y1 is None:
            y1 = y2
            right_int[y1] = right_map[y1]
            continue

        # create right edge node for each left edge node
        for l in left_map:
            if (y1 < l < y2) or (y1 > l > y2):
                
                x1 = right_map[y1]
                x2 = right_map[y2]

                y = l
                x = x1 + ((x2 - x1) * ((y - y1)/(y2 - y1)))

                right_int[y] = x

        y1 = y2
        right_int[y1] = right_map[y1]

    # ensure points are ordered
    left_ordered = dict(sorted(left_int.items(), key = lambda item: -item[0]))
    right_ordered = dict(sorted(right_int.items(), key = lambda item: item[0]))

    # find max width at critical points
    y = max_y
    width = 0
    for l in left_ordered:
        w = abs(left_ordered[l] - right_ordered[l])
        if w > width:
            y = l
            width = w

    return width, y

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

# seperate domain and range (to be used to crop plot)
pinna_x = np.array(pinna_data.T)[0]
pinna_y = np.array(pinna_data.T)[1]

# rotate pinna around arbitrary centre
centre = [average(pinna_x), average(pinna_y)]
pinna_rotation = direction([
    length_data[1][0] - length_data[0][0],
    length_data[1][1] - length_data[0][1],
])
rotated_points = [rotate(centre, p, 90 - pinna_rotation) for p in pinna_points]
rotated_length = [rotate(centre, p, 90 - pinna_rotation) for p in pinna_extrema]

# calculate area
area = Area(pinna_points)
print(area.area)

# calculate width
width, cross_y = max_width(rotated_points)
print(width, cross_y)

# project width cross section (from rotated data)
rotated_width = [
    [min(p[0] for p in rotated_points), cross_y],
    [max(p[0] for p in rotated_points), cross_y]
]
# unrotate width cross section
pinna_width = [rotate(centre, p, pinna_rotation - 90) for p in rotated_width]

# create shapes to be plotted
pinna_poly = Polygon(pinna_points, fill = False, color = 'black')
pinna_poly_rotated = Polygon(rotated_points, fill = False, color = 'grey')
pinna_length = Polygon(pinna_extrema, fill = False, color = 'darkgreen')
pinna_length_rotated = Polygon(rotated_length, fill = False, color = 'lightgreen')
pinna_width = Polygon(pinna_width, fill = False, color = 'darkgreen')
pinna_width_rotated = Polygon(rotated_width, fill = False, color = 'lightgreen')

# create plot
ax = plt.subplot()
# crop plot with 25 pixel margins
ax.set_xlim((min(pinna_x) - 25, max(pinna_x) + 25))
ax.set_ylim((max(pinna_y) + 25), (min(pinna_y) - 25))
# note y axis inverted (max -> min)

# add shapes to plot
ax.add_patch(pinna_poly)
ax.add_patch(pinna_poly_rotated)
ax.add_patch(pinna_length)
ax.add_patch(pinna_length_rotated)
ax.add_patch(pinna_width)
ax.add_patch(pinna_width_rotated)

# render
plt.show()
