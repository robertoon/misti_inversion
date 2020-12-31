'''
Configuration file for pyTephra2 code with umbrella cloud geometry.
All units are in international system units (i.e. kilograms, meters,
meters per second, kilograms per square meter).
The user should input the simulation values for each parameter.
Two run modes are available with this code: 'grid' (calculates f(x,y) at
each grid location on the map) and 'custom points' (calculates f(x,y) at
locations of interest (e.g. field points) given an input file).

AUTHORS: R. CONSTANTINESCU AND H.G. AURELIAN
UPDATE: March 20, 2020
'''

# SELECT RUN MODE: - DEFAULT VALUES ARE 'grid', 'custom points' OR 'dispersal axis points'
RUN_MODE = "custom points"

OUTPUT_FILE_NAME = "output_misti.txt"  # writes (x y thickness) data

# INPUT FILE FOR 'CUSTOM POINTS' RUN MODE (e.g. must contain (x y) data)
CUSTOM_POINTS_FILE = "custom_x_y_misti.txt"

# MAP AREA AND VENT COORDINATES
VENT_EASTING = 242782
VENT_NORTHING = 8196548
MIN_EASTING = VENT_EASTING - 50000
MAX_EASTING = VENT_EASTING + 50000
MIN_NORTHING = VENT_NORTHING - 50000
MAX_NORTHING = VENT_NORTHING + 50000


# GROUND GRID RESOLUTION (DEFAULT MAP HAS 1000 GRID CELLS)
GROUND_GRID_SIZE = 100**2

# UMBRELLA CLOUD GEOMETRY PARAMETERS
# Ellipse full axis length, NOT semi-axis
ELLIPSE_MAJOR_AXIS =  1.150000000000000000e+04 
ELLIPSE_MINOR_AXIS =  1.150000000000000000e+04 
ELLIPSE_GRID_STEP = 250

# ERUPTION PARAMETERS
TOTAL_ERUPTED_MASS =  8.200000000000000000e+10 
COLUMN_HEIGHT =  2.1000000000000e+04 
BULK_DENSITY = 1500
PARTICLE_DENSITY_MAX = 1500
PARTICLE_DENSITY_MIN = 1500
DIFFUSION_COEF =  5.00000000000000e+03 
WIND_SPEED =  5.00000000000000 
WIND_DIRECTION =  2.17000000000000e+02 
ETA = 1             # unitless

# DEFINE PHI INTERVALS (min | max phi according to TGSD)
MAX_PHI = -5
MIN_PHI = 5

TGSD_SIGMA =  1.00000000000000 
TGSD_MEAN =  -1.990000000000 

# SIMULATED PHI INTERVALS (phi interval to be simulated)
SIMULATED_MAX_PHI = -5
SIMULATED_MIN_PHI = 5
STEP_PHI = 1

# CONSTANTS FOR SETTLING VELOCITY CALCULATION
ATMOSPHERE_LEVEL_STEP = 500  # 500 m thick atmosphere slices
RHO_A_S = 1.225  # density at sea level
G = 9.81  # gravitational constant
AIR_VISC = 1.8325e-5  # air viscosity in N s / m^2 at 24C

DISPERSAL_AXIS_POINTS_NUMBER = 500  # total number of points
DISPERSAL_AXIS_POINTS_DISTANCE = 100  # distance between two points
