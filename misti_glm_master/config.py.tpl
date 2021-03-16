ptf ~
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

OUTPUT_FILE_NAME = "output_climax_misti.txt"  # writes (x y thickness) data

# INPUT FILE FOR 'CUSTOM POINTS' RUN MODE (e.g. must contain (x y) data)
CUSTOM_POINTS_FILE = "custom_x_y_climax.txt"

# MAP AREA AND VENT COORDINATES
VENT_HEIGHT = 5800
TOPO_HEIGHT = 2200
VENT_EASTING = 242935
VENT_NORTHING = 8196472
MIN_EASTING = 211000 
MAX_EASTING = 265000 
MIN_NORTHING = 8163000 
MAX_NORTHING = 8219000 

# GROUND GRID RESOLUTION (DEFAULT MAP HAS 1000 GRID CELLS)
GROUND_GRID_SIZE = 100**2

# UMBRELLA CLOUD GEOMETRY PARAMETERS
# Ellipse full axis length, NOT semi-axis
ELLIPSE_MAJOR_AXIS =  ~  ellipse_major_axis  ~ 
ELLIPSE_MINOR_AXIS =  ~  ellipse_minor_axis  ~ 
ELLIPSE_GRID_STEP = 250

# ERUPTION PARAMETERS
TOTAL_ERUPTED_MASS =  ~  total_erupted_mass  ~ 
COLUMN_HEIGHT =  ~  column_height  ~ 
UMBRELLA_HEIGHT = COLUMN_HEIGHT + VENT_HEIGHT
BULK_DENSITY = 1500
PARTICLE_DENSITY_MAX = 1500
PARTICLE_DENSITY_MIN = 1500
DIFFUSION_COEF =  ~  diffusion_coef  ~ 
WIND_SPEED =  ~  wind_speed  ~ 
WIND_DIRECTION =  ~  wind_direction  ~ 
ETA = 1             # unitless

# DEFINE PHI INTERVALS (min | max phi according to TGSD)
MAX_PHI = -6
MIN_PHI = 6

TGSD_SIGMA =  ~  tgsd_sigma  ~ 
TGSD_MEAN =  ~  tgsd_mean  ~ 

# SIMULATED PHI INTERVALS (phi interval to be simulated)
SIMULATED_MAX_PHI = -6
SIMULATED_MIN_PHI = 6
STEP_PHI = 1

# CONSTANTS FOR SETTLING VELOCITY CALCULATION
ATMOSPHERE_LEVEL_STEP = 500  # 500 m thick atmosphere slices
RHO_A_S = 1.225  # density at sea level
G = 9.81  # gravitational constant
AIR_VISC = 1.8325e-5  # air viscosity in N s / m^2 at 24C

DISPERSAL_AXIS_POINTS_NUMBER = 500  # total number of points
DISPERSAL_AXIS_POINTS_DISTANCE = 100  # distance between two points
