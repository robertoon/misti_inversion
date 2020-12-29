import numpy as np
from math import sqrt, cos, sin, radians
from src.tephra import mass_fraction, particle_density
from src.utils import (
    distance_between_points,
    polar2cartesian,
    get_ellipse_center,
    translate_x,
    translate_y,
    is_point_in_ellipse
)
from src.config import (
    VENT_EASTING,
    VENT_NORTHING,
    WIND_DIRECTION,
    ELLIPSE_MAJOR_AXIS,
    ELLIPSE_MINOR_AXIS,
    ELLIPSE_GRID_STEP,
    CUSTOM_POINTS_FILE,
    GROUND_GRID_SIZE,
    DISPERSAL_AXIS_POINTS_DISTANCE,
    DISPERSAL_AXIS_POINTS_NUMBER
)


def generate_ellipse_grid(focus_x, focus_y, rx, ry, theta, grid_step):
    c_x, c_y = get_ellipse_center(
        VENT_EASTING, VENT_NORTHING,
        ELLIPSE_MAJOR_AXIS, ELLIPSE_MINOR_AXIS,
        WIND_DIRECTION
    )

    BOTTOM_LEFT_X = c_x - ELLIPSE_MAJOR_AXIS/2
    BOTTOM_LEFT_Y = c_y - ELLIPSE_MINOR_AXIS/2

    grid_cells = []

    for i in range(0, int(ELLIPSE_MAJOR_AXIS/ELLIPSE_GRID_STEP)):
        for j in range(0, int(ELLIPSE_MINOR_AXIS/ELLIPSE_GRID_STEP)):
            x = BOTTOM_LEFT_X + ELLIPSE_GRID_STEP/2 * (2*i + 1)
            y = BOTTOM_LEFT_Y + ELLIPSE_GRID_STEP/2 * (2*j + 1)

            r_x = translate_x(c_x, c_y, x, y, WIND_DIRECTION)
            r_y = translate_y(c_x, c_y, x, y, WIND_DIRECTION)

            if is_point_in_ellipse(c_x, c_y, x, y, ELLIPSE_MAJOR_AXIS/2, ELLIPSE_MINOR_AXIS/2) <= 1:
                grid_cells.append({'x': r_x, 'y': r_y})

    return grid_cells


def generate_ground_grid(min_x, max_x, min_y, max_y):
    '''generate the ground grid on which tephra load is calculated'''
    step = int((max_x - min_x) / sqrt(GROUND_GRID_SIZE))

    return (
        {'x': x, 'y': y}
        for x in range(min_x, max_x, step)
        for y in range(min_y, max_y, step)
    )


def get_custom_points():
    '''
    read (x y) points from input file of locations of interest
    for tephra load calculation
    '''
    point_x, point_y = np.genfromtxt(CUSTOM_POINTS_FILE, unpack=True)

    return (
        {'x': x, 'y': y} for x, y in zip(point_x, point_y)
    )


def get_dispersal_axis_points(x_vent, y_vent, alpha):
    alpha = polar2cartesian(alpha)
    points = []
    for i in range(DISPERSAL_AXIS_POINTS_NUMBER):
        distance = DISPERSAL_AXIS_POINTS_DISTANCE * (i + 1)
        x = x_vent + (distance * cos(radians(alpha)))
        y = y_vent + (distance * sin(radians(alpha)))
        points.append({'x': x, 'y': y})
    return points
