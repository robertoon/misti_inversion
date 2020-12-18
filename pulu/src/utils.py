'''Utility functions to be used in specific calculations.'''
from math import sqrt, pow, sin, cos, radians


def distance_between_points(point1_x, point1_y, point2_x, point2_y):
    '''calculate distance between two points'''
    return sqrt((point2_x - point1_x)**2 + (point2_y - point1_y)**2)


def polar2cartesian(theta):
    '''convert polar to cartesian wind field data'''
    return (360 - theta + 90) % 360


def particle_diameter(phi):
    '''calculate particle diameter for given phi class'''
    return 2**(-1 * phi) * 10**-3


def get_ellipse_center(focus_x, focus_y, rx, ry, theta):
    theta = polar2cartesian(theta)
    dist_between_focus_and_center = sqrt(pow(rx/2, 2) - pow(ry/2, 2))
    cx = focus_x + dist_between_focus_and_center * cos(radians(theta))
    cy = focus_y + dist_between_focus_and_center * sin(radians(theta))
    return cx, cy


def translate_x(cx, cy, x, y, theta):
    theta = polar2cartesian(theta)
    return cos(radians(theta)) * (x - cx) - sin(radians(theta)) * (y - cy) + cx


def translate_y(cx, cy, x, y, theta):
    theta = polar2cartesian(theta)
    return sin(radians(theta)) * (x - cx) + cos(radians(theta)) * (y - cy) + cy


def is_point_in_ellipse(cx, cy, x, y, rx, ry):
    return ((pow((x - cx), 2) / pow(rx, 2)) +
            (pow((y - cy), 2) / pow(ry, 2)))
