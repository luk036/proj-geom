"""
This module contains functions that are common to projective plane.
"""


# Write a dot product of two vectors.
def dot_product(vector1, vector2):
    """
    The dot_product function calculates the dot product of two vectors.

    :Example:
    >>> dot_product([1, 2, 3], [4, 5, 6])
    32
    >>> dot_product([1, 2, 3], [4, 5, 6, 7])
    Traceback (most recent call last):
    ...
    ValueError: Vectors must be the same length
    """
    if len(vector1) != len(vector2):
        raise ValueError("Vectors must be the same length")
    return sum(vector1[i] * vector2[i] for i in range(len(vector1)))


# Write a cross product of two vectors.
def cross_product(vector1, vector2):
    """
    The function `cross_product` calculates the cross product of two 3D vectors.

    :Example:
    >>> cross_product([1, 2, 3], [4, 5, 6])
    [-3, 6, -3]
    >>> cross_product([1, 2, 3], [4, 5, 6, 7])
    Traceback (most recent call last):
    ...
    ValueError: Vectors must be the same length
    >>> cross_product([1, 2, 3], [4, 5])
    Traceback (most recent call last):
    ...
    ValueError: Vectors must be the same length
    """
    if len(vector1) != len(vector2):
        raise ValueError("Vectors must be the same length")
    if len(vector1) != 3:
        raise ValueError("Vectors must be 3D")
    return [
        vector1[1] * vector2[2] - vector1[2] * vector2[1],
        vector1[2] * vector2[0] - vector1[0] * vector2[2],
        vector1[0] * vector2[1] - vector1[1] * vector2[0],
    ]

# Write a parametric line equation.
def parametric_line_equation(point1, point2):
    """
    The function `parametric_line_equation` calculates the parametric line equation of two points.

    :Example:
    >>> parametric_line_equation([1, 2, 3], [4, 5, 6])
    'x = 4, y = 5, z = 6'
    >>> parametric_line_equation([1, 2, 3], [4, 5, 6, 7])
    Traceback (most recent call last):
    ...
    ValueError: Points must be the same length
    """
    if len(point1) != len(point2):
        raise ValueError("Points must be the same length")
    return "x = " + str(point2[0]) + ", y = " + str(point2[1]) + ", z = " + str(point2[2])