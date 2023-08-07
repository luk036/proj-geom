"""This module contains functions that are common to projective plane."""

from typing import List, Sequence
from typing import Generic, TypeVar
from abc import abstractmethod
from typing_extensions import Self

Dual = TypeVar("Dual", bound="ProjectivePlane")
Measure = TypeVar("Measure", bound=int)


class ProjectivePlane(Generic[Dual, Measure]):
    """
    The `ProjectivePlane` class is a generic class that represents a projective plane.
    """

    @abstractmethod
    def dual(self) -> type:
        """
        The function "dual" is a placeholder function that does not have any implementation.
        """

    @abstractmethod
    def __eq__(self, rhs) -> bool:
        """
        The above function is an abstract method that defines the equality comparison for a class.
        
        :param rhs: The parameter "rhs" stands for "right-hand side" and represents the object that is
        being compared to the current object
        """

    @abstractmethod
    def circ(self, rhs: Self) -> Dual:
        """
        The function circ takes two arguments, self and rhs, and returns a Dual object.
        
        :param rhs: The parameter "rhs" stands for "right-hand side" and it is of type "Self"
        :type rhs: Self
        """

    @abstractmethod
    def aux(self) -> Dual:
        """
        The function `aux` is a method that returns a `Dual` object.
        """

    @abstractmethod
    def dot(self, line_l: Dual) -> Measure:
        """
        The dot function takes in a Dual object and returns a Measure object.
        
        :param line_l: The parameter "line_l" is of type Dual
        :type line_l: Dual
        """

    @abstractmethod
    def plucker(self, lambda_p: Measure, pt_q: Self, mu_q: Measure) -> Self:
        """
        The plucker function takes in three parameters (lambda_p, pt_q, mu_) and returns a value of type Self.
        
        :param lambda_p: The parameter "lambda_p" is of type Measure
        :type lambda_p: Measure
        :param pt_q: The parameter "pt_q" is of type "Self", which means it is an instance of the same class
        as the method it is being used in
        :type pt_q: Self
        :param mu_q: The parameter "mu_q" is of type Measure
        :type mu_q: Measure
        """

    @abstractmethod
    def incident(self, line_l: Dual) -> bool:
        """
        The function checks if a point is incident to a line.
        
        :param line_l: The parameter "line_l" is of type "Dual"
        :type line_l: Dual
        :return: a boolean value. It returns True if the dot product of the self vector and the line
        vector is equal to 0, and False otherwise.
        """
        return self.dot(line_l) == 0

    def coincident(self, pt_q: Self, pt_r: Self) -> bool:
        """
        The `coincident` function checks if three points `pt_p`, `pt_q`, and `pt_r` are collinear.

        :param pt_q: pt_q is an instance of the class ProjectivePlanePrim<Point>
        :type pt_q: Self
        :param pt_r: The parameter `pt_r` is of type `ProjectivePlanePrim<Point>`
        :type pt_r: Self
        :return: A boolean value is being returned.
        """
        return self.circ(pt_q).incident(pt_r)

    def harm_conj(self, pt_a: Self, pt_b: Self) -> Self:
        """
        The `harm_conj` function calculates the harmonic conjugate of two points on a projective plane.

        :param pt_a: The parameter `pt_a` is of type `ProjectivePlane`
        :type pt_a: Self
        :param pt_b: The parameter `pt_b` is of type `ProjectivePlane`
        :type pt_b: Self
        :return: a ProjectivePlane object.
        """
        assert self.coincident(pt_a, pt_b)
        line_ab = pt_a.circ(pt_b)
        line_l = line_ab.aux().circ(self)
        return pt_a.plucker(line_l.dot(pt_b), pt_b, line_l.dot(pt_a))


Point = ProjectivePlane["Line", Measure]
Line = ProjectivePlane["Point", Measure]

# ProjectivePlanePrim = PgObject

# trait ProjectivePlanePrim<Line>: Eq {
#     def circ(self, rhs: Self) -> Line
#     def incident(self, line) -> bool
# }


def check_axiom(pt_p: Point, pt_q: Point, line_l: Line):
    """
    The function `check_axiom` checks various axioms related to a projective plane.

    :param pt_p: pt_p is a ProjectivePlanePrim object, which represents a point in a projective plane
    :type pt_p: Point
    :param pt_q: The parameter `pt_q` is a ProjectivePlanePrim object, which represents a point or a line in a projective plane
    :type pt_q: Point
    :param line: The `line` parameter represents a projective plane line
    :type line: Line
    """
    assert pt_p == pt_p
    assert (pt_p == pt_q) == (pt_q == pt_p)
    assert pt_p.incident(line_l) == line_l.incident(pt_p)
    assert pt_p.circ(pt_q) == pt_q.circ(pt_p)
    line_m = pt_p.circ(pt_q)
    assert line_m.incident(pt_p) and line_m.incident(pt_q)


def coincident(pt_p: Point, pt_q: Point, pt_r: Point) -> bool:
    """
    The `coincident` function checks if three points `pt_p`, `pt_q`, and `pt_r` are collinear in a projective
    plane.

    :param pt_p: pt_p is an object of type ProjectivePlanePrim<Point>. It represents a point in a projective plane
    :type pt_p: Point
    :param pt_q: pt_q is a ProjectivePlanePrim object, which represents a point in a projective plane
    :type pt_q: Point
    :param pt_r: The parameter `pt_r` represents a point in a projective plane
    :type pt_r: Point
    :return: The function `coincident` returns a boolean value.

    Examples:
        >>> from projgeom_py.pg_point import PgLine, PgPoint
        >>> coincident(PgPoint([0, 1, 0]), PgPoint([0, 0, 1]), PgPoint([1, 0, 0]))
        False
    """
    return pt_p.circ(pt_q).incident(pt_r)


def check_pappus(co1: List[Point], co2: List[Point]) -> bool:
    """
    The function `check_pappus` checks if three lines in a projective plane satisfy Pappus' theorem.

    :param co1: The parameter `co1` is a list of `ProjectivePlanePrim` objects
    :type co1: List[Point]
    :param co2: The parameter `co2` is a list of `ProjectivePlanePrim` objects
    :type co2: List[Point]
    :return: a boolean value.

    Examples:
        >>> from projgeom_py.pg_point import PgLine, PgPoint
        >>> co1 = [PgPoint([0, 1, 0]), PgPoint([0, 0, 1]), PgPoint([1, 0, 0])]
        >>> co2 = [PgPoint([0, 0, 1]), PgPoint([0, 1, 0]), PgPoint([1, 0, 0])]
        >>> check_pappus(co1, co2)
        True
    """
    [pt_a, pt_b, pt_c] = co1
    [pt_d, pt_e, pt_f] = co2
    pt_g = (pt_a.circ(pt_e)).circ(pt_b.circ(pt_d))
    pt_h = (pt_a.circ(pt_f)).circ(pt_c.circ(pt_d))
    pt_i = (pt_b.circ(pt_f)).circ(pt_c.circ(pt_e))
    return coincident(pt_g, pt_h, pt_i)


def tri_dual(tri: Sequence) -> List:
    """
    The function `tri_dual` takes a list of three `ProjectivePlanePrim` objects representing a triangle and
    returns a list of three `ProjectivePlanePrim` objects representing the circumcircles of the triangle's
    three edges.

    :param tri: The `tri` parameter is expected to be a sequence (e.g., list, tuple) of three elements. Each element should be an object of type `ProjectivePlanePrim`
    :type tri: Sequence
    :return: The function `tri_dual` returns a list of three `ProjectivePlanePrim` objects.
    """
    [a_1, a_2, a_3] = tri
    assert not coincident(a_1, a_2, a_3)
    return [a_2.circ(a_3), a_1.circ(a_3), a_1.circ(a_2)]


def persp(tri1: List[Point], tri2: List[Point]) -> bool:
    """
    The `persp` function checks whether two triangles are perspective.

    :param tri1: tri1 is a list of three ProjectivePlanePrim objects representing the vertices of the first triangle
    :type tri1: List[Point]
    :param tri2: tri2 is a list of three ProjectivePlanePrim objects representing the vertices of the second triangle
    :type tri2: List[Point]
    :return: a boolean value.

    Examples:
        >>> from projgeom_py.pg_point import PgLine, PgPoint
        >>> tri1 = [PgPoint([0, 1, 0]), PgPoint([0, 0, 1]), PgPoint([1, 0, 0])]
        >>> tri2 = [PgPoint([0, 0, 1]), PgPoint([0, 1, 0]), PgPoint([1, 0, 0])]
        >>> persp(tri1, tri2)
        True
    """
    [pt_a, pt_b, pt_c] = tri1
    [pt_d, pt_e, pt_f] = tri2
    pt_o = pt_a.circ(pt_d).circ(pt_b.circ(pt_e))
    return pt_c.circ(pt_f).incident(pt_o)


def check_desargue(tri1: List[Point], tri2: List[Point]) -> bool:
    """
    The function `check_desargue` checks if two triangles in a projective plane satisfy the Desargue's
    theorem.

    :param tri1: tri1 is a list of ProjectivePlanePrim objects representing the first triangle in the Desargue's theorem
    :type tri1: List[Point]
    :param tri2: The `tri2` parameter is a list of `ProjectivePlanePrim` objects representing the second triangle
    :type tri2: List[Point]
    :return: a boolean value.

    Examples:
        >>> from projgeom_py.pg_point import PgLine, PgPoint
        >>> tri1 = [PgPoint([0, 1, 0]), PgPoint([0, 0, 1]), PgPoint([1, 0, 0])]
        >>> tri2 = [PgPoint([0, 0, 1]), PgPoint([0, 1, 0]), PgPoint([1, 0, 0])]
        >>> check_desargue(tri1, tri2)
        True
    """
    trid1 = tri_dual(tri1)
    trid2 = tri_dual(tri2)
    bool1 = persp(tri1, tri2)
    bool2 = persp(trid1, trid2)
    return (bool1 and bool2) or (not bool1 and not bool2)


# trait ProjectivePlane<Line, Measure: Default + Eq>: ProjectivePlanePrim<Line>:
#     def aux(self) -> Line
#     def dot(self, line) -> Measure; # basic measurement
#     def plucker(lambda_: Measure, pt_p: Self, mu_: Measure, pt_q: Self)
#     def incident(self, line) -> bool:
#         self.dot(line) == Measure::default()


def harm_conj(pt_a: Point, pt_b: Point, pt_c: Point):
    """
    The `harm_conj` function calculates the harmonic conjugate of three points on a projective plane.

    :param pt_a: a is an object of type ProjectivePlane
    :type pt_a: Point
    :param pt_b: The parameter `pt_b` represents a point on the projective plane
    :type pt_b: Point
    :param pt_c: The parameters `pt_a`, `pt_b`, and `pt_c` are of type `ProjectivePlane`
    :type pt_c: Point
    :return: The function `harm_conj` returns a `ProjectivePlane` object.
    """
    assert coincident(pt_a, pt_b, pt_c)
    ab = pt_a.circ(pt_b)
    lc = ab.aux().circ(pt_c)
    # Point = type(pt_a)
    return pt_a.plucker(lc.dot(pt_b), pt_b, lc.dot(pt_a))


def involution(origin: Point, mirror: Line, pt_p: Point):
    """
    The function `involution` performs an involution transformation on a point `pt_p` with respect to an
    origin point `origin` and a mirror line `mirror`.

    :param origin: The `origin` parameter represents a point in a projective plane
    :type origin: Point
    :param mirror: The `mirror` parameter represents a mirror line or mirror plane in a projective plane. It is used to perform a reflection or mirror transformation on a point `pt_p` with respect to the mirror line or plane
    :type mirror: Line
    :param pt_p: The parameter `pt_p` represents a point in a projective plane
    :type pt_p: Point
    :return: a ProjectivePlane<Point> object.
    """
    line_po = pt_p.circ(origin)
    pt_b = line_po.circ(mirror)
    return harm_conj(origin, pt_b, pt_p)