from icebin import element_l1, ibgrid
import unittest

def area_of_polygon(polygon):
    """The Surveyor's Formula"""
    ret = 0.
    v0 = polygon[-1]
    for v1 in polygon:
        ret += v0.x*v1.y - v1.x*v0.y
        v0 = v1

    ret *= .5
    return ret


class ElementTests(unittest.TestCase):
    def setUp(self):
        self.triangle = [
            ibgrid.Vertex(0, 0., 0.),
            ibgrid.Vertex(1, 1., 0.),
            ibgrid.Vertex(2, 0., 1.)]

    def test_eqn_plane(self):
        triangle = self.triangle
        for basis_vertex in range(0,len(triangle)):
            A,B,C = element_l1.eqn_plane_subelement(triangle, basis_vertex)

            for j in range(0,len(triangle)):
                v = triangle[j]
                self.assertAlmostEqual(
                    1. if j == basis_vertex else 0.,
                    A*v.x + B*v.y + C)

    def test_surveyors_formula(self):
        """Check that our Surveyor's Formula above is correct.
        Evaluate the area of some simple polygons."""
        triangle = [
            ibgrid.Vertex(0, 0., 0.),
            ibgrid.Vertex(1, 1., 0.),
            ibgrid.Vertex(2, 0., 1.)]
        self.assertAlmostEqual(.5, area_of_polygon(triangle))

        square = [
            ibgrid.Vertex(0, 0., 0.),
            ibgrid.Vertex(1, 1., 0.),
            ibgrid.Vertex(2, 1., 1.),
            ibgrid.Vertex(2, 0., 1.)]
        self.assertAlmostEqual(1., area_of_polygon(square))


    def test_integrate_subelement(self):
        triangle = self.triangle

        triangle_area = area_of_polygon(triangle)

        # ------ Integrate the sub-element that's 1 at element[0]
        for basis_vertex in range(0, len(triangle)):
            e0 = [triangle[0], triangle[1], triangle[2]]    # Element to integrate
            poly0 = [triangle[1], triangle[2], triangle[0]]
            i0 = element_l1.integrate_subelement(e0, basis_vertex, poly0)
            self.assertAlmostEqual((1./3.)*triangle_area, i0)


    def test_cutoff_element(self):
        element = [
            ibgrid.Vertex(0, -1., 0.),
            ibgrid.Vertex(1, 1., 0.),
            ibgrid.Vertex(2, 0., 1.)]

        # Cut in half
        poly0 = [
            ibgrid.Vertex(0, 0., 0.),
            ibgrid.Vertex(1, 1., 0.),
            ibgrid.Vertex(2, 0., 1.)]
        i0 = element_l1.integrate_subelement(element, 2, poly0)
        self.assertAlmostEqual(.5*(1./3.)*area_of_polygon(element), i0)

        # Cut off just the tip
        poly1 = [
            ibgrid.Vertex(0, -.5, .5),
            ibgrid.Vertex(1, .5, .5),
            ibgrid.Vertex(2, 0., 1.)]
        i1 = element_l1.integrate_subelement(element, 2, poly1)

        area_poly1 = area_of_polygon(poly1)
        vol1_part0 = .5*(1./3.)*area_poly1 # Volume of the tip pyramid, starting at .5 height 
        vol1 = vol1_part0 \
            + (area_poly1 * .5)         # Volume on the "pedastal" beneath
        self.assertAlmostEqual(vol1, i1)

        # Cut off just the base
        poly2 = [
            ibgrid.Vertex(0, -1., 0.),
            ibgrid.Vertex(1, 1., 0.),
            ibgrid.Vertex(1, .5, .5),
            ibgrid.Vertex(0, -.5, .5)]
        i2 = element_l1.integrate_subelement(element, 2, poly1)
        self.assertAlmostEqual(i0*2.,i1+i2)

