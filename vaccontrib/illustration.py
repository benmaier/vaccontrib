"""
Helper classes and functions do illustrate contribution matrices
as segments of circles.
"""
# coding: utf-8


import matplotlib.pyplot as pl

from shapely.geometry import Point, LineString, Polygon, MultiPolygon
from shapely.affinity import translate

import numpy as np
from scipy.optimize import root

import bfmplot as bp
import matplotlib.pyplot as pl

_BASE_CIRCLE_RESOLUTION = 64

class CircleCaps():

    def __init__(self,r,h,w=1/10,circle_resolution=_BASE_CIRCLE_RESOLUTION):
        #self.y = sorted([amount0, amount1],key=lambda x:-x)
        self.r = r
        self.w = w
        self.h = h
        self.circle_resolution = circle_resolution
        self.compute()

    def compute(self):
        point = Point(0,0)
        self.circle = point.buffer(self.r,resolution=self.circle_resolution)
        r = self.r
        h = self.h
        w = self.w
        box0 = Polygon([(-2*r,h+w), (2*r,h+w), (2*r,h+w+2*r),(-2*r,h+w+2*r)])
        box1 = Polygon([(-2*r,h), (2*r,h), (2*r,h-2*r),(-2*r,h-2*r)])
        self.cap0 = self.circle.intersection(box0)
        self.cap1 = self.circle.intersection(box1)

        filtered_polygons = list(filter(lambda x: x.area > 0,  [self.cap0,self.cap1]))
        self.all = MultiPolygon(filtered_polygons)

    def get_areas(self):
        return (self.cap0.area, self.cap1.area)

    def area(self):
        return self.all.area

class CircleCapSegments():

    def __init__(self,r,h,x_hi,x_lo,w=1/10,circle_resolution=_BASE_CIRCLE_RESOLUTION):
        #self.y = sorted([amount0, amount1],key=lambda x:-x)
        self.r = r
        self.w = w
        self.h = h
        self.x_lo = x_lo
        self.x_hi = x_hi
        self.circle_resolution = circle_resolution
        self.compute()

    def compute(self):
        point = Point(0,0)
        self.circle = point.buffer(self.r,resolution=self.circle_resolution)
        r = self.r
        h = self.h
        w = self.w
        x_lo = self.x_lo
        x_hi = self.x_hi
        box0 = Polygon([(-2*r,h+w), (2*r,h+w), (2*r,h+w+2*r),(-2*r,h+w+2*r)])
        box1 = Polygon([(-2*r,h), (2*r,h), (2*r,h-2*r),(-2*r,h-2*r)])
        self.cap0 = self.circle.intersection(box0)
        self.cap1 = self.circle.intersection(box1)

        box_lo_left  = Polygon([(-2*r,h+w/2), (x_lo,h+w/2), (x_lo,h-2*r),(-2*r,h-2*r)])
        box_lo_right = Polygon([(x_lo + w,h+w/2), (2*r,h+w/2), (2*r,h-2*r),(x_lo+w,h-2*r)])

        box_hi_left  = Polygon([(-2*r,h+w/2), (x_hi,h+w/2), (x_hi,h+2*r),(-2*r,h+2*r)])
        box_hi_right = Polygon([(x_hi+w,h+w/2), (2*r,h+w/2), (2*r,h+2*r),(x_hi+w,h+2*r)])

        self.seg00 = self.cap0.intersection(box_hi_left)
        self.seg01 = self.cap0.intersection(box_hi_right)
        self.seg10 = self.cap1.intersection(box_lo_left)
        self.seg11 = self.cap1.intersection(box_lo_right)

        filtered_polygons = list(filter(lambda x: x.area > 0, [self.seg00, self.seg01, self.seg10, self.seg11]))
        self.all = MultiPolygon(filtered_polygons)

    def get_areas(self):
        return [ [self.seg00.area, self.seg01.area], [self.seg10.area, self.seg11.area] ]

    def area(self):
        return self.all.area


class CircleCapPresentation():

    def __init__(self,y,r=1,w=1/10,area=None,circle_resolution=_BASE_CIRCLE_RESOLUTION):

        if area is not None:
            self.initial_r = r = np.sqrt(area/np.pi)
            self.target_area = area
        else:
            self.initial_r = r
            self.target_area = np.pi * r**2

        self.y = np.array(y)
        assert(self.y.shape == (2,))
        self.relative_y = self.y/self.y.sum()
        self.target_areas = self.target_area * self.relative_y
        self.w = w
        self.circle_resolution = circle_resolution

    def get_caps(self,r,h,w):
        caps = CircleCaps(r,h,w,circle_resolution=self.circle_resolution)
        return caps

    def get_areas(self,r,h,w):
        caps = self.get_caps(r,h,w)
        return np.array(caps.get_areas())

    def compute(self,tol=1e-3):
        r0 = self.initial_r
        h0 = (self.relative_y[0] - self.relative_y[1])
        #print(r0,h0)
        func = lambda param: self.get_areas(param[0], param[1], self.w) - self.target_areas
        solution = root(func,[r0,h0],tol=tol)
        self.caps = self.get_caps(solution.x[0],solution.x[1],self.w)
        self.r = solution.x[0]

        return self
        #print(solution)

class CircleCapPresentationConstR():

    def __init__(self,y,r=1,w=1/10,circle_resolution=_BASE_CIRCLE_RESOLUTION):

        self.r = r
        self.target_area = np.pi * r**2

        self.y = np.array(y)
        assert(self.y.shape == (2,))
        self.relative_y = self.y/self.y.sum()
        self.rel_target_areas = self.relative_y
        self.w = w
        self.circle_resolution = circle_resolution


    def get_caps(self,r,h,w):
        caps = CircleCaps(r,h,w,circle_resolution=self.circle_resolution)
        return caps

    def get_relative_areas(self,r,h,w):
        caps = self.get_caps(r,h,w)
        areas = np.array(caps.get_areas())
        return areas / areas.sum()

    def compute(self,tol=1e-3):
        h0 = 0

        def func(param):
            rel = self.get_relative_areas(self.r, param[0], self.w)
            trg = self.rel_target_areas
            return [ rel[0] - trg[0] ]

        solution = root(func,[h0],tol=tol)
        self.caps = self.get_caps(self.r,solution.x[0],self.w)
        return self

class CircleCapSegmentPresentation():

    def __init__(self,C,r=1,w=1/10,area=None,circle_resolution=_BASE_CIRCLE_RESOLUTION):

        if area is not None:
            self.initial_r = r = np.sqrt(area/np.pi)
            self.target_area = area
        else:
            self.initial_r = r
            self.target_area = np.pi * r**2

        self.C = np.array(C)
        assert(self.C.shape == (2,2))
        self.relative_C = self.C/self.C.sum()
        self.target_areas = (self.target_area * self.relative_C)
        self.w = w
        self.circle_resolution = circle_resolution


    def get_segs(self,r,h,xhi,xlo,w):
        segs = CircleCapSegments(r,h,xhi,xlo,w,circle_resolution=self.circle_resolution)
        return segs

    def get_areas(self,r,h,xhi,xlo,w):
        segs = self.get_segs(r,h,xhi,xlo,w)
        return np.array(segs.get_areas())

    def compute(self,tol=1e-3):
        r0 = self.initial_r

        h0 = 0
        xhi0 = 0
        xlo0 = 0

        func = lambda p: self.get_areas(p[0], p[1], p[2], p[3], self.w).ravel() - self.target_areas.ravel()
        solution = root(func,[r0,h0,xhi0,xlo0],tol=tol)
        p = solution.x.tolist() + [self.w]
        self.r = p[0]
        self.segs = self.get_segs(*p)

        return self

class CircleCapSegmentPresentationConstR():

    def __init__(self,C,r=1,w=1/10,circle_resolution=_BASE_CIRCLE_RESOLUTION):

        self.r = r
        self.target_area = np.pi * r**2

        self.C = np.array(C)
        assert(self.C.shape == (2,2))
        self.relative_C = self.C/self.C.sum()
        self.rel_target_areas = self.relative_C
        self.w = w
        self.circle_resolution = circle_resolution

    def get_segs(self,r,h,xhi,xlo,w):
        segs = CircleCapSegments(r,h,xhi,xlo,w,circle_resolution=self.circle_resolution)
        return segs

    def get_relative_areas(self,r,h,xhi,xlo,w):
        segs = self.get_segs(r,h,xhi,xlo,w)
        areas = np.array(segs.get_areas())
        return areas / areas.sum()

    def compute(self,tol=1e-3):
        r0 = self.r

        h0 = 0
        xhi0 = 0
        xlo0 = 0

        def func(p):
            areas = self.get_relative_areas(self.r, p[0], p[1], p[2], self.w).ravel()
            targets = self.rel_target_areas.ravel()
            return areas[1:] - targets[1:]

        solution = root(func,[h0,xhi0,xlo0],tol=tol)

        p = [self.r] + solution.x.tolist() + [self.w]
        self.segs = self.get_segs(*p)

        return self

class JoinedVectorAndMatrixPresentation():

    def __init__(self,
                 vector_presentation,
                 matrix_presentation,
                 translate_x=None,
                ):

        caps = vector_presentation.caps
        segs = matrix_presentation.segs

        cap_width = 2 * vector_presentation.r
        seg_r = matrix_presentation.r

        if translate_x is None:
            translate_x = cap_width + seg_r

        self.cap0 = caps.cap0
        self.cap1 = caps.cap1

        self.seg00 = translate(segs.seg00, translate_x)
        self.seg01 = translate(segs.seg01, translate_x)
        self.seg10 = translate(segs.seg10, translate_x)
        self.seg11 = translate(segs.seg11, translate_x)

        self.caps = [self.cap0, self.cap1]
        self.segs = [self.seg00, self.seg01, self.seg10, self.seg11]
        self.seg_matrix = [ [self.seg00, self.seg01], [self.seg10, self.seg11] ]

        filtered_polygons = list(filter(lambda x: x.area > 0, self.caps+self.segs))
        self.all = MultiPolygon(filtered_polygons)

    def plot(self,ax=None,upper_color=None,lower_color=None,brighter_base=2,ec=None):
        if upper_color is None:
            upper_color = bp.epipack[1]
        if lower_color is None:
            lower_color = bp.epipack[2]
        upper_brighter = bp.brighter(upper_color, brighter_base)
        lower_brighter = bp.brighter(lower_color, brighter_base)

        if ax is None:
            fig, ax = pl.subplots(1,1)

        def _plot(geom, color, ec=None):
            xs, ys = geom.exterior.xy
            ax.fill(xs, ys, fc=color, ec=ec)

        for geom, color in [
                (self.cap0, upper_color),
                (self.cap1, lower_color),
                (self.seg00, upper_color),
                (self.seg01, upper_brighter),
                (self.seg10, lower_color),
                (self.seg11, lower_brighter),
            ]:
            if geom.area > 0:
                _plot(geom, color, ec=ec)

        ax.axis('equal')

        return ax

def get_circular_vector_and_matrix_presentation(y, C, r=1, w=0.1):

    ypres = CircleCapPresentationConstR(y,r=r,w=w).compute()
    Cpres = CircleCapSegmentPresentation(C, area=ypres.caps.area()*C.sum(), w=w).compute()
    joined = JoinedVectorAndMatrixPresentation(ypres, Cpres)

    return joined

if __name__ == "__main__":

    import vaccontrib as vc

    C = vc.covid.get_reduced_vaccinated_susceptible_contribution_matrix_covid([1.,6,6,6,6],variant='delta')
    K = vc.covid.get_next_generation_matrix_covid([1.,6,6,6,6],variant='delta')

    C = vc.covid.get_reduced_vaccinated_susceptible_contribution_matrix_covid([1.,4,4,4,4],variant='alpha')
    K = vc.covid.get_next_generation_matrix_covid([1.,4,4,4,4],variant='alpha')

    C = vc.covid.get_reduced_vaccinated_susceptible_contribution_matrix_covid([1.,1,1,1,1],variant='delta')
    K = vc.covid.get_next_generation_matrix_covid([1.,1,1,1,1],variant='delta')

    y = vc.get_eigenvector(K)
    y = y.sum(axis=0)
    y = np.array([y[0],y[1:].sum()])

    pres = get_circular_vector_and_matrix_presentation(y, C)

    pres.plot()

    pl.show()
