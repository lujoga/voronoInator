# Copyright (C) 2017  Luca S.
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

from random import random
from quickhull import QuickHull
#import cairocffi as cairo

#MM_TO_DOTS = 72 / 25.4

def voronoi(w, h, n):
    P = [ (random(), random()) for i in range(n) ]
    CH = QuickHull([ (p[0], p[1], p[0]**2 + p[1]**2) for p in P ])
    CH.generate()
    CH.export('convexhull')

#   svg = cairo.SVGSurface('voronoi.svg', w * MM_TO_DOTS, h * MM_TO_DOTS)
#   ctx = cairo.Context(svg)

if __name__ == '__main__':
    voronoi(0, 0, 20)
