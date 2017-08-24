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

from math import sqrt

class Vector:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def __sub__(self, other):
        return Vector(self.x-other.x, self.y-other.y, self.z-other.z)

    def __mul__(self, other):
        return self.x*other.x + self.y*other.y + self.z*other.z

    def __neg__(self):
        return Vector(-self.x, -self.y, -self.z)

    def cross(self, other):
        return Vector(self.y*other.z-self.z*other.y, self.z*other.x-self.x*other.z, self.x*other.y-self.y*other.x)

    def length(self):
        return sqrt(self*self)

    def normalize(self):
        l = self.length()
        return Vector(self.x/l, self.y/l, self.z/l)

def make_vector(p):
    return Vector(p[0], p[1], p[2])

def dist_point_line(p, q, r):
    u = r - q
    return (p-q).cross(u).length() / u.length()

def dist_point_plane(q, n_0, d):
    return q * n_0 - d

def th_face_normal(p, q, r, v):
    n = (q-p).cross(r-p).normalize()
    return n if (v-p)*n < 0.0 else -n

def orient_normal(n, p):
    return (n, 1.0) if p*n >= 0 else (-n, -1.0)
