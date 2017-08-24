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

from vecmath import *
from operator import itemgetter

class QuickHull:
    def __init__(self, points):
        self.points = points
        self.normals = []
        self.faces = []
        self.shared_edges = {}

    def add_shared_edge(self, edge, face):
        if edge not in self.shared_edges.keys():
            self.shared_edges[edge] = []
        self.shared_edges[edge].append(face)

    def find_hull(self, I, i_p, i_q, i_r, n, n_0, d, sign):
        if len(I) == 0:
            face = len(self.faces)
            self.normals.append(n)
            self.faces.append((i_p, i_q, i_r))
            self.add_shared_edge((i_p, i_q) if i_p < i_q else (i_q, i_p), face)
            self.add_shared_edge((i_q, i_r) if i_q < i_r else (i_r, i_q), face)
            self.add_shared_edge((i_p, i_r) if i_p < i_r else (i_r, i_p), face)
            return

        d_max = 0.0
        i_d_max = -1
        for i in I:
            d_p = sign * dist_point_plane(make_vector(self.points[i]), n_0, d)
            if d_p > d_max:
                i_d_max = i
                d_max = d_p

        p = make_vector(self.points[i_p])
        q = make_vector(self.points[i_q])
        r = make_vector(self.points[i_r])
        s = make_vector(self.points[i_d_max])

        n_f1 = th_face_normal(s, p, q, r)
        n_0_f1, sign_f1 = orient_normal(n_f1, s)
        d_f1 = s * n_0_f1

        n_f2 = th_face_normal(s, p, r, q)
        n_0_f2, sign_f2 = orient_normal(n_f2, s)
        d_f2 = s * n_0_f2

        n_f3 = th_face_normal(s, r, q, p)
        n_0_f3, sign_f3 = orient_normal(n_f3, s)
        d_f3 = s * n_0_f3

        I_f1 = []
        I_f2 = []
        I_f3 = []
        for i in I:
            if i == i_d_max:
                continue

            v = make_vector(self.points[i])
            d_p_f1 = sign_f1 * dist_point_plane(v, n_0_f1, d_f1)
            d_p_f2 = sign_f2 * dist_point_plane(v, n_0_f2, d_f2)
            d_p_f3 = sign_f3 * dist_point_plane(v, n_0_f3, d_f3)

            d_p_min = 0
            for d_p_f in sorted([(1, d_p_f1), (2, d_p_f2), (3, d_p_f3)], key=itemgetter(1)):
                if d_p_f[1] > 0.0:
                    d_p_min = d_p_f[0]
                    break

            if d_p_min == 1:
                I_f1.append(i)
            elif d_p_min == 2:
                I_f2.append(i)
            elif d_p_min == 3:
                I_f3.append(i)

        self.find_hull(I_f1, i_d_max, i_p, i_q, n_f1, n_0_f1, d_f1, sign_f1)
        self.find_hull(I_f2, i_d_max, i_p, i_r, n_f2, n_0_f2, d_f2, sign_f2)
        self.find_hull(I_f3, i_d_max, i_r, i_q, n_f3, n_0_f3, d_f3, sign_f3)

    def generate(self):
        x_min = 1.0
        x_max = -1.0
        i_x_min = -1
        i_x_max = -1
        for i, p in enumerate(self.points):
            if p[0] < x_min:
                i_x_min = i
                x_min = p[0]
            if p[0] > x_max:
                i_x_max = i
                x_max = p[0]

        q = make_vector(self.points[i_x_min])
        r = make_vector(self.points[i_x_max])

        d_max = 0.0
        i_d_max = -1
        for i, p in enumerate(self.points):
            if i in [i_x_min, i_x_max]:
                continue

            d_p = dist_point_line(make_vector(p), q, r)
            if d_p > d_max:
                i_d_max = i
                d_max = d_p

        s = make_vector(self.points[i_d_max])

        n = (r-q).cross(s-q).normalize()
        n_0, _ = orient_normal(n, q)
        d = q * n_0

        d_max_abs = 0.0
        i_d_max_abs = -1
        for i, p in enumerate(self.points):
            if i in [i_x_min, i_x_max, i_d_max]:
                continue

            d_p = dist_point_plane(make_vector(p), n_0, d)
            if abs(d_p) > abs(d_max_abs):
                i_d_max_abs = i
                d_max_abs = d_p

        t = make_vector(self.points[i_d_max_abs])

        n_f1 = n_0
        n_0_f1 = n_0
        d_f1 = d
        sign_f1 = 1.0
        if d_max_abs > 0.0:
            n_f1 = -n_f1
            sign_f1 = -1.0

        n_f2 = th_face_normal(t, q, r, s)
        n_0_f2, sign_f2 = orient_normal(n_f2, t)
        d_f2 = t * n_0_f2

        n_f3 = th_face_normal(s, q, t, r)
        n_0_f3, sign_f3 = orient_normal(n_f3, s)
        d_f3 = s * n_0_f3

        n_f4 = th_face_normal(r, s, t, q)
        n_0_f4, sign_f4 = orient_normal(n_f4, r)
        d_f4 = r * n_0_f4

        I_f1 = []
        I_f2 = []
        I_f3 = []
        I_f4 = []
        for i, p in enumerate(self.points):
            if i in [i_x_min, i_x_max, i_d_max, i_d_max_abs]:
                continue

            v = make_vector(p)
            d_p_f1 = sign_f1 * dist_point_plane(v, n_0_f1, d_f1)
            d_p_f2 = sign_f2 * dist_point_plane(v, n_0_f2, d_f2)
            d_p_f3 = sign_f3 * dist_point_plane(v, n_0_f3, d_f3)
            d_p_f4 = sign_f4 * dist_point_plane(v, n_0_f4, d_f4)

            d_p_min = 0
            for d_p_f in sorted([(1, d_p_f1), (2, d_p_f2), (3, d_p_f3), (4, d_p_f4)], key=itemgetter(1)):
                if d_p_f[1] > 0.0:
                    d_p_min = d_p_f[0]
                    break

            if d_p_min == 1:
                I_f1.append(i)
            elif d_p_min == 2:
                I_f2.append(i)
            elif d_p_min == 3:
                I_f3.append(i)
            elif d_p_min == 4:
                I_f4.append(i)

        self.find_hull(I_f1, i_x_min, i_x_max, i_d_max, n_f1, n_0_f1, d_f1, sign_f1)
        self.find_hull(I_f2, i_x_min, i_x_max, i_d_max_abs, n_f2, n_0_f2, d_f2, sign_f2)
        self.find_hull(I_f3, i_x_min, i_d_max, i_d_max_abs, n_f3, n_0_f3, d_f3, sign_f3)
        self.find_hull(I_f4, i_x_max, i_d_max, i_d_max_abs, n_f4, n_0_f4, d_f4, sign_f4)

    def export(self, name):
        with open(name + '.obj', 'w') as fp:
            for p in self.points:
                fp.write('v ' + str(p[0]*10.0) + ' ' + str(p[1]*10.0) + ' ' + str(p[2]*10.0) + '\n')
            for n in self.normals:
                fp.write('vn ' + str(n.x) + ' ' + str(n.y) + ' ' + str(n.z) + '\n')
            for i, f in enumerate(self.faces):
                fp.write('f ' + str(f[0]+1) + '//' + str(i+1) + ' ' + str(f[1]+1) + '//' + str(i+1) + ' ' + str(f[2]+1) + '//' + str(i+1) + '\n')
