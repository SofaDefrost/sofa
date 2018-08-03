# -*- coding: ASCII -*-
from libcpp cimport bool
from cpp_vector cimport Vec3d as _Vec3d
from cpp_Ray cimport Ray as _Ray

cdef extern from "SofaGeometry/Plane.h" namespace "sofageometry":
    cdef cppclass Plane:
        Plane(const _Vec3d& normal, const double& distance) except +
        Plane(const _Vec3d& normal, const _Vec3d& point) except +

        _Vec3d normal
        double distance

        bool raycast(_Ray& r, double& p) 

