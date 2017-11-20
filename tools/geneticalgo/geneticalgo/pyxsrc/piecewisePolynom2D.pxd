#!/usr/bin/env python
# -*- coding: utf-8 -*-
####################################################################################################
## Copyright 2017 INRIA
##
## This file is part of the ShapeGenerator project.
##
## Contributors:
##     - damien.marchal@inria.fr
##     - thomas.morzadec@inria.fr
##
####################################################################################################

# distutils: language=c++
# cython: profile=True


cimport primitives2D

cdef class Tangency2D(object):

    cdef public float firstCoord, secondCoord

cpdef Det(u , v)

cdef class WeightedTangencedPoint2D(object):

    cdef public primitives2D.Point2D p
    cdef public Tangency2D t
    cdef public double wLeft, wRight

cdef class Polynom(primitives2D.Shape2D):
    cdef public str side
    cdef public double orientation
    cdef public WeightedTangencedPoint2D X1
    cdef public WeightedTangencedPoint2D X2
    cdef public primitives2D.Vector2D vect
    cdef public list identifier
