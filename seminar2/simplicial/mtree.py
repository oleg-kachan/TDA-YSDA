# A data structure for approximate search
#
# Copyright (C) 2017 Simon Dobson
# 
# This file is part of simplicial, simplicial topology in Python.
#
# Simplicial is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Simplicial is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Simplicial. If not, see <http://www.gnu.org/licenses/gpl.html>.

import math
import heapq

class MTree(object):
    '''A data structure for approximate search.

    An M-tree is an associative structure that uses points in space
    as keys for values. The structure can be queried by providing a
    point in the space and looking for the nearest neighbouring
    points, or looking for all points within a given range. The
    dimensions of the key space, the co-ordinates of each dimension,
    and the notion ot distance used are all abstracted. The tree is
    self-balancing, creating intermediate "routing" nodes to structure
    the real data points and reduce the amount of search needed.

    The M-tree supports a dict-like interface as well as a procedural
    one.

    :param dim: dimension of the key space (defaults to 2)
    :param levelSize: optional tuning parameter for the sizes of routing nodes

    '''

    DEFAULT_LEVEL_SIZE = 2     #: Default size of a routing node

    
    def __init__( self, dim = 2, levelSize = None ):
        self._dimension = dim
        self._root = MTreeNode(None, self.origin())
        if levelSize is None:
            levelSize = self.DEFAULT_LEVEL_SIZE
        self._levelSize = levelSize
        

    # ----- Distance -----

    def dimension( self ):
        '''Return the dimensions of the key space.

        :returns: the dimensions'''
        return self._dimension
    
    def origin( self ):
        '''Return the origin of the key space. By default this is
        a list of zeros of the dimensionality of the space.

        :returns: the origin'''
        return [ 0.0 ] * self.dimension()
    
    def distance(self, p, q ):
        '''Compute the distance between two key points. This method
        may be overridden by sub-classes. The default assumes that the
        key space is Euclidian, keyed by numbers, and computes the
        usual straight-line distance between the points.

        Any re-definition of this method must be a valid metric:

        - distance(p, p) == 0
        - distance(p, q) == distance(q, p)
        - distance(p, r) <= distance(p, q) + distance(q, r)

        :param p: one point
        :param q: the other point
        :returns: the distance between them'''
        sumsq = 0.0
        for d in range(self.dimension()):
            sumsq = sumsq + math.pow(q[d] - p[d], 2)
        return math.sqrt(sumsq)


    # ----- Insertion and deletion -----

    def insert( self, o, v ):
        '''Insert value v at a point o.

        :param o: the key point in the key space
        :param v: the value'''
        return self._insert(self._root, o, v)

    def _insert(self, n, o, v ):
        if not n.isLeaf():
            # n is a routing node, compute all children whose distance
            # to o is less than n's radius
            nin = [ r for r in n.children() if self.distance(r, o) < n.radius() ]
            if len(nin) > 0:
                # find the closest of these children to o
                ds = [ (r, self.distance(r, o)) for r in nin.children() ]
                sds = sorted(ds, key = (lambda (a, b): b))
                (rmin, dmin) = sds[0]
            else:
                # no such child, find the child whose radius will expand
                # the least if we add o as its child
                ds = [ (r, self.distance(r, o) - r.radius()) for r in n.children() ]
                sds = sorted(ds, key = (lambda (a, b): b))
                (rmin, dmin) = sds[0]

            # recursively insert into chosen branch
            return self._insert(rmin, o, v)
            
        else:
            # n is a leaf
            if n.size() < self._levelSize:
                # child is not full, add the new node
                n.addChild(MTreeNode(v, o))
            else:
                self._split(n, v, o)

    def _split(self, n, v, o ):
        if not n.isRoot():
            # splitting a non-root node
            p = n.parent()
            

    # ----- Nearest-neighbour queries -----

    def withinRange( self, q, rq ):
        '''Return all values lying within a given range of the search point.

        :param q: the search point
        :param rq: the range
        :returns: a list of (point, value) pairs, in order of increasing range'''
        return self._rs(self._root, q, rq )

    def _rs( self, n, q, rq ):
        rs = []
        p = n.parent()
        if not n.isLeaf():
            for r in n.children():
                if self.distance(n, q) - self.distance(r, n) <= rq + r.distanceToParent():
                    dr = self.distance(r, q)
                    if dr <= rq + r.radius():
                        rs.append(self._rs(r.children(), q, rq))
        else:
            for j in n.children():
                if self.distance(p, q) - self.distance(j, p) < rq:
                    dj = self.distance(j, q)
                    if dj <= rq:
                        rs.append(j.value())
                        
        return rs


    # ----- k-nearest neighbour queries -----

    def nearestNeighbours( self, q, k ):
        '''Return the k nearest neighbours to the search point. There needn't
        be (and typically won't be) a node at the actual search point. There may
        not be k elements in the returned list if k is larger than the number
        of points in the structure.

        :param q: the search point
        :param k: the number of neighbours sought
        :returns: a list of points in increasing order of distance'''
        pass

    #def _chooseNode( self, pr ):
        
    # ----- Dict-like interface-----

    def __getitem__( self, o ):
        '''Return the value closest to the given key point. This is equivalent
        to a nearest-neighbour query :meth:`MTree.nearestNeighbour` looking
        for exactly one neighbour.

        :param o: the search point
        :returns: the nearest value to the search point'''
        return self.nearestNeighbours(o, 1)
    
    def __setitem__( self, o, v ):
        '''Insert the given value at the given key point. This is equivalent
        to :meth:`MTree.insert`.

        :param o: the key point
        :param v: the value'''
        return self.insert(o, v)
    

class MTreeNode(object):
    '''Internal structure for representing MTree nodes. Never exposed to clients.'''
    
    def __init__( self, v, pos ):
        self._value = v
        self._position = pos
        self._children = []
        self._radius = 0.0
        self._parent = None
        self._distance = 0.0


    # ----- Access -----

    def value( self ):
        return self._value
    
    def parent( self ):
        return self._parent
    
    def children( self ):
        return self._children

    def radius( self ):
        return self._radius

    def distanceToParent( self ):
        return self._distance

    def isLeaf( self ):
        return len(self.children()) == 0

    def size( self ):
        return len(self.children())

    def isRoot( self ):
        return (self.parent() is None)

    def position( self ):
        return self.position()
    
    def __getitem__( self, d ):
        '''Return the position of the node in dimension d.

        :param d: the dimension
        :returns: the position along that dimension'''
        return self._position[d]


    # ----- Adding and removing children -----

    def addChild( self, mn ):
        self._children.append(mn)
        mn._parent = self
        
