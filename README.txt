======================================
Matching catalogues using Groth method
======================================

A Pattern-Matching Algorithm for Two-Dimensional Coordinate Lists
Eduard J. Groth 1986 AJ 19, 5

http://adsabs.harvard.edu/abs/1986AJ.....91.1244G

Gmatch is released under the GPL (version 3 or any later version)


Sample usage
------------

    >>> from gmatch import gmatch
    >>> # load cat1 and cat2 from somewhere
    >>> # they are numpy arrays
    >>> # N rows, 2 columns
    >>> cat1.shape
    (20, 2)
    >>> cat2.shape
    (20, 2)
    >>> matches = gmatch(cat1, cat2, eps=1e-3)
    >>> # eps is the *relative* tolerance in xy
    >>> matches
    None
    >>> # if the return value is None, the catalogues
    >>> # do not match


    >>> from gmatch import gmatch
    >>> # load cat1 and cat2 from somewhere
    >>> # they are numpy arrays
    >>> matches = gmatch(cat1, cat2, eps=1e-3)
    >>> # eps is the *relative* tolerance in xy
    >>> if matches:
    ...     matches[0]
    ...     matches[1]
    >>> # if the return value is not None
    >>> # then is a tuple
    >>> # first element matched objects in cat1
    >>> # second element matched objects in cat2


