/*
    win32-compat.h -- header file with definitions for compatability with Windows.

    fgen -- Library for optimization using a genetic algorithm or particle swarm optimization.
    Copyright 2012-13, Harm Hanemaaijer <fgenfb at yahoo.com>

    This file is part of fgen.

    fgen is free software: you can redistribute it and/or modify it
    under the terms of the GNU Lesser General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    fgen is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
    License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with fgen.  If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef __GNUC__

#define POSITIVE_INFINITY_DOUBLE DBL_MAX
#define NEGATIVE_INFINITY_DOUBLE DBL_MIN

#define isnan(x) _isnan(x)

#else

#define POSITIVE_INFINITY_DOUBLE INFINITY
#define NEGATIVE_INFINITY_DOUBLE (- INFINITY)

#endif

