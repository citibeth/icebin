/*
 * GLINT2: A Coupling Library for Ice Models and GCMs
 * Copyright (c) 2013 by Robert Fischer
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
C++ API for proj.4 Projection Library.
By Robert Fischer: robert.fischer@nasa.gov
April 5, 2012

This file is in the public domain.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 THE AUTHORS OR AUTHOR'S EMPLOYERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 DEALINGS IN THE SOFTWARE.
*/

#ifndef PROJPP_HPP
#define PROJPP_HPP

#include <proj_api.h>
#include <string>

#if 0
class ProjContext;

class ProjStatic {
public:
    ProjContext const defaultContext;

    ProjStatic();
};
extern ProjStatic projStatic;


class ProjContext {
    projCtx ctx;

public :


};

#endif

namespace giss {

/** Memory-safe peer class for projections in the proj.4 C library.
@see http://trac.osgeo.org/proj */
class Proj {
    projPJ pj;

    explicit Proj(projPJ _pj) : pj(_pj) {}

public:
    Proj() : pj(0) {}

	bool is_valid() const { return (pj != 0); }

    friend int transform(Proj const &src, Proj const &dest,
        long point_count, int point_offset,
        double *x, double *y, double *z);


    // ------------------ Five Standard constructors/methods for C++
    // See: http://www2.research.att.com/~bs/C++0xFAQ.html

	/** Create a projection.
	@param definition The proj.4 string describing the projection. */
    explicit Proj(std::string const &definition)
    {
        pj = pj_init_plus(definition.c_str());
        // pj_def = 0;
    }

	/** Create a projection.
	@param definition The proj.4 string describing the projection. */
    explicit Proj(char const *definition)
    {
        pj = pj_init_plus(definition);
    }

	void clear() {
		if (pj) pj_free(pj);
		pj = 0;
	}


    ~Proj()
    {
        if (pj) {
			pj_free(pj);
			pj = 0;
		}
        // if (pj_def) pj_dalloc(pj_def);
    }

    /** Transfer ownership (move) */
    Proj(Proj&& h) : pj{h.pj} //, pj_def{h.pj_def}
    {
        h.pj = 0;
        // h.pj_def = 0;
    }

    /** Transfer value */
    Proj& operator=(Proj&& h)
    {
        if (pj) pj_free(pj);
        pj = h.pj;
        h.pj = 0;
		return *this;
    }

    /** Copy constructor */
    Proj(const Proj &h)
    {
        char *pj_def = pj_get_def(h.pj, 0);
        pj = pj_init_plus(pj_def);
        pj_dalloc(pj_def);
    }

    /** Copying of Proj not allowed.
	No copy with operator=() */
    Proj& operator=(const Proj&) = delete;

    // --------------------------- Other Stuff


    /** Returns TRUE if the passed coordinate system is geographic
    (proj=latlong). */
    int is_latlong() const
        { return pj_is_latlong(pj); }


    /** Returns TRUE if the coordinate system is geocentric (proj=geocent). */
    int is_geocent() const
        { return pj_is_geocent(pj); }

    /** Returns the PROJ.4 initialization string suitable for use with
    pj_init_plus() that would produce this coordinate system, but with the
    definition expanded as much as possible (for instance +init= and
    +datum= definitions).
    @param options Unused at this point
    */
    std::string get_def(int options=0) const
    {
        char *pj_def = 0;
        pj_def = pj_get_def(pj, options);

        std::string ret = std::string(pj_def);
        pj_dalloc(pj_def);
        return ret;
    }


    /** Returns a new coordinate system definition which is the geographic
    coordinate (lat/long) system underlying pj_in.  This is essential in
	creating TWO Proj objects to be used in the transform() subroutines below. */
    Proj latlong_from_proj() const
    {
        return Proj(pj_latlong_from_proj(pj));
    }

};


/** Transforms a set of coordinate pairs.
@param src Source coordinate system
@param dest Destination coordinate system.
@param point_count The number of points to be processed (the size of the x/y/z arrays).
@param point_offset The step size from value to value (measured in doubles) within the x/y/z arrays - normally 1 for a packed array. May be used to operate on xyz interleaved point arrays.
@param x IN/OUT: The x (or longitude) array of points.
@param y IN/OUT: The y (or longitude) array of points.
@param z IN/OUT, OPTIONAL: The z (or longitude) array of points.
*/
inline int transform(Proj const &src, Proj const &dest,
    long point_count, int point_offset, double *x, double *y, double *z=0)
{
    return pj_transform(src.pj, dest.pj,
        point_count, point_offset, x, y, z);
}

/** Transforms a single coordinate pair
@param src Source coordinate system
@param dest Destination coordinate system.
@param x0 Source x (or longitude) coordinate (radians)
@param y0 Source y (or latitude) coordinate (radians)
@param x1 Destination x (or longitude) coordinate (radians)
@param y1 Destination y (or latitude) coordinate (radians) */
inline int transform(Proj const &src, Proj const &dest,
    double x0, double y0, double &x1, double &y1)
{
    x1 = x0;
    y1 = y0;
    int ret = transform(src, dest, 1, 1, &x1, &y1);
    return ret;
}


}

#endif
