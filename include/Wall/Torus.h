#ifndef TORUS_H // include guard
#define TORUS_H

#include "voro++.hh"
using namespace voro;


// This class creates a custom toroidal wall object that is centered on the
// origin and is aligned with the xy plane. It is derived from the pure virtual
// "wall" class. The "wall" class contains virtual functions for cutting the
// Voronoi cell in response to a wall, and for telling whether a given point is
// inside the wall or not. In this derived class, specific implementations of
// these functions are given.
class wall_torus : public wall {
	public:
    // The wall constructor initializes constants for the major and
    // minor axes of the torus. It also initializes the wall ID
	// number that is used when the plane cuts are made. This is
	// only tracked with the voronoicell_neighbor class and is
	// ignored otherwise. It can be omitted, and then an arbitrary
	// value of -99 is used.
	wall_torus(double imjr,double imnr,int iw_id=-99)
		: w_id(iw_id), mjr(imjr), mnr(imnr) {};

	// This returns true if a given vector is inside the torus, and
	// false if it is outside. For the current example, this
	// routine is not needed, but in general it would be, for use
	// with the point_inside() routine in the container class.
	bool point_inside(double x,double y,double z) {
		double temp=sqrt(x*x+y*y)-mjr;
		return temp*temp+z*z<mnr*mnr;
	}

	// This template takes a reference to a voronoicell or
	// voronoicell_neighbor object for a particle at a vector
	// (x,y,z), and makes a plane cut to to the object to account
	// for the toroidal wall
	template<class vc_class>
	inline bool cut_cell_base(vc_class &c,double x,double y,double z) {
		double orad=sqrt(x*x+y*y);
		double odis=orad-mjr;
		double ot=odis*odis+z*z;
 
		// Unless the particle is within 1% of the major
		// radius, then a plane cut is made
		if(ot>0.01*mnr) {
			ot=2*mnr/sqrt(ot)-2;
			z*=ot;
			odis*=ot/orad;
			x*=odis;
			y*=odis;
			return c.nplane(x,y,z,w_id);
		}
		return true;
	}

	// These virtual functions are called during the cell
	// computation in the container class. They call instances of
	// the template given above.
	bool cut_cell(voronoicell &c,double x, double y,double z) {return cut_cell_base(c,x,y,z);}
	bool cut_cell(voronoicell_neighbor &c,double x, 
				double y,double z) {return cut_cell_base(c,x,y,z);}
private:
	// The ID number associated with the wall
	const int w_id;
	// The major radius of the torus
	const double mjr;
	// The minor radius of the torus
	const double mnr;
};

#endif