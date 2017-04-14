#ifndef DATADEF_H
#define DATADEF_H
#include <vector> 
#include <string> 

/**
 * \struct wtransform datadef.h
 * \brief contains parameters of a wtransform, a transform used to create a new instance of a wgemetry object
 * \details cell number and material, transform coordinates and angles
 */
struct wtransform {
	unsigned	cellnum;	/**< cell number */
	unsigned	cellmat;	/**< cell material */
	int 		tally_index;/**< index of tally associated with this instance*/
	float		dx;			/**< displacement in x */
	float		dy;			/**< displacement in y */
	float		dz;			/**< displacement in z */
	float		theta;		/**< roation in polar angle */
	float		phi;		/**< roation in azimuthal angle */
};

/**
 * \struct geom_data datadef.h
 * \brief contains parameters of a wgeometry
 * \details extrema arrays, location array, cell and material numbers, tally number, fissile flag
 */
struct geom_data{
	float	min[3];		/**< array of coordinate (x,y,z) minima */
	float	max[3];		/**< array of coordinate (x,y,z) maxima */
	float	loc[3];		/**< array of coordinate (x,y,z) locations */
	int		cellnum;	/**< cell number */
	int 	talnum;		/**< tally index*/
	int		matnum;		/**< material number */
	int		is_fissile;	/**< fissile flag */
};

/**
 * \struct spatial_data datadef.h
 * \brief contains the spatial parameters of the neutron
 * \details 3D Cartesian coordinates, direction vector, distance to nearest surface, total macroscopic
 *  cross section, surface normal vector, boundary condition enforcement flag, neutron statsitical weight
 */
struct spatial_data{ 
	float		x;			/**< x-coordinate */
	float		y;			/**< y-coordinate */
	float		z;			/**< z-coordinate */
	float		xhat;		/**< direction vector x-component */ 
	float		yhat;		/**< direction vector y-component */
	float		zhat;		/**< direction vector z-component */
	float		surf_dist;	/**< distance to nearest surface */
	float		norm[3];	/**< normal vector of surface intersection */
	unsigned	enforce_BC;	/**< boundary condition enforcement flag */
	unsigned	weight;		/**< particle statistical weight */
	unsigned	cell[10];
	float		dist[10];
	unsigned	mat[10];
};

/**
 * \struct intersection_point datadef.h
 * \brief contains information pertinent to an intersection point, used in OptiX
 * \details 3D cartesian coordinates of intersection point, distance to nearest surface, first cell 
 * potentially hit, material of hit cell, fissile flag of hit cell, normal vector of hit surface, 
 * hit cell sense value, index of tally for hit cell
 */
struct intersection_point {
	float		x; 			/**< x-coordinate */
	float		y; 			/**< y-coordinate */
	float		z; 			/**< z-coordinate */
	float		surf_dist;  /**< distance to nearest surface */
	int		cell; 		/**< cell number */
	int		mat;  		/**< material number */
	int		fiss; 		/**< fissile flag */
	float		norm[3];    /**< most recent normal */
	int		sense;   	/**< most recent cell sense */
	int 		tally_index;/**< tally index of most recent cell */
	int		buff_index;
	unsigned	launch_dex;
	int 		cont;
};

/**
 * \struct material_def datadef.h
 * \brief contains information that defines a material
 * \details material index, label number, fissile flag, number of isotopes, density, isotope list, isotope 
 * atom fraction list
 */
struct material_def {
	unsigned					id;				/**< material index */
	unsigned					matnum;			/**< material label number */
	unsigned					is_fissile;		/**< fissile flag */
	unsigned					num_isotopes;	/**< number of isotopes */
	float						density;		/**< density [g/cc] */
	std::vector<std::string>	isotopes;		/**< isotope list */
	float*  					fractions;		/**< isotope atom fractions */
};

/**
 * \struct dist_data datadef.h
 * \brief contains information that defines an ENDF cross section distribution
 * \details ernergy, length, law, interpolation type, variable/pdf/cdf arrays
 */
struct dist_data {
	float 			erg;	/**< energy point of this distribution */
	unsigned 		len;	/**< length of the arrays in this distribution */
	unsigned 		law;	/**< sampling law of this distribution */
	unsigned  		intt;	/**< interpolation flag of this distribution */
	float* 			var;	/**< independent variable  (mu, E, etc.) */
	float* 			pdf;	/**< probability density function array */
	float* 			cdf;	/**< cumulative density function array */
};

/**
 * \struct dist_container datadef.h
 * \brief    container for pointers that map the nearest distributions to the energy grid point where it resides
 * \details  pointers to the nearest distribution lower and greater in energy 
 */
struct dist_container {
	dist_data*		lower;			/**< pointer to distribution data of grid point below current energy */
	dist_data*		upper;			/**< pointer to distribution data of grid point above current energy */
};

/**
 * \struct cross_section_data datadef.h
 * \brief    structure that holds the topmost level of cross section data
 * \details  contains pointers and parameters to compute any reaction in this requested isotope set - number of isotopes, 
 * length of main energy grid, total number of reactions channels, reaction number vector, 
 * total reaction channels for each isotope, unionized energy grid vector, reaction Q values, cross section data, 
 * isotope atomic weight ratios, isotope temperatures, energy and scattering data distributions
 */
struct cross_section_data {
	unsigned			n_isotopes;					/**< number of isotopes */
	unsigned			energy_grid_len;			/**< length of main energy grid*/
	unsigned			total_reaction_channels;	/**< total number of reactions channels*/
	unsigned*			rxn_numbers;				/**< reaction number vector */
	unsigned*			rxn_numbers_total;			/**< total reaction channels for each isotope */
	float*				energy_grid;				/**< unionized energy grid vector */
	float*				Q;							/**< reaction Q values */
	float*				xs;							/**< cross section data matrix */
	float*				awr;						/**< isotope atomic weight ratio (AWR) list */
	float*				temp;						/**< isotope temperature list (MeV) */
	dist_container*		dist_scatter;				/**< scattering distribution data redirection matrix */
	dist_container*		dist_energy;				/**< energy distribution data redirection matrix */
};


/**
 * \struct particle_data datadef.h
 * \brief    structure that holds all the arrays that define a particle's state
 * \details  Hold arrays that define a neutron's state and/or need to be passed between kernels.  Data locality efficiency 
 * dictates that this must be a structure of arrays (SoA) and not be a structure built into arrays (array of structures - AoS).  
 * This container structure is passed to almost all kernels so they can access neutron state data.
 */
struct particle_data {
	spatial_data*	space;			/**< spatial data array */
	unsigned*		rxn;			/**< current reaction array */
	float*			E;				/**< energy array */
	float*			Q;				/**< current reaction Q value array */
	unsigned*		rn_bank;		/**< random number seed array */
	unsigned*		cellnum;		/**< current cell number array */
	unsigned*		matnum;			/**< current material number array */
	unsigned*		isonum;			/**< current isotope number array */
	int*			talnum;			/**< current tally number array */
	unsigned*		yield;			/**< total yield of history array */
	float*			weight;			/**< statistical weight array */
	unsigned*		index;			/**< current energy grid index array */
};

/**
 * \struct tally_data_host datadef.h
 * \brief    Tally data that lives on the host side.
 * \details  Tally data that lives on the host side.  Basically holds the same thing as the device array but also contains 64-bit
 * arrays that the tallies are accumulated into to avoid too much roundoff error.
 */
struct tally_data_host {
	float*			score;			/**< tally score */
	float*			square;			/**< tally square */
	unsigned*		count;			/**< tally count */
	double*			score_total;	/**< tally score accumulated total */
	double*			square_total;	/**< tally square accumulaed total */
	long unsigned*	count_total;	/**< tally count accumulated total */
	unsigned		cell;			/**< tally cell (input) */
	unsigned		length;			/**< tally length, edges are equi-log (input) */
	float			E_min;			/**< minimum energy (input) */
	float			E_max;			/**< maximum energy (input) */
};

/**
 * \struct tally_data datadef.h
 * \brief    Tally data that lives on the device side.
 * \details  Tally data that lives on the device side.  Everything needed to properly index a tally and the vectors to store 
 * the scores and keep track of statistics
 */
struct tally_data {
	float*			score;			/**< tally score */
	float*			square;			/**< tally square */
	unsigned*		count;			/**< tally count */
	unsigned		cell;			/**< tally cell (input) */
	unsigned		length;			/**< tally length, edges are equi-log (input) */
	float			E_min;			/**< minimum energy (input) */
	float			E_max;			/**< maximum energy (input) */
};

#endif
