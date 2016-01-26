#ifndef DATADEF_H
#define DATADEF_H
#include <vector> 
#include <string> 

/**
 * \struct wtransform datadef.h
 * \brief contains parameters of a wtransform
 * \details cell number and material, transform coordinates and angles
 */
struct wtransform {
	unsigned	cellnum;	/**< cell number */
	unsigned	cellmat;	/**< cell material */
	int 		tally_index;/**< tally index*/
	float		dx;			/**< dx */
	float		dy;			/**< dy */
	float		dz;			/**< dz */
	float		theta;		/**< polar angle */
	float		phi;		/**< azimuthal angle */
};

/**
 * \struct geom_data datadef.h
 * \brief contains parameters of a wgeometry
 * \details extrema arrays, location array, cell and material numbers, fissile flag
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
 * \brief contains the parameters of the neutron source point
 * \details 3D Cartesian coordinates, shifted coordinates, distance to nearest surface, total macroscopic
 *  cross section, boundary condition enforcement flag
 */
struct spatial_data{ 
	float		x;			/**< x-coordinate */
	float		y;			/**< y-coordinate */
	float		z;			/**< z-coordinate */
	float		xhat;		/**< shifted x-coordinate */ 
	float		yhat;		/**< shifted y-coordinate */
	float		zhat;		/**< shifted z-coordinate */
	float		surf_dist;	/**< distance to nearest surface */
//	float		macro_t;	/**< total macroscopic cross section */
	float		norm[3];	/**< normal vector of intersection */
	unsigned	enforce_BC;	/**< boundary condition enforcement flag */
	unsigned	weight;		/**< particle statistical weight */
};

/**
 * \struct intersection_point datadef.h
 * \brief contains information pertinent to an intersection point
 * \details 3D cartesian coordinates of intersection point, distance to nearest surface, first cell 
 * potentially hit, continuity flag, hit buffer, and hit buffer index
 */
struct intersection_point {
	float		x; 			/**< x-coordinate */
	float		y; 			/**< y-coordinate */
	float		z; 			/**< z-coordinate */
	float		surf_dist;  /**< distance to nearest surface */
	int			cell; 		/**< cell number */
	int			mat;  		/**< material number */
	int			fiss; 		/**< fissile flag */
	float		norm[3];    /**< most recent normal */
	int			sense;   	/**< most recent cell sense */
	int 		tally_index;/**< tally index of most recent cell */
};

/**
 * \struct material_def datadef.h
 * \brief contains information that defines a material
 * \details material ID, number, fissile flag, number of isotopes, density, isotope list, isotope 
 * fraction list
 */
struct material_def {
	unsigned					id;				/**< material ID */
	unsigned					matnum;			/**< material number */
	unsigned					is_fissile;		/**< fissile flag */
	unsigned					num_isotopes;	/**< number of isotopes */
	float						density;		/**< density [g/cc] */
	std::vector<std::string>	isotopes;		/**< isotope list */
	float*  					fractions;		/**< isotope fractions */
};

/**
 * \struct dist_data datadef.h
 * \brief contains information that defines a material
 * \details material ID, number, fissile flag, number of isotopes, density, isotope list, isotope 
 * fraction list
 */
struct dist_data {
	float 			erg;	/**< energy point of this distribution */
	unsigned 		len;	/**< length of the arrays in this ditribution */
	unsigned 		law;	/**< sampling law of this distribution */
	unsigned  		intt;	/**< interpolation flag of this distribution */
	float* 			var;	/**< independent variable  (mu, E, etc.) */
	float* 			pdf;	/**< probability density function array */
	float* 			cdf;	/**< cumulative density function array */
};

/**
 * \struct dist_container datadef.h
 * \brief contains information that defines a material
 * \details material ID, number, fissile flag, number of isotopes, density, isotope list, isotope 
 * fraction list
 */
struct dist_container {
	dist_data*		lower;			/**< pointer to distribution data of grid point below current energy */
	dist_data*		upper;			/**< pointer to distribution data of grid point above current energy */
};

/**
 * \struct cross_section_data datadef.h
 * \brief contains information that defines a material
 * \details material ID, number, fissile flag, number of isotopes, density, isotope list, isotope 
 * fraction list
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
	dist_container*		dist_scatter;				/**< distribution data redirection matrix */
	dist_container*		dist_energy;				/**< distribution data redirection matrix */
};


/**
 * \struct particle_data datadef.h
 * \brief contains information that defines a material
 * \details material ID, number, fissile flag, number of isotopes, density, isotope list, isotope 
 * fraction list
 */
struct particle_data {
	spatial_data*	space;			/**< distribution data redirection matrix */
	unsigned*		rxn;			/**< distribution data redirection matrix */
	float*			E;				/**< distribution data redirection matrix */
	float*			Q;				/**< distribution data redirection matrix */
	unsigned*		rn_bank;		/**< distribution data redirection matrix */
	unsigned*		cellnum;		/**< distribution data redirection matrix */
	unsigned*		matnum;			/**< distribution data redirection matrix */
	unsigned*		isonum;			/**< distribution data redirection matrix */
	int*			talnum;
	unsigned*		yield;			/**< distribution data redirection matrix */
	float*			weight;			/**< distribution data redirection matrix */
	unsigned*		index;			/**< distribution data redirection matrix */
};

/**
 * \struct tally_data datadef.h
 * \brief contains information that defines a material
 * \details material ID, number, fissile flag, number of isotopes, density, isotope list, isotope 
 * fraction list
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
