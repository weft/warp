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
    unsigned    cellnum; /**< cell number */
    unsigned    cellmat; /**< cell material */
    float       dx;	 /**< dx */
    float       dy;	 /**< dy */
    float       dz;	 /**< dz */
    float       theta; /**< polar angle */
    float       phi;   /**< azimuthal angle */
};

/**
 * \struct geom_data datadef.h
 * \brief contains parameters of a wgeometry
 * \details extrema arrays, location array, cell and material numbers, fissile flag
 */
struct geom_data{
    float       min[3]; /**< array of coordinate (x,y,z) minima */
    float       max[3]; /**< array of coordinate (x,y,z) maxima */
    float       loc[3]; /**< array of coordinate (x,y,z) locations */
    int         cellnum; /**< cell number */
    int         matnum;  /**< material number */
    int         is_fissile; /**< fissile flag */
};

/**
 * \struct source_point datadef.h
 * \brief contains the parameters of the neutron source point
 * \details 3D Cartesian coordinates, shifted coordinates, distance to nearest surface, total macroscopic
 *  cross section, boundary condition enforcement flag
 */
struct source_point{ 
    float x;    /**< x-coordinate */
    float y; 	/**< y-coordinate */
    float z;    /**< z-coordinate */
    float xhat; /**< shifted x-coordinate */ 
    float yhat; /**< shifted y-coordinate */
    float zhat; /**< shifted z-coordinate */
    float surf_dist; /**< distance to nearest surface */
    float macro_t;   /**< total macroscopic cross section */
    float norm[3];
    unsigned enforce_BC; /**< boundary condition enforcement flag */
    unsigned weight;
};

/**
 * \struct qnode datadef.h
 * \brief quaternary search node
 */
struct qnode{
    float  values[4]; /**< array of values */
    qnode* leaves[4]; /**< array of node leaves */
};

/**
 * \struct qnode_host datadef.h
 * \brief quaternary search node host container
 */
struct qnode_host{
    qnode* cuda_pointer; /**< CUDA pointer */
    qnode  node;	 /**< quaternary node */
};

/**
 * \struct hit_buffer datadef.h
 * \brief contains information for the hit buffer
 * \details structure that contains cell number, material number, and fissile flag for the hit 
 * buffer.
 */
struct hit_buffer {
    int     cell; /**< cell number */
    int     mat;  /**< material number */
    int     fiss; /**< fissile flag */
};

/**
 * \struct intersection_point datadef.h
 * \brief contains information pertinent to an intersection point
 * \details 3D cartesian coordinates of intersection point, distance to nearest surface, first cell 
 * potentially hit, continuity flag, hit buffer, and hit buffer index
 */
struct intersection_point {
    float       x; /**< x-coordinate */
    float       y; /**< y-coordinate */
    float       z; /**< z-coordinate */
    float       surf_dist;  /**< distance to nearest surface */
    int         cell_first; /**< number of first cell to potentially enter */
    int         cont;       /**< continuity flag */
    hit_buffer  hitbuff[10]; /**< hit buffer array */
    unsigned    buff_index;  /**< index in hit buffer array */
    float       norm[3];
};

/**
 * \struct material_def datadef.h
 * \brief contains information that defines a material
 * \details material ID, number, fissile flag, number of isotopes, density, isotope list, isotope 
 * fraction list
 */
struct material_def {
    unsigned    id;	      /**< material ID */
    unsigned    matnum;       /**< material number */
    unsigned    is_fissile;   /**< fissile flag */
    unsigned    num_isotopes; /**< number of isotopes */
    float       density;      /**< density [g/cc] */
    std::vector<std::string>   isotopes;     /**< isotope list */
    float    *  fractions;    /**< isotope fractions */
};

#endif
