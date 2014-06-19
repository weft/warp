#ifndef DATADEF_H
#define DATADEF_H

/**
 * \struct wtransform datadef.h
 *
 * contains parameters of a wtransform.
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
 *
 * contains parameters of a wgeometry.
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
 *
 * contains the parameters of the neutron source point
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
    unsigned enforce_BC; /**< boundary condition enforcement flag */
};

/**
 * \struct qnode datadef.h
 *
 * quaternary search node
 */
struct qnode{
    float  values[4]; /**< array of values */
    qnode* leaves[4]; /**< array of node leaves */
};

/**
 * \struct qnode_host datadef.h
 *
 * quaternary search node host container
 */
struct qnode_host{
    qnode* cuda_pointer; /**< CUDA pointer */
    qnode  node;	 /**< quaternary node */
};

/**
 * \struct hit_buffer datadef.h
 *
 * contains information for the hit buffer.
 */
struct hit_buffer {
    int     cell; /**< cell number */
    int     mat;  /**< material number */
    int     fiss; /**< fissile flag */
};

/**
 * \struct instersection_point datadef.h
 *
 * contains information pertinent to an intersection point
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
};

/**
 * \struct material_def datadef.h
 *
 * contains information that defines a material
 */
struct material_def {
    unsigned    id;	      /**< material ID */
    unsigned    matnum;       /**< material number */
    unsigned    is_fissile;   /**< fissile flag */
    unsigned    num_isotopes; /**< number of isotopes */
    float       density;      /**< density [g/cc] */
    unsigned *  isotopes;     /**< isotope list */
    float    *  fractions;    /**< isotope fractions */
};

#endif
