#ifndef OPTIX_STUFF_H
#define OPTIX_STUFF_H

#include <optix_world.h>

/**
 * \brief OptiX stuff class
 */

class optix_stuff{
	optix::Context 	context; /**< OptiX context */
	std::string accel_type, traverse_type, image_type; /**< types of acceleration, traverse, image */
	unsigned mincell; /**< minimum (usually innermost) cell */
	unsigned maxcell; /**< maximum (usually outermost) cell */
	unsigned outer_cell; /**< outermost cell */
	unsigned n_materials; /**< number of materials */
	unsigned compute_device; /**< compute device number (always zero) */
	unsigned GEOM_FLAG; /**< geometry flag */
	/**
	 * \brief makes a geometry with a transform
	 * \details makes the top level group/acceleration as children of the top level object. for 
	 * each primitive in the geometry, creates the geometry type, sets the intersection and 
	 * bounding box programs, sets the hit programs to the geometry material, sets the program 
	 * variables for the instance, creates the instances, sets cell-specific variables, makes the 
	 * geometry group for the primitive, puts the geometry instance into its corresponding group, 
	 * makes any necessary transforms, and attaches to the root node.
	 * @param[in] problem_geom - problem geometry
	 */
	void make_geom_xform(wgeometry);
        /**
         * \brief makes a geometry with a transform
         * \details makes the top level group/acceleration as children of the top level object. for 
         * each primitive in the geometry, creates the geometry type, sets the intersection and 
         * bounding box programs, sets the hit programs to the geometry material, sets the program 
         * variables for the instance, creates the instances, sets cell-specific variables, makes the 
         * geometry group for the primitive, puts the geometry instance into its corresponding group, 
         * makes any necessary transforms, and attaches to the root node.
         * @param[in] problem_geom - problem geometry
         */
	void make_geom_xform_common(wgeometry);
	/**	
	 * \brief makes a primitive geometry (no transform)
	 * \details makes the top level group/acceleration as children of the top level object. for
	 * each primitive in the geometry, creates the geometry types, sets the intersection and
	 * bounding box programs, sets the hit programs to the geometry material, sets the program 
	 * variables for the instance, creates the instances, sets cell-specific variables, makes the
	 * geometry group fro the primitive, and puts the geometry instance into its group.
	 * @param[in] problem_geom - problem geometry
	 */
	void make_geom_prim(wgeometry);
	/**
	 * \brief initializes internal variables needed for OptiX stuff
	 * \details sets compute device and acceleration type; sets geometry and image types; creates 
	 * OptiX context; gets compute device information; sets up scene information; sets stack size;
	 * renders the buffers for particles, reactions, done flags, cell numbers, material numbers, 
	 * and remaps; attaches all buffers to variables; gets CUDA pointers for buffer variables; 
	 * creates programs for ray generation, exceptions, and misses; sets boundary condition for 
	 * outer cell; sets trace type; sets the outer cell and gets its dimensions; creates all 
	 * geometry instances; and validates and compiles the context.
	 * @param[in] problem_geom - problem geometry
	 * @param[in] compute_device_in - compute device to use (always zero)
	 * @param[in] accel_type_in - acceleration type
	 */ 
	void init_internal(wgeometry, unsigned, std::string);
public:
	CUdeviceptr 	positions_ptr; 
	CUdeviceptr 	      rxn_ptr; 
	CUdeviceptr 	     done_ptr;
	CUdeviceptr 	  cellnum_ptr;
	CUdeviceptr 	   matnum_ptr;
	CUdeviceptr 	    remap_ptr;
	unsigned 			stack_size_multiplier;
	unsigned 			N;
	float 				outer_cell_dims[6];
	unsigned 			outer_cell_type;
	optix_stuff(unsigned,unsigned);
	optix_stuff();
	~optix_stuff();
	void init(wgeometry, unsigned, std::string);
	void trace();
	void trace(unsigned);
	void trace(unsigned, unsigned);
	void set_trace_type(unsigned);
	void print();
	void trace_geometry(unsigned,unsigned,std::string,std::string);
	float trace_test();
	void make_color(float*,unsigned,unsigned,unsigned);
	float get_rand();
	void set_image_type(std::string);
	unsigned get_outer_cell();
};

#endif
