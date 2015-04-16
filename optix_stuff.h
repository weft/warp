#ifndef OPTIX_STUFF_H
#define OPTIX_STUFF_H
#include <optix_world.h>

/**
 * \brief OptiX stuff class
 */

class optix_stuff{
	optix::Context 	context; /**< OptiX context */
	std::string accel_type; /**< acceleration type */ 
	std::string traverse_type; /**< traverse type */
	std::string image_type; /**< image type */
	unsigned mincell; /**< minimum (usually innermost) cell */
	unsigned maxcell; /**< maximum (usually outermost) cell */
	unsigned outer_cell; /**< outermost cell */
	unsigned boundary_condition; /**< boundary condition of outermost cell */
	unsigned outer_cell_type; /**< outermost cell type */
	unsigned n_materials; /**< number of materials */
	unsigned compute_device; /**< compute device number */
	unsigned optix_device; /**< optix device number, always zero since the optix device list should only have the specified cuda device in it*/
	unsigned GEOM_FLAG; /**< geometry flag: 0 = primitive instancing, 1 = transform instancing, 2 = transform instancing with common primitives */
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
	CUdeviceptr 	positions_ptr; /**< CUDA positions pointer */ 
	CUdeviceptr 	      rxn_ptr; /**< CUDA reactions pointer */
	CUdeviceptr 	     done_ptr; /**< CUDA done flags pointer */
	CUdeviceptr 	  cellnum_ptr; /**< CUDA cell numbers pointer */
	CUdeviceptr 	   matnum_ptr; /**< CUDA material numbers pointer */
	CUdeviceptr 	    remap_ptr; /**< CUDA remaps pointer */
	unsigned 		stack_size_multiplier; /**< stack size multiplier */
	unsigned 		N;		       /**< number of histories */
	float 			outer_cell_dims[6];    /**< outermost cell dimensions */
	/**
	 * \brief constructor
	 * \details sets stack size multiplier and number of histories.
	 * @param[in] Nin - number of histories
	 * @param[in] mult - stack size multiplier
	 */ 
	optix_stuff(unsigned,unsigned);
	/**
	 * \brief "default" constructor
	 * \details empty
	 */ 
	optix_stuff();
	/**
	 * \brief destructor
	 */ 
	~optix_stuff();
	/**
	 * \brief initialization function
	 * \details sets minimum and maximum cell numbers, gets material numbers, tries to 
	 * initialize OptiX and throws an error if not.
	 * @param[in] problem_geom - problem geometry
	 * @param[in] compute_device_in - compute device to use (always 0)
	 * @param[in] accel_type_in - acceleration type
	 */ 
	void init(wgeometry, unsigned, std::string);
	/**
	 * \brief creates a trace
	 * \details launches the trace on the compute device with N histories.
	 */ 
	void trace();
	/**
	 * \brief creates a trace
	 * \details sets the trace type, then launches the trace on the compute device 
	 * with N histories.
	 * @param[in] trace_type - trace type for OptiX context
	 */ 
	void trace(unsigned);
	/**
	 * \brief creates a trace
	 * \details sets the trace type, then launches the trace on the compute device 
	 * with n_active histories.
	 * @param[in] trace_type - trace type for OptiX context
	 * @param[in] n_active - number of active histories
	 */ 
	void trace(unsigned, unsigned);
	/**
	 * \brief sets trace type in the OptiX context
	 */ 
	void set_trace_type(unsigned);
	/**
	 * \brief prints a summary of OptiX information
	 * \details prints out instancing, image type, compute device, acceleration type, 
	 * traverse type, stack size, and print buffer size.
	 */
	void print();
	/**
	 * \brief traces and plots a geometry
	 * \details gets the aspect ratio of the geometry, initializes starting points, 
	 * copies starting positions to a pointer, traces the geometry, and generates 
	 * an image. currently all commented out.
	 * @param[in] width_in,height_in - width and height of geometry object
	 * @param[in] type - type, doesn't seem to be used
	 * @param[in] filename - filename of image created
	 */ 
	void trace_geometry(unsigned,unsigned,std::string,std::string);
	/**
	 * \brief function to test tracing
	 * \details gets cell minimza and maxima, randomizes starting positions, 
	 * copies starting positions to a pointer, traces a place to generate an image,
	 * copies the data to a local buffer, creates the images, makes the distribution 
	 * random, copies the data to a pointer, executes and times the trace, and returns
	 * the time it took do to the trace.
	 * \returns time_out - time taken to do the trace
	 */
	float trace_test();
	/**
	 * \brief creates a color map
	 * @param[in] color - color map
	 * @param[in] x - used to check for a miss or normalize the color 
	 * @param[in] min,max - values used to normalize the color  
	 */
	void make_color(float*,unsigned,unsigned,unsigned);
	/**
	 * \brief returns a random float
	 */ 
	float get_rand();
	/**
	 * \brief sets image type
	 * @param[in] string_in - image type to be set
	 */
	void set_image_type(std::string);
	/**
	 * \brief returns the outermost cell
	 * \returns outer_cell - number of outermost cell
	 */  
	unsigned get_outer_cell();
	/**
	 * \brief returns the outermost cell type
	 * \returns outer_cell_type - geometrical primitive tpye of the outermost cell
	 */  
	unsigned get_outer_cell_type();
};

#endif
