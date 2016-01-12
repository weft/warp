#ifndef WGEOMETRY_H
#define WGEOMETRY_H
#include "wprimitive.h"

/**
 * \class wgeometry wgeometry.h
 * \brief wgeometry class
 */

class wgeometry {
public:
	unsigned 	n_box; /**<  number of boxes */
	unsigned 	n_cyl; /**<  number of cylinders */
	unsigned 	n_hex; /**<  number of hexagons */
	unsigned 	n_sph; /**<  number of spheres */
	unsigned 	n_primitives; /**<  number of primitives */
	unsigned 	n_transforms; /**<  number of transforms */
	unsigned 	outer_cell; /**<  outermost cell (usually used for tallying) */
	unsigned 	n_materials; /**<  number of materials */
	unsigned 	n_isotopes; /**<  number of isotopes */
	unsigned 	fissile_flag; /**<  indicates whether or not a material is fissile */
	unsigned	boundary_condition; /**<  flag for the cell's boundary condition */
	unsigned * 	material_num_list; /**<  list of material numbers */
	unsigned * 	cell_num_list; /**<  list of cell numbers */
	std::string datapath; /**< path to xsdir and data */
	/**
	 *  wgeometry constructor
	 */ 
	 wgeometry();
	/**
	 *  wgeometry destructor
	 */ 
	~wgeometry();
	/**
	 * returns the smallest cell number, typically the innermost cell
	 */
	unsigned get_minimum_cell();
	/**
	 * returns the largest cell number, typically the outermost cell
	 */ 
	unsigned get_maximum_cell();
	/**
	 * returns the smallest material number
	 */
	unsigned get_minimum_material();
	/**
	 * returns the largest material number
	 */ 
	unsigned get_maximum_material();
	/**
	 * returns the number of primitves in the geometry object
	 */ 
	unsigned get_primitive_count();
	/**
	 * returns the number of transforms in the geometry object
	 */ 
	unsigned get_transform_count();
	/**
	 * adds a primitive to the geometry object
	 */ 
	unsigned add_primitive();
	unsigned add_primitive(int , unsigned  ,std::vector<float>, std::vector<float>, std::vector<float> );
	/**
	 * updates the numbers of all shapes, compiles the list of all isotopes, creates
	 * an isotope table.
	 */ 
	void update();
	/**
	 * prints a summary of the geometry object: numbers of the different kinds of
	 * shapes in the geometry, total numbers of primitives and transforms, outer
	 * cell, numbers of materials and isotopes, isotope list, properties (density,
	 * fissile flag, isotopes) of each material.
	 */
	void print_summary();
	/**
	 * prints all of the transforms of all the primitives, then prints a geometry
	 * summary.
	 */ 
	void print_all();
	/**
	 * goes through all the cells of all of the primitives and checks that the
	 * outer cell is set
	 *
	 * @param[in] ocell - the outermost cell
	 */ 
	void set_outer_cell(unsigned, unsigned);
	/**
	 * returns the outermost cell and its boundary conition
	 *
	 * \returns outer_cell
	 */ 
	unsigned get_outer_cell();
	/**
	 * returns the outermost cell type
	 *
	 * \returns outer_cell_type
	 */ 
	unsigned get_boundary_condition();
	/**
	 * returns the boundary condition flag of outermost cell type
	 *
	 * \returns boundary_condition
	 */ 
	unsigned get_outer_cell_type();
	/**
	 * adds a material and its properties to the geometry, allocates space
	 * for all of the material information.
	 *
	 * @param[in] matnum - material number
	 * @param[in] is_fissile - fissile flag
	 * @param[in] num_topes - number of isotopes in material
	 * @param[in] density - density of material
	 * @param[in] isotopes - list of isotopes
	 * @param[in] fractions - fractions of the constituent isotopes
	 */ 
	void add_material(unsigned , unsigned, unsigned , float, std::vector<std::string> , std::vector<float> );
	/**
	 * checks that all cells have unique IDs, checks that there are materials for
	 * each number specified in the geometry, checks to make sure that the outer 
	 * cell exists, checks to see if there are any fissile isotopes.
	 */ 
	int check();
	/**
	 * returns the dimensions of the outermost cell.
	 * 
	 * @param[in] input_array
	 * \returns primitives[k].type
	 */ 
	unsigned get_outer_cell_dims(float*);
	/**
	 * returns the number of materials.
	 *
	 * \returns n_materials
	 */
	unsigned get_material_count();
	/**
	 * makes a table of all of the materials.
	 *
	 * allocates and copies the isotope and material number lists to their respsective
	 * arrays, allocates and copies the isotope fractions to the concentration
	 * matrix, converts the fractions into number densities, normalizes the
	 * fractions, gets the average number density, prints each isotope's material,
	 * isotope, and density.
	 */
	void make_material_table();
	/**
	 * creates material and isotope arrays, creates concentration matrix. copies
	 * memory for all of those arrays.
	 *
	 * @param[in] n_mat_in - number of input materials
	 * @param[in] n_tope_in - number of input isotopes 
	 * @param[in] material_list_in - list of input materials
	 * @param[in] isotope_list_in - list of input isotopes
	 * @param[in] conc_mat_in - input concentration matrix
	 */
	//void get_material_table(unsigned*,unsigned*,unsigned**,unsigned**,float**);
	void get_material_table(unsigned*,unsigned*,float**);
	/**
	 * prints out all materials, including each material's constituent isotopes 
	 * and their number densities.
	 */ 
	void print_materials_table();
	/**
	 * checks whether or not the geometry contains a fissile material.
	 *
	 * \returns fissile_flag
	 */
	void set_datapath(std::string);
	/**
	 * sets the data path
	 *
	 * \returns void
	 */	
	unsigned check_fissile();
	unsigned add_transform(unsigned);
	unsigned add_transform(unsigned, unsigned, float, float, float, float, float);
	unsigned add_transform(unsigned, unsigned, unsigned, float, float, float, float, float);
	void make_hex_array(unsigned, int, float, float, float, unsigned);
	void delete_primitive(unsigned);
	void delete_transform(unsigned,unsigned);
	std::vector<primitive>   	primitives; /**< primitives vector */
	std::vector<material_def>	materials;  /**< materials vector */
	std::vector<std::string>	isotopes;   /**< isotopes vector */
	std::string 			    isotope_list; /**< isotope list */
	unsigned *	isotope_list_array; /**< isotope list array */
	unsigned *	material_list_array; /**< material list array */
	float * 	concentrations_matrix; /**< concentrations matrix */
	float * 	awr_list; /**< atomic weight ratio (AWR) list */
};

#endif
