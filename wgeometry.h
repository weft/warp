#ifndef WGEOMETRY_H
#define WGEOMETRY_H

/**
 *  wgeometry class
 */

class wgeometry {
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
	unsigned * 	material_num_list; /**<  list of material numbers */
	unsigned * 	cell_num_list; /**<  list of cell numbers */
public:
	/**
	 *  wgeometry constructor
	 *  \todo make a google test for construction
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
	void add_primitive();
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
	void set_outer_cell(unsigned);
	/**
	 * returns the outermost cell
	 *
	 * \returns outer_cell
	 */ 
	unsigned get_outer_cell();
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
	void add_material(unsigned , unsigned, unsigned , float, unsigned * , float * );
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
	 */ 
	unsigned get_outer_cell_dims(float*);
	unsigned get_material_count();
	void make_material_table();
	void get_material_table(unsigned*,unsigned*,unsigned**,unsigned**,float**);
	void print_materials_table();
	unsigned check_fissile();
	std::vector<primitive>   	primitives;
	std::vector<material_def>	materials;
	std::vector<unsigned>		isotopes;
	std::string 				isotope_list;
	unsigned *	isotope_list_array;
	unsigned *	material_list_array;
	float * 	concentrations_matrix;
	float * 	awr_list;
};

#endif
