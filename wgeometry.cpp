#include <vector> 
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <cmath>
#include <assert.h>
#include <time.h>
#include <string.h>
#include "datadef.h"
#include "primitive.h"
#include "wgeometry.h"

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
//							Geometry stuff
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

wgeometry::wgeometry(){
	n_box        = 0;
	n_cyl        = 0;
	n_hex        = 0;
	n_primitives = 0;
	n_transforms = 0;
	outer_cell   = 0;
	n_materials  = 0;
	n_isotopes   = 0;
	fissile_flag = 0;
}
wgeometry::~wgeometry(){
	//material destructor
	//for(int k=0;k<n_materials;k++){
	//	delete materials[k].fractions;
	//	delete materials[k].isotopes;
	//}
	//delete cell_num_list;
	//delete material_num_list;
}
unsigned wgeometry::add_primitive(){

	primitive this_primitive;
	primitives.push_back(this_primitive);
	n_primitives++;



}
unsigned wgeometry::add_primitive(int ptype, unsigned cellmat , std::vector<float> mins, std::vector<float> maxs, std::vector<float> origin){

	primitive this_primitive(ptype, cellmat, mins, maxs, origin);
	primitives.push_back(this_primitive);
	n_primitives++;

	return (n_primitives-1);

}
void wgeometry::update(){
	n_box        = 0;
	n_cyl        = 0;
	n_hex        = 0;
	n_sph 		 = 0;
	n_transforms = 0;
	for(int k=0;k<n_primitives;k++){
		if (primitives[k].n_transforms==0){
			std::cout << "No transforms for primitive id = " << primitives[k].primitive_id << ", it will not be included in problem geometry" << "\n";
		}
		if(primitives[k].type==0){
				n_box+=primitives[k].n_transforms;
		}
		else if(primitives[k].type==1){
				n_cyl+=primitives[k].n_transforms;
		}
		else if(primitives[k].type==2){
				n_hex+=primitives[k].n_transforms;
		}
		else if(primitives[k].type==3){
				n_sph+=primitives[k].n_transforms;
		}
		n_transforms+=primitives[k].n_transforms;
	}

	// allocate arrays for lists
	cell_num_list = new unsigned [n_transforms];
	material_num_list = new unsigned [n_transforms]; // allocate enough for every cell to have its own material

	// compile list of all isotopes
	unsigned this_mat  = 0;
	unsigned n_topes   = 0;
	unsigned this_tope = 0;
	std::vector<unsigned>  all_isotopes;
	for(this_mat=0 ; this_mat<n_materials ; this_mat++){
		n_topes = materials[this_mat].num_isotopes;
		for(int k=0;k<n_topes;k++){
			all_isotopes.push_back(materials[this_mat].isotopes[k]);
		}
	}
	// go through list, get rid of extra copies
	n_isotopes = 0;
	unsigned notfound=0;
	//std::cout << "all_isotopes.size() = " << all_isotopes.size() << "\n";
	for(int k=0;k<all_isotopes.size();k++){
		notfound=1;
		for(int j=0;j<isotopes.size();j++){
			if(isotopes[j]==all_isotopes[k])
				notfound=0; 
		}
		if(notfound){
			isotopes.push_back(all_isotopes[k]);  //add if not found already
		}
	}
	n_isotopes = isotopes.size();

	//make string from isotope table
	char numstr[16];
	for(int k =0;k<n_isotopes;k++){
		sprintf(numstr,"%u",isotopes[k]);
		isotope_list += numstr;
		if(k<n_isotopes-1){
			isotope_list += ",";
		}
	}

}
void wgeometry::print_summary(){
	std::cout << "\e[1;32m" << "--- GEOMETRY SUMMARY ---" << "\e[m \n";
	std::cout << "rectangular prisms = " << n_box << "\n";
	std::cout << "cylinders          = " << n_cyl << "\n";
	std::cout << "hexagons           = " << n_hex << "\n";
	std::cout << "spheres            = " << n_sph << "\n";
	std::cout << "total primitives   = " << n_primitives << "\n";
	std::cout << "total transforms   = " << n_transforms << "\n";
	std::cout << "outer cell         = " << outer_cell << "\n";
	std::cout << "\e[1;32m" << "--- INPUT MATERIAL SUMMARY ---" << "\e[m \n";
	std::cout << "materials          = " << n_materials << "\n";
	std::cout << "isotopes           = " << n_isotopes << "\n";
	std::cout << "isotope list:    " << isotope_list << "\n";
	std::cout << "  --------------   " << n_materials << "\n";
	for(int k=0;k<n_materials;k++){
		std::cout << "material #       = " << materials[k].matnum << "\n";
		std::cout << "density (g/cc)   = " << materials[k].density << "\n";
		std::cout << "is fissile       = " << materials[k].is_fissile << "\n";
		std::cout << "isotopes         = " << materials[k].num_isotopes << "\n";
		for(int j=0;j<materials[k].num_isotopes;j++){
			std::cout << "  number "<< j << ":  isotope " << materials[k].isotopes[j] << " frac = " << materials[k].fractions[j] << "\n";
		}
		if(k!=n_materials-1){
			std::cout << "     ***\n";
		}
	}
}
void wgeometry::print_all(){
	for(int k=0;k<n_primitives;k++){
		primitives[k].print_transform();
	}
	print_summary();
}
unsigned wgeometry::get_primitive_count(){
	return(n_primitives);
}
unsigned wgeometry::get_transform_count(){
	return(n_transforms);
}
void wgeometry::set_outer_cell(unsigned ocell){
	outer_cell = ocell;
	unsigned this_cell;
	for(int j=0;j<n_primitives;j++){
		for(int k=0;k<primitives[j].n_transforms;k++){
			this_cell = primitives[j].transforms[k].cellnum;
			if (this_cell==outer_cell){break;}
		}
		if (this_cell==outer_cell){break;}
	}
	if (this_cell!=outer_cell) {
		std::cout << "Cell " << ocell << " not found, outer cell not set!!!" << "\n";
	}
}
unsigned wgeometry::get_outer_cell(){
	return outer_cell;
}
unsigned wgeometry::get_minimum_cell(){
	unsigned mincell=-1;
	for(int j=0;j<n_primitives;j++){
		for(int k=0;k<primitives[j].n_transforms;k++){
			if (primitives[j].transforms[k].cellnum<mincell){mincell=primitives[j].transforms[k].cellnum;}
		}
	}
	return mincell;
}
unsigned wgeometry::get_maximum_cell(){
	unsigned maxcell=0;
	for(int j=0;j<n_primitives;j++){
		for(int k=0;k<primitives[j].n_transforms;k++){
			if (primitives[j].transforms[k].cellnum>maxcell){maxcell=primitives[j].transforms[k].cellnum;}
		}
	}
	return maxcell;
}
void wgeometry::add_material(unsigned matnum, unsigned is_fissile, unsigned num_topes, float density, std::vector<unsigned> isotopes, std::vector<float> fractions){
	
	// get current material index
	unsigned dex = materials.size(); 

	material_def this_material_def;

	this_material_def.fractions = new float    [num_topes];
	this_material_def.isotopes  = new unsigned [num_topes];

	this_material_def.num_isotopes  = num_topes;
	this_material_def.matnum        = matnum;
	this_material_def.id 			= dex;
	this_material_def.density       = density;
	this_material_def.is_fissile    = is_fissile;
	for (int i=0; i<num_topes;i++){
		this_material_def.isotopes[i] = isotopes[i];
		this_material_def.fractions[i] = fractions[i];
	}

	// memcpy(this_material_def.fractions,  &fractions[0],   num_topes*sizeof(float));
	// memcpy(this_material_def.isotopes,   &isotopes[0],    num_topes*sizeof(unsigned));

	materials.push_back(this_material_def);

	n_materials++;
}
int wgeometry::check(){

	std::cout << "\e[1;32m" << "Checking cell numbers and materials..." << "\e[m \n";

	unsigned cellnum,matnum;
	unsigned cell_list_index = 0;
	unsigned mat_list_index  = 0;
	unsigned z,notfound,found_cell;
	// check that all cells have their own ID
	for (int k=0;k<n_primitives;k++){
		for (int j=0;j<primitives[k].n_transforms;j++){	
			cellnum = primitives[k].transforms[j].cellnum;
			matnum  = primitives[k].transforms[j].cellmat;
			// scan the cell list 
			for (z = 0 ; z<cell_list_index; z++){
				if (cell_num_list[z]==cellnum){
					std::cout << "cell number " << cellnum << " has duplicate entries!\n";
					return 1;
				}
			}
			cell_num_list[z]=cellnum; //append this cell number
			cell_list_index++;

			// scan the material list
			notfound=1;
			for (z = 0 ; z<mat_list_index ; z++){
				if (material_num_list[z]==matnum){
					notfound=0;
					break;   //break on this index if found
				}
			}
			if(notfound){
				material_num_list[mat_list_index]=matnum;  // append this material and increment index counter
				mat_list_index++;
			}
		}
	}

	// check that there are materials for each number specified in the geom
	for(int k=0;k<mat_list_index;k++){
		notfound=1;
		for(int j=0;j<n_materials;j++){
			if(material_num_list[k]==materials[j].matnum){
				notfound=0;
				break;
			}
		}
		if(notfound){
			std::cout << "material " << material_num_list[k] << " not defined!\n";
			return 1;
		}
	}

	// check to make sure the outer cell exists
	notfound = 1;
	for (int k=0;k<n_primitives;k++){
		for (int j=0;j<primitives[k].n_transforms;j++){	
			if(primitives[k].transforms[j].cellnum==outer_cell & notfound){
				notfound=0;
			}
		}
	}
	if(notfound){
		std::cout << "Cell " << outer_cell << " not found, cannot set it as the outer cell!\n";
		return 1;
	}

	//see if there are any fissile isotopes
	for(int k=0;k<n_materials;k++){
		fissile_flag += materials[k].is_fissile;
	}

	std::cout << "They check out.\n";
	return 0;

}
unsigned wgeometry::get_outer_cell_dims(float * input_array){

	float this_min[3];
	float this_max[3];

	for (int k=0;k<n_primitives;k++){
		for (int j=0;j<primitives[k].n_transforms;j++){	
			if(primitives[k].transforms[j].cellnum==outer_cell){
				// apply transform to base primitive, just translation now, maybe add rotation later?  no this is a maximum extent projection onto the axes, should always be a box.
				memcpy(this_min , primitives[k].min , 3*sizeof(float));
				memcpy(this_max , primitives[k].max , 3*sizeof(float));
				this_min[0] += primitives[k].transforms[j].dx;
				this_min[1] += primitives[k].transforms[j].dy;
				this_min[2] += primitives[k].transforms[j].dz;
				this_max[0] += primitives[k].transforms[j].dx;
				this_max[1] += primitives[k].transforms[j].dy;
				this_max[2] += primitives[k].transforms[j].dz;
				// copy and return type
				memcpy(&input_array[0] , this_min , 3*sizeof(float));
				memcpy(&input_array[3] , this_max , 3*sizeof(float));
				return primitives[k].type;
			}
		}
	}

}
unsigned wgeometry::get_material_count(){
	return n_materials;
}
unsigned wgeometry::check_fissile(){

	return fissile_flag;

}
void wgeometry::make_material_table(){

	// allocate and copy the insotope list to the array
	isotope_list_array = new unsigned [n_isotopes];
	memcpy(isotope_list_array,isotopes.data(),n_isotopes*sizeof(unsigned));

	// allocate and copy the material number list to the array
	material_list_array = new unsigned [n_materials];
	for(int k=0;k<n_materials;k++){
		material_list_array[k]=materials[k].matnum;
	}

	// allocate and copy the fractions to the matrix
	unsigned notfound=1;
	int z=0;
	concentrations_matrix = new float [n_materials*n_isotopes];
	for(int j=0;j<n_materials;j++){     // isotope in a column
		for(int k=0;k<n_isotopes;k++){  // material in a row
			
			notfound=1;
			//scan the material object to see if the isotope is there
			for(z=0;z<materials[j].num_isotopes;z++){
				if(materials[j].isotopes[z] == isotope_list_array[k]){
					notfound=0;
					break;
				}
			}

			// use the internal index to copy to matrix
			if(notfound){
				concentrations_matrix[j*n_isotopes + k] = 0.0;
			}
			else{
				concentrations_matrix[j*n_isotopes + k] = materials[j].fractions[z];
			}
		}
	}

	// now convert fractions into number densities
	float frac   = 0.0;
	float m_avg  = 0.0;
	float N_avg  = 0.0;
	float awr    = 0.0;
	float dens   = 0.0;
	float u_to_g = 1.66053892e-24; // grams
	float m_n    = 1.008664916;    // u
	float barns  = 1e24;

	for(int j=0;j<n_materials;j++){

		m_avg = 0.0;
		frac  = 0.0;

		//normalize the fractions for this material and calculate average mass
		for(int k=0;k<n_isotopes;k++){
			frac += concentrations_matrix[j*n_isotopes+k];
		}
		for(int k=0;k<n_isotopes;k++){
			concentrations_matrix[j*n_isotopes+k] = concentrations_matrix[j*n_isotopes+k]/frac;
			m_avg += concentrations_matrix[j*n_isotopes+k] * awr_list[k] * m_n;
			//std::cout << "awr["<<k<<"] = "<<awr_list[k]<<"\n";
		}

		//get density
		dens = materials[j].density;

		// average num density
		N_avg = dens/(m_avg * u_to_g * barns);

		//  multiply normalized fractions by average number density to get topes number density
		for(int k=0;k<n_isotopes;k++){
			concentrations_matrix[j*n_isotopes+k] = concentrations_matrix[j*n_isotopes+k] * N_avg;
			printf("material = %d, isotope %d, dens = %6.5E\n",j,k,concentrations_matrix[j*n_isotopes+k]);
		}
	}
}
void wgeometry::get_material_table(unsigned* n_mat_in, unsigned * n_tope_in, unsigned** material_list_in, unsigned** isotope_list_in, float** conc_mat_in){

	*n_mat_in  = n_materials;
	*n_tope_in = n_isotopes;

	*material_list_in 	= new unsigned [n_materials];
	*isotope_list_in 	= new unsigned [n_isotopes];
	*conc_mat_in 		= new float    [n_materials*n_isotopes];

	memcpy(*material_list_in,  material_list_array,    n_materials*sizeof(unsigned)         );
	memcpy(*isotope_list_in,   isotope_list_array,     n_isotopes *sizeof(unsigned)         );
	memcpy(*conc_mat_in,       concentrations_matrix,  n_materials*n_isotopes*sizeof(float) );
}
void wgeometry::print_materials_table(){

	std::cout << "\e[1;32m" << "--- MATERIALS SUMMARY ---" << "\e[m \n";

	for(int j=0;j<n_materials;j++){

		assert(j==materials[j].id);
		std::cout <<  "material index " << j << " = material " << material_list_array[j] << "\n";
		std::cout <<  " (isotope index, ZZZAAA) \n";
		std::cout <<  " (number density #/bn-cm) \n";
		
		for(int k=0;k<n_isotopes;k++){

			if (k==n_isotopes-1){
				std::cout << "( "<< k << " , "<< isotope_list_array[k] << " ) \n";
			}
			else{
				std::cout << "  ( "<< k << " , "<< isotope_list_array[k] << " )     ";
			}
		}

		for(int k=0;k<n_isotopes;k++){

			if (k==n_isotopes-1){
				std::cout << "( " <<concentrations_matrix[j*n_isotopes+k] << " )\n";
			}
			else{
				std::cout << "  ( " <<concentrations_matrix[j*n_isotopes+k] << " )     ";
			}
		}

		
	}

}
void wgeometry::add_transform(unsigned index ){

	primitives[index].add_transform();

}

void wgeometry::add_transform(unsigned index, unsigned cellnum , float dx , float dy , float dz , float theta , float phi ){

	primitives[index].add_transform(cellnum, dx, dy, dz, theta, phi);

}

void wgeometry::add_transform(unsigned index, unsigned cellnum ,unsigned cellmat, float dx , float dy , float dz , float theta , float phi ){

	primitives[index].add_transform(cellnum, cellmat, dx, dy, dz, theta, phi);

}

void wgeometry::make_hex_array(unsigned index, int n, float x, float y, float phi, unsigned starting_index){

	primitives[index].make_hex_array(n, x, y, phi, starting_index);

}
