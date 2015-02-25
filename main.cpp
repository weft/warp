#include "warp.h"

int main(int argc, char* argv[]){

	//using namespace std;

	///////////////////
	// BEGIN PROGRAM //
	///////////////////

	//print banner
	print_banner();

	// names
	unsigned tallycell = 999;
	unsigned N = 0;
	std::string tallyname, filename;
	std::string assemblyname = "assembly";
	std::string homfuelname  = "homfuel";
	std::string godivaname   = "godiva";
	std::string pincellname  = "pincell";
	std::string pincellname2  = "pincell-square";


	// check
	if(argc<=2){
		printf("MUST ENTER A RUN TYPE : %s, %s, %s, or %s; and a number of particles to run!\n",assemblyname.c_str(),homfuelname.c_str(), godivaname.c_str(),  pincellname.c_str() );
		exit(0);
	}

	// get number of histories to do
	N = atoi(argv[2]);

	// set up geometry
	std::vector<float> mins  (3);
	std::vector<float> maxs  (3);
	std::vector<float> origin(3);
	unsigned type, material, n_topes;
	unsigned prim_id    = 0;
	wgeometry geom;

	if(assemblyname.compare(argv[1])==0){
		//assembly mats
		n_topes    = 4;
		std::vector<unsigned> 	topes      (n_topes);
		std::vector<float> 		fracs_fuel (n_topes);
		std::vector<float> 		fracs_water(n_topes);
		topes[0] = 92235;
		topes[1] = 92238;
		topes[2] = 8016;
		topes[3] = 1001;
		fracs_fuel[0] = 0.1;  
		fracs_fuel[1] = 0.9;   
		fracs_fuel[2] = 2;   
		fracs_fuel[3] = 0;
		fracs_water[0] = 0;  
		fracs_water[1] = 0;   
		fracs_water[2] = 1;   
		fracs_water[3] = 2;
		float    dens_fuel = 15;
		float 	 dens_water = 3;
		geom.add_material(1,1,n_topes,dens_fuel, topes,fracs_fuel);
		geom.add_material(2,0,n_topes,dens_water,topes,fracs_water);

		// run stuff
		tallycell = 316;   //center pin
		filename = assemblyname;
		tallyname = assemblyname;
		tallyname.append(".tally");
	
		// assembly geom
		type=1;
		material=1;
		mins[0]=-1.0;
		mins[1]=-1.0;
		mins[2]=-20.0;
		maxs[0]= 1.0;
		maxs[1]= 1.0;
		maxs[2]= 20.0;
		origin[0]=0.0;
		origin[1]=0.0;
		origin[2]=0.0;
		prim_id=geom.add_primitive(type,material,mins,maxs,origin);
		geom.make_hex_array(prim_id,15,0.0,0.0,1.164,1); 

		type=0;
		material=2;
		mins[0]=-48;
		mins[1]=-48;
		mins[2]=-48;
		maxs[0]=48;
		maxs[1]=48;
		maxs[2]=48;
		prim_id=geom.add_primitive(type,material,mins,maxs,origin);
		geom.add_transform(prim_id,999,0,0,0,0,0);
		//geom.print_all();
		//geom.update();
		//geom.print_summary();
	}
	else if(homfuelname.compare(argv[1])==0){
		// homogenized UO2 
		n_topes = 4;
		std::vector<unsigned> topes (n_topes);
		std::vector<float>    fracs (n_topes);
		topes[0]=92235,
		topes[1]=92238;
		topes[2]=8016;
		topes[3]=1001;
		fracs[0]=0.1;
		fracs[1]=0.9;  
		fracs[2]=3.0;  
		fracs[3]=2.0;
		float 	 dens = 10;
		geom.add_material(1,1,n_topes,dens,topes,fracs);
		
		// run stuff
		tallycell = 999;
		filename = homfuelname;
		tallyname = homfuelname;
		tallyname.append(".tally");
		
		//simple geom
		type=0;
		material=1;
		mins[0]=-30.0;
		mins[1]=-30.0;
		mins[2]=-30.0;
		maxs[0]= 30.0;
		maxs[1]= 30.0;
		maxs[2]= 30.0;
		origin[0]=0.0;
		origin[1]=0.0;
		origin[2]=0.0;
		prim_id=geom.add_primitive(type,material,mins,maxs,origin);
		geom.add_transform(prim_id,999,0,0,0,0,0);
	}
	else if(godivaname.compare(argv[1])==0){
		// godiva mats
		n_topes = 1;
		std::vector<unsigned> topes (n_topes);
		std::vector<float>    fracs (n_topes);

		topes[0] = 94239;
		fracs[0] = 1;      
		float    dens = 19.816;
		geom.add_material(1,1,n_topes,dens,topes,fracs);
		
		// run stuff
		tallycell = 999;
		filename = godivaname;
		tallyname = godivaname;
		tallyname.append(".tally");
	
		//godiva geom
		type=3;
		material=1;
		mins[0]= -5.1;
		mins[1]= -5.1;
		mins[2]= -5.1;
		maxs[0]=  5.1;
		maxs[1]=  5.1;
		maxs[2]=  5.1;
		origin[0]=0.0;
		origin[1]=0.0;
		origin[2]=0.0;
		prim_id=geom.add_primitive(type,material,mins,maxs,origin);
		geom.add_transform(prim_id,999,0,0,0,0,0);
		
	}
	else if(pincellname.compare(argv[1])==0){
		// pincell mats
		n_topes = 4;
		std::vector<unsigned> topes (n_topes);
		std::vector<float>    fracs_fuel  (n_topes);
		std::vector<float>    fracs_water (n_topes);
		topes[0]=92235,
		topes[1]=92238;
		topes[2]=8016;
		topes[3]=1001;
		fracs_fuel[0] = 0.1;  
		fracs_fuel[1] = 0.9;   
		fracs_fuel[2] = 2;   
		fracs_fuel[3] = 0;
		fracs_water[0] = 0;  
		fracs_water[1] = 0;   
		fracs_water[2] = 1;   
		fracs_water[3] = 2;
		float    dens_fuel = 15;
		float 	 dens_water = 3;
		geom.add_material(1,1,n_topes,dens_fuel, topes,fracs_fuel);
		geom.add_material(2,0,n_topes,dens_water,topes,fracs_water);
		
		// run stuff
		tallycell = 1;
		filename = pincellname;
		tallyname = pincellname;
		tallyname.append(".tally");
	
		//pin cell
		type=1;
		material=1;
		mins[0]=-1;
		mins[1]=-1;
		mins[2]=-20;
		maxs[0]= 1; 
		maxs[1]= 1; 
		maxs[2]= 20;
		origin[0]=0.0;
		origin[1]=0.0;
		origin[2]=0.0;
		prim_id=geom.add_primitive(type,material,mins,maxs,origin);
		geom.add_transform(prim_id,1,0,0,0,0,0);

		//water 
		type=0;
		material=2;
		mins[0]=-5.0;
		mins[1]=-5.0;
		mins[2]=-25.0;
		maxs[0]= 5.0;
		maxs[1]= 5.0;
		maxs[2]= 25.0;
		origin[0]=0.0;
		origin[1]=0.0;
		origin[2]=0.0;
		prim_id=geom.add_primitive(type,material,mins,maxs,origin);
		geom.add_transform(prim_id,999,0,0,0,0,0);
	}
	else if(pincellname2.compare(argv[1])==0){
		// pincell mats
		n_topes = 4;
		std::vector<unsigned> topes (n_topes);
		std::vector<float>    fracs_fuel  (n_topes);
		std::vector<float>    fracs_water (n_topes);
		topes[0]=92235,
		topes[1]=92238;
		topes[2]=8016;
		topes[3]=1001;
		fracs_fuel[0] = 0.1;  
		fracs_fuel[1] = 0.9;   
		fracs_fuel[2] = 2;   
		fracs_fuel[3] = 0;
		fracs_water[0] = 0;  
		fracs_water[1] = 0;   
		fracs_water[2] = 1;   
		fracs_water[3] = 2;
		float    dens_fuel = 15;
		float 	 dens_water = 3;
		geom.add_material(1,1,n_topes,dens_fuel, topes,fracs_fuel);
		geom.add_material(2,0,n_topes,dens_water,topes,fracs_water);
		
		// run stuff
		tallycell = 1;
		filename = pincellname2;
		tallyname = pincellname2;
		tallyname.append(".tally");
	
		//pin cell
		type=0;
		material=1;
		mins[0]=-1;
		mins[1]=-1;
		mins[2]=-20;
		maxs[0]= 1; 
		maxs[1]= 1; 
		maxs[2]= 20;
		origin[0]=0.0;
		origin[1]=0.0;
		origin[2]=0.0;
		prim_id=geom.add_primitive(type,material,mins,maxs,origin);
		geom.add_transform(prim_id,1,0,0,0,0,0);

		//water 
		type=0;
		material=2;
		mins[0]=-5.0;
		mins[1]=-5.0;
		mins[2]=-25.0;
		maxs[0]= 5.0;
		maxs[1]= 5.0;
		maxs[2]= 25.0;
		origin[0]=0.0;
		origin[1]=0.0;
		origin[2]=0.0;
		prim_id=geom.add_primitive(type,material,mins,maxs,origin);
		geom.add_transform(prim_id,999,0,0,0,0,0);
	}
	else{
		printf("MUST ENTER A *VALID* RUN TYPE : %s, %s, %s, or %s\n",assemblyname.c_str(),homfuelname.c_str(), godivaname.c_str(),  pincellname.c_str() );
		exit(0);
	}

	// finalize geom
	geom.set_outer_cell(999);
	geom.update();
	if(geom.check()){std::cout << "geometry failed check!\n"; return 1;}
	//geom.print_all();
	geom.print_summary();

	///////////////////////////////////
	// INIT OptiX STUFF for plotting //
	///////////////////////////////////

	// trace geom if requested
	// make new context that fits the reqested image size, trace, then destroy to free resources
	//unsigned geom_width  = 1024; 
	//unsigned geom_height = 1024;
	//unsigned N_geom = geom_width*geom_height;
	//optix_stuff geom_optix ( N_geom , 4 );
	//geom_optix.init(geom,0,"Sbvh");
	//geom_optix.trace_geometry(geom_width,geom_height,"geom.png");
	//geom_optix.~optix_stuff();


	/////////////////////////////////////////////////////////////////
	// INIT CUDA and HISTORY STUFF and LOAD/UNIONIZE CROS SECTIONS //
	/////////////////////////////////////////////////////////////////

	whistory hist ( N , geom );
	hist.set_device(0);
	hist.init();
	hist.print_xs_data();
	hist.print_materials_table();

	/////////////////////////////////////
	// converge fission source and run //
	/////////////////////////////////////

	hist.set_run_type("criticality");
	hist.set_tally_cell(tallycell);
	hist.set_run_param(40,20);  //run, skip
	hist.set_filename(filename);
	hist.run();
	hist.write_tally(0);
	//hist.write_xs_data("xsdata");

	return 0;

}

