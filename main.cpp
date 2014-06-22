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
	std::string homfuelname = "homfuel";
	std::string godivaname   = "godiva";
	std::string pincellname  = "pincell";

	// check
	if(argc<=2){
		printf("MUST ENTER A RUN TYPE : %s, %s, %s, or %s; and a number of particles to run!\n",assemblyname.c_str(),homfuelname.c_str(), godivaname.c_str(),  pincellname.c_str() );
		exit(0);
	}

	// get number of histories to do
	N = atoi(argv[2]);

	// set up geometry
	wgeometry geom;

	if(assemblyname.compare(argv[1])==0){
		//assembly mats
		unsigned n_topes    = 4;
		unsigned topes      [n_topes];
		float    fracs_fuel [n_topes];
		float 	 fracs_water[n_topes];
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
		geom.add_primitive();
		geom.primitives[0].type=1;
		geom.primitives[0].material=1;
		geom.primitives[0].min[0]=-1.0;
		geom.primitives[0].min[1]=-1.0;
		geom.primitives[0].min[2]=-20.0;
		geom.primitives[0].max[0]= 1.0;
		geom.primitives[0].max[1]= 1.0;
		geom.primitives[0].max[2]= 20.0;
		geom.primitives[0].make_hex_array(15,0.0,0.0,1.164,0.0,1); 
		geom.add_primitive();
		geom.primitives[1].type=0;
		geom.primitives[1].material=2;
		geom.primitives[1].min[0]=-48;
		geom.primitives[1].min[1]=-48;
		geom.primitives[1].min[2]=-48;
		geom.primitives[1].max[0]=48;
		geom.primitives[1].max[1]=48;
		geom.primitives[1].max[2]=48;
		geom.primitives[1].add_transform();
		geom.primitives[1].transforms[0].cellnum = 999;
		geom.primitives[1].transforms[0].dx      = 0;
		geom.primitives[1].transforms[0].dy      = 0;
		geom.primitives[1].transforms[0].dz      = 0;
		geom.primitives[1].transforms[0].theta   = 0;
		geom.primitives[1].transforms[0].phi     = 0;
	}
	else if(homfuelname.compare(argv[1])==0){
		// homogenized UO2 
		unsigned topes[4]={92235,92238,8016,1001};
		float    fracs[4]={.1 , 0.9,   3,   2};
		float 	 dens = 10;
		geom.add_material(1,1,4,dens,topes,fracs);
		
		// run stuff
		tallycell = 999;
		filename = homfuelname;
		tallyname = homfuelname;
		tallyname.append(".tally");
		
		//simple geom
		geom.add_primitive();
		geom.primitives[0].type=0;
		geom.primitives[0].material=1;
		geom.primitives[0].min[0]=-30.0;
		geom.primitives[0].min[1]=-30.0;
		geom.primitives[0].min[2]=-30.0;
		geom.primitives[0].max[0]= 30.0;
		geom.primitives[0].max[1]= 30.0;
		geom.primitives[0].max[2]= 30.0;
		geom.primitives[0].add_transform();
		geom.primitives[0].transforms[0].cellnum = 999;
		geom.primitives[0].transforms[0].dx      = 0;
		geom.primitives[0].transforms[0].dy      = 0;
		geom.primitives[0].transforms[0].dz      = 0;
		geom.primitives[0].transforms[0].theta   = 0;
		geom.primitives[0].transforms[0].phi     = 0;
	}
	else if(godivaname.compare(argv[1])==0){
		// godiva mats
		unsigned n_topes    = 1;
		unsigned topes      [n_topes];
		float    fracs_fuel [n_topes];
		float 	 fracs_water[n_topes];
		topes[0] = 94239;
		fracs_fuel[0] = 1;      
		float    dens_fuel = 19.816;
		geom.add_material(1,1,n_topes,dens_fuel, topes,fracs_fuel);
		
		// run stuff
		tallycell = 999;
		filename = godivaname;
		tallyname = godivaname;
		tallyname.append(".tally");
	
		//godiva geom
		geom.add_primitive();
		geom.primitives[0].type=3;
		geom.primitives[0].material=1;
		geom.primitives[0].min[0]= -5.1;
		geom.primitives[0].min[1]= -5.1;
		geom.primitives[0].min[2]= -5.1;
		geom.primitives[0].max[0]=  5.1;
		geom.primitives[0].max[1]=  5.1;
		geom.primitives[0].max[2]=  5.1;
		geom.primitives[0].add_transform();
		geom.primitives[0].transforms[0].cellnum = 999;
		geom.primitives[0].transforms[0].dx      = 0;
		geom.primitives[0].transforms[0].dy      = 0;
		geom.primitives[0].transforms[0].dz      = 0;
		geom.primitives[0].transforms[0].theta   = 0;
		geom.primitives[0].transforms[0].phi     = 0;
	}
	else if(pincellname.compare(argv[1])==0){
		// pincell mats
		unsigned n_topes    = 4;
		unsigned topes      [n_topes];
		float    fracs_fuel [n_topes];
		float 	 fracs_water[n_topes];
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
		tallycell = 1;
		filename = pincellname;
		tallyname = pincellname;
		tallyname.append(".tally");
	
		//pin cell
		geom.add_primitive(); //pin
		geom.primitives[0].type=1;
		geom.primitives[0].material=1;
		geom.primitives[0].min[0]=-1;
		geom.primitives[0].min[1]=-1;
		geom.primitives[0].min[2]=-20;
		geom.primitives[0].max[0]= 1; 
		geom.primitives[0].max[1]= 1; 
		geom.primitives[0].max[2]= 20;
		geom.primitives[0].add_transform();
		geom.primitives[0].transforms[0].cellnum = 1;
		geom.primitives[0].transforms[0].dx      = 0;
		geom.primitives[0].transforms[0].dy      = 0;
		geom.primitives[0].transforms[0].dz      = 0;
		geom.primitives[0].transforms[0].theta   = 0;
		geom.primitives[0].transforms[0].phi     = 0;
		geom.add_primitive();  //clad 
		geom.primitives[1].type=0;
		geom.primitives[1].material=2;
		geom.primitives[1].min[0]=-5.0;
		geom.primitives[1].min[1]=-5.0;
		geom.primitives[1].min[2]=-25.0;
		geom.primitives[1].max[0]= 5.0;
		geom.primitives[1].max[1]= 5.0;
		geom.primitives[1].max[2]= 25.0;
		geom.primitives[1].add_transform();
		geom.primitives[1].transforms[0].cellnum = 999;
		geom.primitives[1].transforms[0].dx      = 0;
		geom.primitives[1].transforms[0].dy      = 0;
		geom.primitives[1].transforms[0].dz      = 0;
		geom.primitives[1].transforms[0].theta   = 0;
		geom.primitives[1].transforms[0].phi     = 0;
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

