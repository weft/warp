#include "warp.h"

int main(int argc, char* argv[]){

	//using namespace std;

	///////////////////
	// BEGIN PROGRAM //
	///////////////////

	// names
	unsigned tallycell = 999, outer_cell=999;
	unsigned N = 0;
	std::string tallyname, filename, runtype;
	std::string homfuelname		= "homfuel";
	std::string assemblyname  	= "assembly-lw";
	std::string flibename   	= "flibe";
	std::string fusionname    	= "fusion";
	std::string guidename   	= "guide";
	std::string jezebelname  	= "jezebel";
	std::string pincellname  	= "pincell";
	std::string sodiumpinname  	= "sodiumpin";
	std::string testname 	 	= "test";
	std::string card0name 	 	= "debug";
	std::string card1name 	 	= "k20";
	std::string card2name 	 	= "k80";
	std::string card3name 	 	= "titan";


	// check
	if(argc<=2){
		printf("MUST ENTER A RUN TYPE : ");
		printf("%s, ",assemblyname.c_str());
		printf("%s, ",flibename.c_str());
		printf("%s, ",fusionname.c_str());
		printf("%s, ",guidename.c_str());
		printf("%s, ",jezebelname.c_str());
		printf("%s, ",pincellname.c_str());
		printf("%s, ",sodiumpinname.c_str());
		printf("%s, ",testname.c_str());
		printf("and a number of neutrons per cycle!\n");
		//printf("and a card name (%s %s %s %s)!\n",card0name.c_str(),card1name.c_str(),card2name.c_str(),card3name.c_str());
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
	unsigned bc 		= 1;
	unsigned dev		= 0;
	wgeometry geom;

	// set datapath
	//geom.set_datapath("/usr/local/LANL/MCNP6_DATA/xsdir_mcnp6.1");
	//geom.set_datapath("/usr/local/SERPENT/xsdata/endfb7/sss_endfb7u.xsdir");
	geom.set_datapath("/scratch/bergmann_r/SERPENT/xsdata/endfb7/sss_endfb7u.xsdir");

	if(assemblyname.compare(argv[1])==0){
		//assembly mats
		n_topes = 9;
		std::vector<std::string> topes (n_topes);
		std::vector<float>    fracs_fuel  (n_topes);
		std::vector<float>    fracs_water (n_topes);
		std::vector<float>    fracs_clad  (n_topes);
		topes[0] = "92235.03c";
		topes[1] = "92238.03c";
		topes[2] =  "8016.03c" ;
		topes[3] =  "1001.03c" ;
		topes[4] = "40090.03c";
 		topes[5] = "40091.03c";
 		topes[6] = "40092.03c";
 		topes[7] = "40094.03c";
 		topes[8] = "40096.03c";
		fracs_fuel[0] = 0.1;  
		fracs_fuel[1] = 0.9;   
		fracs_fuel[2] = 2;   
		fracs_fuel[3] = 0;
		fracs_fuel[4] = 0;
		fracs_fuel[5] = 0;
		fracs_fuel[6] = 0;  
		fracs_fuel[7] = 0;
		fracs_fuel[8] = 0;
	  	fracs_water[0] = 0; 
		fracs_water[1] = 0; 
		fracs_water[2] = 1;   
		fracs_water[3] = 2;
		fracs_water[4] = 0;
		fracs_water[5] = 0;
		fracs_water[6] = 0;  
		fracs_water[7] = 0;
		fracs_water[8] = 0;
	    fracs_clad[0] = 0; 
		fracs_clad[1] = 0; 
		fracs_clad[2] = 0;   
		fracs_clad[3] = 0;
		fracs_clad[4] = 0.5145;
		fracs_clad[5] = 0.1122;
		fracs_clad[6] = 0.1715;  
		fracs_clad[7] = 0.1738;
		fracs_clad[8] = 0.0280;
	   
	
		float    dens_fuel = 10.97;
		float 	 dens_water = 1.00;
		float 	 dens_clad = 6.52;
		geom.add_material(1,1,n_topes,dens_fuel, topes,fracs_fuel);
		geom.add_material(2,0,n_topes,dens_water,topes,fracs_water);
		geom.add_material(3,0,n_topes,dens_clad,topes,fracs_clad);

		// run stuff
		tallycell = 316;   //center pin
		outer_cell = 3000;
		filename  = assemblyname;
		tallyname = assemblyname;
		tallyname.append(".tally");
		bc = 1;
		runtype = "criticality";
	
		// assembly geom, fuel
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
		geom.make_hex_array(prim_id,15,0.0,0.0,1.3,1); 

		// assembly geom, clad
		type=1;
		material=3;
		mins[0]=-1.2;
		mins[1]=-1.2;
		mins[2]=-20.2;
		maxs[0]= 1.2;
		maxs[1]= 1.2;
		maxs[2]= 20.2;
		origin[0]=0.0;
		origin[1]=0.0;
		origin[2]=0.0;
		prim_id=geom.add_primitive(type,material,mins,maxs,origin);
		geom.make_hex_array(prim_id,15,0.0,0.0,1.3*2.0/2.4,1000); 

		// water 
		type=0;
		material=2;
		mins[0]=-48;
		mins[1]=-48;
		mins[2]=-48;
		maxs[0]=48;
		maxs[1]=48;
		maxs[2]=48;
		prim_id=geom.add_primitive(type,material,mins,maxs,origin);
		geom.add_transform(prim_id,outer_cell,0,0,0,0,0);
		//geom.print_all();
		//geom.update();
		//geom.print_summary();
	}
	else if(flibename.compare(argv[1])==0){
		//flibe pebble mats
		n_topes    = 8;
		std::vector<std::string> 	topes      (n_topes);
		std::vector<float> 			fracs_fuel (n_topes);
		std::vector<float> 			fracs_flibe(n_topes);
		topes[0] = "92235.12c";
		topes[1] = "92238.12c";
		topes[2] =  "8016.12c";
		topes[3] =  "6000.12c";
		topes[4] =  "3006.09c";
		topes[5] =  "3007.09c";
		topes[6] =  "4009.09c";
		topes[7] =  "9019.09c";

		fracs_fuel[0] = 0.10 ;  
		fracs_fuel[1] = 0.90;   
		fracs_fuel[2] = 2;   
		fracs_fuel[3] = 2;
		fracs_fuel[4] = 0;  
		fracs_fuel[5] = 0;   
		fracs_fuel[6] = 0;   
		fracs_fuel[7] = 0;

		fracs_flibe[0] = 0;  
		fracs_flibe[1] = 0;   
		fracs_flibe[2] = 0;   
		fracs_flibe[3] = 0;
		fracs_flibe[4] = 0.15;  
		fracs_flibe[5] = 1.85;   
		fracs_flibe[6] = 1;   
		fracs_flibe[7] = 4;

		float    dens_fuel  = 8.75;
		float 	 dens_flibe = 1.94;

		geom.add_material(1,1,n_topes,dens_fuel, topes,fracs_fuel);
		geom.add_material(2,0,n_topes,dens_flibe,topes,fracs_flibe);

		// run stuff
		tallycell = 1;   //center pin
		filename  = flibename;
		tallyname = flibename;
		tallyname.append(".tally");
		bc = 2;
		runtype = "criticality";
	
		// flibe geom
		type=3;
		material=1;
		mins[0]=-5.0;
		mins[1]=-5.0;
		mins[2]=-5.0;
		maxs[0]= 5.0;
		maxs[1]= 5.0;
		maxs[2]= 5.0;
		origin[0]=0.0;
		origin[1]=0.0;
		origin[2]=0.0;
		prim_id=geom.add_primitive(type,material,mins,maxs,origin);
		geom.add_transform(prim_id,1,0,0,0,0,0);

		type=2;
		material=2;
		mins[0]=-5.1;
		mins[1]=-5.1;
		mins[2]=-5.1;
		maxs[0]= 5.1;
		maxs[1]= 5.1;
		maxs[2]= 5.1;
		prim_id=geom.add_primitive(type,material,mins,maxs,origin);
		geom.add_transform(prim_id,999,0,0,0,0,0);
	}
	else if(fusionname.compare(argv[1])==0){
		// fusion mats
		n_topes = 18;
		std::vector<std::string> topes  (n_topes);
		std::vector<float>  void_fracs  (n_topes);
		std::vector<float>   sic_fracs  (n_topes);
		std::vector<float>    li_fracs  (n_topes);
		std::vector<float>    ss_fracs  (n_topes);
		topes[0]  =  "6000.03c";
		topes[1]  = "14028.03c";
		topes[2]  = "14029.03c";
		topes[3]  = "14030.03c";
		topes[4]  =  "3006.03c";
		topes[5]  =  "3007.03c";
		topes[6]  = "26054.03c";
		topes[7]  = "26056.03c";
		topes[8]  = "26057.03c";
		topes[9]  = "26058.03c";
		topes[10] = "24050.03c";
		topes[11] = "24052.03c";
		topes[12] = "24053.03c";
		topes[13] = "24054.03c";
		topes[14] = "28058.03c";
		topes[15] = "28060.03c";
		topes[16] = "28062.03c";
		topes[17] = "28064.03c";

		void_fracs[0]  = 0.0;
		void_fracs[1]  = 0.0;
		void_fracs[2]  = 0.0;
		void_fracs[3]  = 0.0;
		void_fracs[4]  = 0.0;
		void_fracs[5]  = 0.0;
		void_fracs[6]  = 0.0;
		void_fracs[7]  = 0.0;
		void_fracs[8]  = 0.0;
		void_fracs[9]  = 0.0;
		void_fracs[10] = 0.0;
		void_fracs[11] = 0.0;
		void_fracs[12] = 0.0;
		void_fracs[13] = 0.0;
		void_fracs[14] = 0.0;
		void_fracs[15] = 0.0;
		void_fracs[16] = 0.0;
		void_fracs[17] = 0.0;

		sic_fracs[0]  = 1.0;
		sic_fracs[1]  = 0.9223;
		sic_fracs[2]  = 0.0467;
		sic_fracs[3]  = 0.0310;
		sic_fracs[4]  = 0.0;
		sic_fracs[5]  = 0.0;
		sic_fracs[6]  = 0.0;
		sic_fracs[7]  = 0.0;
		sic_fracs[8]  = 0.0;
		sic_fracs[9]  = 0.0;
		sic_fracs[10] = 0.0;
		sic_fracs[11] = 0.0;
		sic_fracs[12] = 0.0;
		sic_fracs[13] = 0.0;
		sic_fracs[14] = 0.0;
		sic_fracs[15] = 0.0;
		sic_fracs[16] = 0.0;
		sic_fracs[17] = 0.0;

		li_fracs[0]  = 0.0;
		li_fracs[1]  = 0.0;
		li_fracs[2]  = 0.0;
		li_fracs[3]  = 0.0;
		li_fracs[4]  = 0.075;
		li_fracs[5]  = 0.925;
		li_fracs[6]  = 0.0;
		li_fracs[7]  = 0.0;
		li_fracs[8]  = 0.0;
		li_fracs[9]  = 0.0;
		li_fracs[10] = 0.0;
		li_fracs[11] = 0.0;
		li_fracs[12] = 0.0;
		li_fracs[13] = 0.0;
		li_fracs[14] = 0.0;
		li_fracs[15] = 0.0;
		li_fracs[16] = 0.0;
		li_fracs[17] = 0.0;

		ss_fracs[0]  = 0.0;
		ss_fracs[1]  = 0.0;
		ss_fracs[2]  = 0.0;
		ss_fracs[3]  = 0.0;
		ss_fracs[4]  = 0.0;
		ss_fracs[5]  = 0.0;
		ss_fracs[6]  = 0.0435;
		ss_fracs[7]  = 0.6879;
		ss_fracs[8]  = 0.0165;
		ss_fracs[9]  = 0.0021;
		ss_fracs[10] = 0.0065;
		ss_fracs[11] = 0.1257;
		ss_fracs[12] = 0.0143;
		ss_fracs[13] = 0.0035;
		ss_fracs[14] = 0.0681;
		ss_fracs[15] = 0.0262;
		ss_fracs[16] = 0.0036;
		ss_fracs[17] = 0.0009;

		float   void_dens= 0.000;
		float   sic_dens = 3.210;
		float    li_dens = 0.534;
		float    ss_dens = 7.990;
		geom.add_material(0,0,n_topes,void_dens,topes,void_fracs);
		geom.add_material(1,0,n_topes, sic_dens,topes, sic_fracs);
		geom.add_material(2,0,n_topes,  li_dens,topes,  li_fracs);
		geom.add_material(3,0,n_topes,  ss_dens,topes,  ss_fracs);
		
		// run stuff
		tallycell = 999;
		filename  = fusionname;
		tallyname = fusionname;
		tallyname.append(".tally");
		bc = 1;
		runtype = "fixed";
	
		//fusion geom
		type=3;
		material=0;
		mins[0]= -150.0;
		mins[1]= -150.0;
		mins[2]= -150.0;
		maxs[0]=  150.0;
		maxs[1]=  150.0;
		maxs[2]=  150.0;
		origin[0]=0.0;
		origin[1]=0.0;
		origin[2]=0.0;
		prim_id=geom.add_primitive(type,material,mins,maxs,origin);
		geom.add_transform(prim_id,0,0,0,0,0,0);

		type=3;
		material=1;
		mins[0]= -151.0;
		mins[1]= -151.0;
		mins[2]= -151.0;
		maxs[0]=  151.0;
		maxs[1]=  151.0;
		maxs[2]=  151.0;
		origin[0]=0.0;
		origin[1]=0.0;
		origin[2]=0.0;
		prim_id=geom.add_primitive(type,material,mins,maxs,origin);
		geom.add_transform(prim_id,1,0,0,0,0,0);

		type=3;
		material=2;
		mins[0]= -156.0;
		mins[1]= -156.0;
		mins[2]= -156.0;
		maxs[0]=  156.0;
		maxs[1]=  156.0;
		maxs[2]=  156.0;
		origin[0]=0.0;
		origin[1]=0.0;
		origin[2]=0.0;
		prim_id=geom.add_primitive(type,material,mins,maxs,origin);
		geom.add_transform(prim_id,2,0,0,0,0,0);

		type=3;
		material=3;
		mins[0]= -160.0;
		mins[1]= -160.0;
		mins[2]= -160.0;
		maxs[0]=  160.0;
		maxs[1]=  160.0;
		maxs[2]=  160.0;
		origin[0]=0.0;
		origin[1]=0.0;
		origin[2]=0.0;
		prim_id=geom.add_primitive(type,material,mins,maxs,origin);
		geom.add_transform(prim_id,3,0,0,0,0,0);

	}
	else if(guidename.compare(argv[1])==0){
		
	}
	else if(jezebelname.compare(argv[1])==0){
		// jezebel mats
		n_topes = 6;
		std::vector<std::string> topes (n_topes);
		std::vector<float>    fracs (n_topes);
		topes[0] = "94239.03c";
		topes[1] = "94240.03c";
		topes[2] = "94241.03c";
		topes[3] = "94242.03c";
		topes[4] = "31069.03c";
		topes[5] = "31071.03c";
		fracs[0] = 0.029934;
		fracs[1] = 0.0078754;
		fracs[2] = 0.0012146;
		fracs[3] = 0.00015672;
		fracs[4] = 0.00082482942;
		fracs[5] = 0.00054737058;
		float    dens = 15.730;
		geom.add_material(1,1,n_topes,dens,topes,fracs);
		
		// run stuff
		tallycell = 999;
		filename  = jezebelname;
		tallyname = jezebelname;
		tallyname.append(".tally");
		bc = 1;
		runtype = "criticality";
	
		//jezebel geom
		type=3;
		material=1;
		mins[0]= -6.6595;
		mins[1]= -6.6595;
		mins[2]= -6.6595;
		maxs[0]=  6.6595;
		maxs[1]=  6.6595;
		maxs[2]=  6.6595;
		origin[0]=0.0;
		origin[1]=0.0;
		origin[2]=0.0;
		prim_id=geom.add_primitive(type,material,mins,maxs,origin);
		geom.add_transform(prim_id,999,0,0,0,0,0);
		
	}
	else if(pincellname.compare(argv[1])==0){
		// pincell mats
		n_topes = 9;
		std::vector<std::string> topes (n_topes);
		std::vector<float>    fracs_fuel  (n_topes);
		std::vector<float>    fracs_water (n_topes);
		std::vector<float>    fracs_clad  (n_topes);
		topes[0] = "92235.03c";
		topes[1] = "92238.03c";
		topes[2] =  "8016.03c" ;
		topes[3] =  "1002.03c" ;
		topes[4] = "40090.03c";
 		topes[5] = "40091.03c";
 		topes[6] = "40092.03c";
 		topes[7] = "40094.03c";
 		topes[8] = "40096.03c";
		fracs_fuel[0] = 0.1;  
		fracs_fuel[1] = 0.9;   
		fracs_fuel[2] = 2;   
		fracs_fuel[3] = 0;
		fracs_fuel[4] = 0;
		fracs_fuel[5] = 0;
		fracs_fuel[6] = 0;  
		fracs_fuel[7] = 0;
		fracs_fuel[8] = 0;
	  	fracs_water[0] = 0; 
		fracs_water[1] = 0; 
		fracs_water[2] = 1;   
		fracs_water[3] = 2;
		fracs_water[4] = 0;
		fracs_water[5] = 0;
		fracs_water[6] = 0;  
		fracs_water[7] = 0;
		fracs_water[8] = 0;
		fracs_clad[0] = 0; 
		fracs_clad[1] = 0; 
		fracs_clad[2] = 0;   
		fracs_clad[3] = 0;
		fracs_clad[4] = 0.5145;
		fracs_clad[5] = 0.1122;
		fracs_clad[6] = 0.1715;  
		fracs_clad[7] = 0.1738;
		fracs_clad[8] = 0.0280;
	   
	
		float    dens_fuel = 10.97;
		float 	 dens_water = 1.11;
		float 	 dens_clad = 6.52;
		geom.add_material(1,1,n_topes,dens_fuel, topes,fracs_fuel);
		geom.add_material(2,0,n_topes,dens_water,topes,fracs_water);
		geom.add_material(3,0,n_topes,dens_clad,topes,fracs_clad);
		
		// run stuff
		tallycell = 1;
		filename  = pincellname;
		tallyname = pincellname;
		tallyname.append(".tally");
		bc = 1;
		runtype = "criticality";
	
		//pin cell
		type=1;
		material=1;
		mins[0]=-2;
		mins[1]=-2;
		mins[2]=-20;
		maxs[0]= 2; 
		maxs[1]= 2; 
		maxs[2]= 20;
		origin[0]=0.0;
		origin[1]=0.0;
		origin[2]=0.0;
		prim_id=geom.add_primitive(type,material,mins,maxs,origin);
		geom.add_transform(prim_id,1,0,0,0,0,0);

		//clad
		type=1;
		material=3;
		mins[0]=-2.2;
		mins[1]=-2.2;
		mins[2]=-20.2;
		maxs[0]= 2.2; 
		maxs[1]= 2.2; 
		maxs[2]= 20.2;
		origin[0]=0.0;
		origin[1]=0.0;
		origin[2]=0.0;
		prim_id=geom.add_primitive(type,material,mins,maxs,origin);
		geom.add_transform(prim_id,200,0,0,0,0,0);

		//water 
		type=0;
		material=2;
		mins[0]=-50.0;
		mins[1]=-50.0;
		mins[2]=-25.0;
		maxs[0]= 50.0;
		maxs[1]= 50.0;
		maxs[2]= 25.0;
		origin[0]=0.0;
		origin[1]=0.0;
		origin[2]=0.0;
		prim_id=geom.add_primitive(type,material,mins,maxs,origin);
		geom.add_transform(prim_id,999,0,0,0,0,0);
	}
	else if(homfuelname.compare(argv[1])==0){
		// homfuel mats
		n_topes = 9;
		std::vector<std::string> topes (n_topes);
		std::vector<float>    fracs_fuel  (n_topes);
		topes[0] = "92235.03c";
		topes[1] = "92238.03c";
		topes[2] =  "8016.03c" ;
		topes[3] =  "1002.03c" ;
		topes[4] = "40090.03c";
 		topes[5] = "40091.03c";
 		topes[6] = "40092.03c";
 		topes[7] = "40094.03c";
 		topes[8] = "40096.03c";
		fracs_fuel[0] = 0.1;  
		fracs_fuel[1] = 0.9;   
		fracs_fuel[2] = 3;   
		fracs_fuel[3] = 2;
		fracs_fuel[4] = 0.5145;
		fracs_fuel[5] = 0.1122;
		fracs_fuel[6] = 0.1715;  
		fracs_fuel[7] = 0.1738;
		fracs_fuel[8] = 0.0280;
	   
		float    dens_fuel = 5.50;
		geom.add_material(1,1,n_topes,dens_fuel, topes,fracs_fuel);
		
		// run stuff
		tallycell = 999;
		filename  = homfuelname;
		tallyname = homfuelname;
		tallyname.append(".tally");
		bc = 1;
		runtype = "criticality";
	
		//water 
		type=0;
		material=1;
		mins[0]=-50.0;
		mins[1]=-50.0;
		mins[2]=-25.0;
		maxs[0]= 50.0;
		maxs[1]= 50.0;
		maxs[2]= 25.0;
		origin[0]=0.0;
		origin[1]=0.0;
		origin[2]=0.0;
		prim_id=geom.add_primitive(type,material,mins,maxs,origin);
		geom.add_transform(prim_id,999,0,0,0,0,0);
	}
	else if(testname.compare(argv[1])==0){
		// pincell mats
		n_topes = 4;
		std::vector<std::string> topes (n_topes);
		std::vector<float>    fracs_fuel  (n_topes);
		std::vector<float>    fracs_water  (n_topes);
		topes[0] = "92235.03c";
		topes[1] = "40090.03c";
		topes[2] = "8016.03c";
		topes[3] = "1002.03c";
		fracs_fuel[0] = 0.1;  
		fracs_fuel[1] = 0.9;
		fracs_fuel[2] = 0.0;  
		fracs_fuel[3] = 0.0;
		fracs_water[0] = 0.0;  
		fracs_water[1] = 0.0;
		fracs_water[2] = 1.0;  
		fracs_water[3] = 2.0;
	
		float    dens_fuel = 10.97;
		float 	 dens_water = 1;

		geom.add_material(1,1,n_topes,dens_fuel, topes,fracs_fuel);
		geom.add_material(2,0,n_topes,dens_water, topes,fracs_water);
		
		// run stuff
		tallycell = 1;
		filename  = testname;
		tallyname = testname;
		tallyname.append(".tally");
		bc = 1;
		runtype = "criticality";

		//pin cell
		type=1;
		material=1;
		mins[0]=-1.0;
		mins[1]=-1.0;
		mins[2]=-20;
		maxs[0]= 1.0; 
		maxs[1]= 1.0; 
		maxs[2]= 20;
		origin[0]=0.0;
		origin[1]=0.0;
		origin[2]=0.0;
		prim_id=geom.add_primitive(type,material,mins,maxs,origin);
		geom.add_transform(prim_id,1,0,0,0,0,0);

		//clad
		type=1;
		material=2;
		mins[0]=-1.2;
		mins[1]=-1.2;
		mins[2]=-20.2;
		maxs[0]= 1.2; 
		maxs[1]= 1.2; 
		maxs[2]= 20.2;
		origin[0]=0.0;
		origin[1]=0.0;
		origin[2]=0.0;
		prim_id=geom.add_primitive(type,material,mins,maxs,origin);
		geom.add_transform(prim_id,2,0,0,0,0,0);

		//water 
		type=2;
		material=2;
		mins[0]=-1.8;
		mins[1]=-1.8;
		mins[2]=-25.0;
		maxs[0]= 1.8;
		maxs[1]= 1.8;
		maxs[2]= 25.0;
		origin[0]=0.0;
		origin[1]=0.0;
		origin[2]=0.0;
		prim_id=geom.add_primitive(type,material,mins,maxs,origin);
		geom.add_transform(prim_id,999,0,0,0,0,0);

	}
	else if(sodiumpinname.compare(argv[1])==0){
		// sodium pincell mats
		n_topes = 15;
		std::vector<std::string> topes    (n_topes);
		std::vector<float>    fracs_fuel  (n_topes);
		std::vector<float>    fracs_water (n_topes);
		std::vector<float>    fracs_clad  (n_topes);
		topes[0]= "92235.09c";
		topes[1]= "92238.09c";
		topes[2]= "11023.06c" ;
		topes[3]= "26054.09c" ;
		topes[4]= "26056.09c";
 		topes[5]= "26057.09c";
 		topes[6]= "26058.09c";
 		topes[7]= "24050.09c";
 		topes[8]= "24052.09c";
 		topes[9]= "24053.09c";
 		topes[10]="24054.09c";
 		topes[11]="28058.09c";
 		topes[12]="28060.09c";
 		topes[13]="28062.09c";
 		topes[14]="28064.09c";

		fracs_fuel[0]  = 0.1;  
		fracs_fuel[1]  = 0.9;   
		fracs_fuel[2]  = 0;   
		fracs_fuel[3]  = 0;
		fracs_fuel[4]  = 0;
		fracs_fuel[5]  = 0;
		fracs_fuel[6]  = 0;  
		fracs_fuel[7]  = 0;
		fracs_fuel[8]  = 0;
		fracs_fuel[9]  = 0;
		fracs_fuel[10] = 0;
		fracs_fuel[11] = 0;
		fracs_fuel[12] = 0;
		fracs_fuel[13] = 0;
		fracs_fuel[14] = 0;

		fracs_water[0]  = 0;  
		fracs_water[1]  = 0;   
		fracs_water[2]  = 1;   
		fracs_water[3]  = 0;
		fracs_water[4]  = 0;
		fracs_water[5]  = 0;
		fracs_water[6]  = 0;  
		fracs_water[7]  = 0;
		fracs_water[8]  = 0;
		fracs_water[9]  = 0;
		fracs_water[10] = 0;
		fracs_water[11] = 0;
		fracs_water[12] = 0;
		fracs_water[13] = 0;
		fracs_water[14] = 0;

		fracs_clad[0]  = 0;  
		fracs_clad[1]  = 0;   
		fracs_clad[2]  = 0;   
		fracs_clad[3]  = 0.0435;
		fracs_clad[4]  = 0.6879;
		fracs_clad[5]  = 0.0165;
		fracs_clad[6]  = 0.0021;  
		fracs_clad[7]  = 0.0065;
		fracs_clad[8]  = 0.1257;
		fracs_clad[9]  = 0.0143;
		fracs_clad[10] = 0.0035;
		fracs_clad[11] = 0.0681;
		fracs_clad[12] = 0.0262;
		fracs_clad[13] = 0.0036;
		fracs_clad[14] = 0.0009;
	   
	
		float    dens_fuel  = 19.100;
		float 	 dens_water =  0.927;
		float 	 dens_clad  =  7.990;
		geom.add_material(1,1,n_topes,dens_fuel, topes,fracs_fuel);
		geom.add_material(2,0,n_topes,dens_water,topes,fracs_water);
		geom.add_material(3,0,n_topes,dens_clad,topes,fracs_clad);
		
		// run stuff
		tallycell = 1;
		filename  = sodiumpinname;
		tallyname = sodiumpinname;
		tallyname.append(".tally");
		bc = 2;
		runtype = "criticality";
	
		//pin cell
		type=1;
		material=1;
		mins[0]=-1.0;
		mins[1]=-1.0;
		mins[2]=-20;
		maxs[0]= 1.0; 
		maxs[1]= 1.0; 
		maxs[2]= 20;
		origin[0]=0.0;
		origin[1]=0.0;
		origin[2]=0.0;
		prim_id=geom.add_primitive(type,material,mins,maxs,origin);
		geom.add_transform(prim_id,1,0,0,0,0,0);

		//clad
		type=1;
		material=3;
		mins[0]=-1.2;
		mins[1]=-1.2;
		mins[2]=-20.2;
		maxs[0]= 1.2; 
		maxs[1]= 1.2; 
		maxs[2]= 20.2;
		origin[0]=0.0;
		origin[1]=0.0;
		origin[2]=0.0;
		prim_id=geom.add_primitive(type,material,mins,maxs,origin);
		geom.add_transform(prim_id,2,0,0,0,0,0);

		//water 
		type=2;
		material=2;
		mins[0]=-1.8;
		mins[1]=-1.8;
		mins[2]=-25.0;
		maxs[0]= 1.8;
		maxs[1]= 1.8;
		maxs[2]= 25.0;
		origin[0]=0.0;
		origin[1]=0.0;
		origin[2]=0.0;
		prim_id=geom.add_primitive(type,material,mins,maxs,origin);
		geom.add_transform(prim_id,999,0,0,0,0,0);
		
	}
	else{
		printf("MUST ENTER A *VALID* RUN TYPE : ");
		printf("%s, ",homfuelname.c_str());
		printf("%s, ",assemblyname.c_str());
		printf("%s, ",flibename.c_str());
		printf("%s, ",fusionname.c_str());
		printf("%s, ",guidename.c_str());
		printf("%s, ",jezebelname.c_str());
		printf("%s, ",pincellname.c_str());
		printf("%s, ",sodiumpinname.c_str());
		printf("and a number of particles to run!\n");
		exit(0);
	}

	// finalize geom
	geom.set_outer_cell(outer_cell,bc);  // cell, BC  1=black, 2=specular
	geom.add_tally(tallycell);
	geom.update();
	if(geom.check()){std::cout << "geometry failed check!\n"; return 1;}
	//geom.print_all();
	geom.print_summary();

	/////////////////////////////////////////////////////////////////
	// INIT CUDA and HISTORY STUFF and LOAD/UNIONIZE CROS SECTIONS //
	/////////////////////////////////////////////////////////////////

	whistory hist ( N , geom );
	hist.set_print_level(2);
	hist.set_dump_level(1);
	hist.set_device(dev);
	hist.init();
	hist.print_xs_data();
	hist.print_materials_table();

	/////////////////////////////////////
	// converge fission source and run //
	/////////////////////////////////////

	hist.set_run_type(runtype);
	hist.set_run_param(40,20);  //run, skip
	hist.set_filename(filename);
	hist.plot_geom("cell");  // **MUST** be called after init.
	hist.run();
	hist.write_tally();

	return 0;

}

