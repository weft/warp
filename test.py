#! /usr/bin/env python
import warp

# init setup container
geom = warp.wgeometry()

# make materials
n_topes    = 4
prim_id    = 0
assemblyname = "assembly"

topes= warp.Unsigned(n_topes)
fracs_fuel = warp.Float(n_topes)
fracs_water = warp.Float(n_topes)
mins    = warp.Float(3)
maxs    = warp.Float(3)
origin  = warp.Float(3)
topes[0] = 92235
topes[1] = 92238
topes[2] = 8016
topes[3] = 1001
fracs_fuel[0] = 0.1 
fracs_fuel[1] = 0.9   
fracs_fuel[2] = 2
fracs_fuel[3] = 0
fracs_water[0] = 0  
fracs_water[1] = 0   
fracs_water[2] = 1   
fracs_water[3] = 2
dens_fuel = 15
dens_water = 3
geom.add_material(1,1,n_topes,dens_fuel, topes,fracs_fuel)
geom.add_material(2,0,n_topes,dens_water,topes,fracs_water)

# run name stuff
tallycell = 316   #center pin
filename = assemblyname
tallyname = assemblyname
tallyname = assemblyname+".tally"

# assembly geom
typ=1
material=1
mins[0]=-1.0
mins[1]=-1.0
mins[2]=-20.0
maxs[0]= 1.0
maxs[1]= 1.0
maxs[2]= 20.0
origin[0]=0.0
origin[1]=0.0
origin[2]=0.0
prim_id=geom.add_primitive(typ,material,mins,maxs,origin)
print prim_id
geom.make_hex_array(prim_id,15,0.0,0.0,1.164,1)

typ=0
material=2
mins[0]=-48
mins[1]=-48
mins[2]=-48
maxs[0]=48
maxs[1]=48
maxs[2]=48
prim_id=geom.add_primitive(typ,material,mins,maxs,origin)
print prim_id
idx=geom.add_transform(prim_id,999,0,0,0,0,0)

# finalize geom and check
geom.set_outer_cell(999)
geom.update()
geom.check()
#geom.print_all()
geom.print_summary()


# init hist and run
hist = warp.whistory(1000,geom)
hist.init()
hist.print_xs_data()
hist.print_materials_table()
hist.set_run_type("criticality")
hist.set_tally_cell(tallycell)
hist.set_run_param(40,20)
hist.set_filename(filename)
hist.run()
hist.write_tally(0)

