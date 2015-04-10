".onLoad" <- function(libname, pkgname) {

	# Try to load the OpenCL package and select a device to use
	.OpenCL_init() # Not yet, the OpenCL package needs a new version
	
	#load('/home/stef/School/TUe_CSE/2IM91_Master_project/SVN/src/scripts_vervet_monkeys/visibility_DEM_reduced_4.Rda')
	#load('/home/stef/School/TUe_CSE/2IM91_Master_project/SVN/src/scripts_vervet_monkeys/tr.Rda')
	#print(encounter.visible(vis, tr, as.POSIXct(c('2011-01-19 17:00:00','2011-01-19 18:00:00'))))
}

".onUnload" <- function(libpath) {
	library.dynam.unload("movementAnalysis", libpath);
}
