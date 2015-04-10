".onLoad" <- function(libname, pkgname) {
	# Try to load the OpenCL package and select a device to use
	#.OpenCL_init() # Not yet, the OpenCL package needs a new version
}

".onUnload" <- function(libpath) {
	library.dynam.unload("movementAnalysis", libpath);
}
