".onLoad" <- function(libname, pkgname) {
	# Try to load the OpenCL package and select a device to use
	.OpenCL_init() # Not yet, the OpenCL package needs a new version
	
	# Set the default number of processes to use for parallel computations
	if (is.null(options("mc.cores")$mc.cores)) {
		options(mc.cores=detectCores(logical=TRUE))
	}
}

".onUnload" <- function(libpath) {
	library.dynam.unload("moveBB", libpath);
}
