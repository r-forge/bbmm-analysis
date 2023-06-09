".MoveBB_OpenCL_dev" <- NA
".OpenCL_kernels" <- NA

# Selects an OpenCL device if possible, based on heuristics
".OpenCL_init" <- function() {
	if (try(packageVersion('OpenCL'), TRUE) == '0.1.3-1' && require(OpenCL)) {
		devices <- sapply(tryCatch(oclPlatforms(), error=function(e) {
			print(e$message)
			NULL
		}), oclDevices)
		devType <- lapply(devices, function (d) { oclInfo(d, full=TRUE)$TYPE })

		# Policy for selecting devices: There is an order of preference
		# based on the reported device type. Use the first device of the most
		# preferred type that has any devices available.
		dev <- NULL
		for (t in c("GPU", "ACCELERATOR", "CPU", "DEFAULT")) {
			t <- paste("CL_DEVICE_TYPE_", t, sep="")
			devs <- devices[vapply(devType, function (dt) { any(dt == t) }, logical(1))]
			if (length(devs) > 0) {
				d <- devs[[1]]
				dev <- d
				break # no need to check less preferred dev types
			}
		}
		.MoveBB_OpenCL_dev <<- dev # remember the selected device
		cat("Selected OpenCL device ", ifelse(is.null(dev), "None", oclInfo(dev)$name), '\n', sep="")
	}
}

## Let the user override the device selected above
#"OpenCL_setDevice" <- function(dev) {
##	if (class(dev) != "clDeviceID") {
##		stop("dev must be an instance of clDeviceID, obtained from oclDevices()")
##	}
#	unlockBinding(".MoveBB_OpenCL_dev", loadNamespace(getPackageName()))
#	unlockBinding(".OpenCL_kernels", loadNamespace(getPackageName()))
#	.MoveBB_OpenCL_dev <<- dev
#	.OpenCL_kernels <<- NA # kernels will need to be rebuilt for this device
#	lockBinding(".MoveBB_OpenCL_dev", loadNamespace(getPackageName()))
#	lockBinding(".OpenCL_kernels", loadNamespace(getPackageName()))
#}

# get an OpenCL kernel with the provided name, if it exists
".OpenCL_getKernel" <- function(name) {
	if (any(is.na(.OpenCL_kernels))) {
		unlockBinding(".OpenCL_kernels", loadNamespace(getPackageName()))
	
		if (class(.MoveBB_OpenCL_dev) == "clDeviceID") {
			oldWd <- getwd()
		
			# Load all files with a '.cl' extension (case insensitive).
			# These files are all compiled in one compilation unit.
			# Unfortunately, one cannot use #include with this scheme;
			# The files are loaded in the lexicographic order of their full path and name,
			# so make sure the order of declarations and definitions is correct.
			oclDir <- system.file("OpenCL", package=getPackageName())
			setwd(oclDir)
			files <- dir(oclDir, recursive=TRUE, full.names=TRUE)
			files <- files[grep("\\.cl$", files, ignore.case=TRUE)]
		
			source <- sapply(files, function (f) {
				paste(readLines(f), collapse="\n")
			}, USE.NAMES=FALSE)
			source <- paste(c("#pragma OPENCL EXTENSION cl_khr_fp64 : enable", source), collapse="\n\n")

			.OpenCL_kernels <<- oclSimpleKernel(.MoveBB_OpenCL_dev, NA, source, "double") # Set precision here
			
			if (!is.null(oldWd)) setwd(oldWd)
		} else {
			# OpenCL not available, so we have no kernels
			.OpenCL_kernels <<- list()
		}
		
		lockBinding(".OpenCL_kernels", loadNamespace(getPackageName()))
		print(names(.OpenCL_kernels))
	}
	return(tryCatch(.OpenCL_kernels[[name]], error=function(e) { NULL }))
}

".OpenCL" <- function(kernelName, size, ...) {
	# OpenCL is disabled for now
	kernel <- .OpenCL_getKernel(kernelName)
	if (!is.null(kernel)) {
		return(oclRun(kernel, size, ...))
	} else {
		cResult <- .C(kernelName, double(size), size, ..., PACKAGE=getPackageName())
		return(cResult[[1]])
	}
}

#OpenCL_test <- function(nu, sd, x) {
#	print("test")
#	r <- .OpenCL("riceTest", length(x), as.double(nu), as.double(sd), as.double(x))
#	print(r)
#}

