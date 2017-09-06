## Clean up workspace, load required packages
	rm(list=ls())
	library(ggplot2)
	library(rasterVis)
	library(moveBB)
	

## Import data
	## The package contains a small GPS data set collected from  two vervet monkeys over two days.
	## We now load these data.
	data("vervet_monkeys", package="moveBB")
	
	## We look at some sample records. The data set is stored in a data.frame. The most important fields are:
	##  - GroupID:    An identifier for the monkey (representing a group) from which the record was obtained.
	##  - GroupDayNo: The data set is divided in bursts: subsets of consecutive measurements between which
	##                we should interpolate. This field identifies the burst to which the record belongs.
	##  - DateTime:   The (local) date and time at which the measurement was obtained.
	##  - X, Y:       The coordinates.
	##  - StdDev:     We assume that the true location is normally distributed around the measured location.
	##                This is the standard deviation of the distribution, based on the reported GPS error.
	head(monkey.data, 10)

	## Plot the measured locations, using one colour for each monkey
	plot(monkey.data$X, monkey.data$Y, pch= 21, bg= (1 + as.numeric(monkey.data$GroupID)), cex= .95)

## Calculate daily trajectories
	## We transform the data.frame holding the measurements into an object of type moveBBStack.
	## This is what the package internally works with. moveBBStack is based on the class moveStack from the package move,
	## Which is based on a SpatialPointsDataFrame holding the relocations and other data.
	## Compared to MoveStack, we added a field for the variances of the measured locations.
	## In addition, there's a field for the diffusion coefficient of Brownian motion, which is estimated automatically.
	tr <- moveBB(monkey.data$X, monkey.data$Y, monkey.data$StdDev^2, monkey.data$DateTime,
			animal=monkey.data$GroupDayNo, data=monkey.data, proj="+proj=utm +zone=36 +south ellps=WGS84")
	# tr <- monkey.tr # Already in the example data set

	## Project the positions to a more convenient coordinate reference system
	# Already projected, otherwise used: tr <- spTransform(tr, "+proj=utm +zone=36 +south ellps=WGS84")
	## The diffusion coefficient is now in wrong units, correct it by setting the optimal values.
	## Compute the most likely diffusion coefficient for each animal in the trajectory.
	## groupBy indicates the attribute by which we group the trajectories such that we get one result per animal.
	## This is done automatically on creation (without grouping), but can be run manually if desired.
	diff <- diffusionCoefficient(tr, groupBy="GroupID")
	diff
	diffusion(tr) <- diff
	
	## Printing tr prints some general summary information about the data set,
	## such as the extent of the coordinates and the ranges of other variables.
	tr
	
	## Printing tr lists the bursts with the group they belong to, the number of measured (and missing) locations,
	## The beginning and end times and some info about the trajectory (day journey length, displacement length and direction).
	summary(tr)
	
	## Collect per-relocation statistics for the trajectory
	## By default they are added to the "data" slot of the MoveBB(Stack),
	## using add=FALSE you get a data.frame with the statistics.
	tr <- statistics(tr)
	head(tr@data)
	
	## Compute the utilization distribution for each burst separately
	## By default, this function determines a default grid on which to compute the UD, which can be adjusted.
	## The user may also provide a grid on which the UD is to be calculated.
	dayUDs <- utilizationDistribution(tr, grid.dim=150)
	## Plot the daily UDs using shading and draw their 99% home range (default level for the contour function).
	plot(dayUDs)
# 	contour(dayUDs)
	
	## Get plots and contours in a single figure
	par(mfrow=c(2,2))
	lapply(unstack(dayUDs), function(ud) { plot(ud); contour(getVolumeUD(ud), add=T, levels=c(0.5, 0.99)) })
	
	
	## Compute one UD for each group
	groupUDs <- utilizationDistribution(tr, groupBy="GroupID")
	## Plot these UDs using a the same visualization, but draw the 50% home ranges as well
	plot(groupUDs)
	#contour(groupUDs, add=T, levels=c(0.5, 0.99))
	
	cp <- contourPolygons(groupUDs)
	cp$NH
	
	## Compute a lists of days and group pairs where encounters may occur.
	## For each pair of bursts that overlap in time, the expected encounter time is reported in the BBMM and linear model
	encListFull <- encounterDuration(tr, 100, timestepSize=60)
	encListFull
	
	## Keep only the interesting days:
	## Linear model detects encounters, or BBMM predicts at least 15 minutes
	encList <- encListFull[encListFull$BBMM >= 900 | encListFull$Linear > 0,]
	
	## Give the entries better labels: "<ID1>-<ID2> day <dayNo>" and sort by these labels
	#rownames(encList) <- sprintf("%s-%s day %s", encList$id1, encList$id2, 
	#		# take the part after the "_" for each burst name
	#		substr(encList$burst1,regexpr("_", encList$burst1)+1, nchar(encList$burst1)))
	encList <- encList[order(rownames(encList)),] # Sort by ascending row names

	## Plot the predicted encounter times in a bar chart
	x11()
	barplot(encList, col=hcl(c(0,120)), units="hours",
			main="Daily duration of encounters between groups",
			xlab="Group IDs and day",
			ylab="Encounter duration (hours)")
	
	
	## Convert the daily encounter durations to total duration of encounters for each pair of groups.
	## These values are stored in an array, containing one layer for each movement model.
	encMatrix <- encounterDurationById(encListFull, groupBy="GroupID")
	encMatrix
	
	## Plot the total encounter durations
	x11()
	barplot(encMatrix, col=hcl(c(0,120)), units="hours",
			main="Total duration of encounters between groups",
			xlab="Group IDs", ylab="Encounter duration (hours)")

## Compute encounters between trajectories for one specific day

	## Extract the bursts for a single day
	trd <- moveBBStack(split(tr)[c('BD_1417', 'NH_1417')])
	trd
	
	## At which times do we evaluate distance and encounter probility?
	## Take 200 regularly spaced time stamps on the intersection of both bursts' time intervals.
	tRange <- sapply(sapply(split(trd), "slot", 'timestamps'), range) # Date range for each burst
	# This is ugly, tRange is converted to type double by the outer call to sapply()
	distance.times <- as.POSIXct( # Convert back to the POSIXct type
		seq(max(tRange[1,]), min(tRange[2,]), length.out=200),
		origin="1970-01-01",tz=attr(trd@timestamps, 'tzone') # Make sure these get the time zone from the data set
	)
	
	x11()
	## Plot the probability of encounter (dist <= 100m) using both the BBMM and linear model
	## pdistance evaluates the CDF of the distance distribution
	## The result of the call to pdistance is a 4D array, indexed by two IDs, distance and time.
	## We select the IDs and distance we are interested in to get a time series of probabilities.
	encounter.plot <- qplot(distance.times,
		pdistance(100, trd, distance.times)[1,2,"100",],
		geom="line", xlab="Time", ylab="Probability of encounter")
	## Compute when encounters occur according to the linear model.
	## Since this can be solved analytically,
	## we have a function that reports the time intervals over which encounters occur
	#encounter.intervals <- encounterIntervals(100, trd)[[1,2]]
	#encounter.intervals
	## For each detected interval, add a gray rectangle to the plot
	#if (ncol(encounter.intervals) > 0) {
	#	encounter.plot <- encounter.plot + geom_rect(aes(NULL, NULL,
	#		xmin=encounter.intervals["start",], xmax=encounter.intervals["end",]),
	#		ymin=-Inf, ymax=Inf, fill=grey(0.8))
	#}
	#encounter.plot <- encounter.plot + geom_line() # re-add the line, it was overwritten by the rectangles
	encounter.plot + theme_bw()
	
	x11()
	## draw distance according to linear model and 5th and 95th percentile for BBMM
	distance.plot <- qplot(distance.times,
		distance(trd, distance.times)[1,2,], geom="line",
		xlab="Time", ylab="Distance (meters)")
	## qdistance computes quantiles for the distance distribution.
	## The result of the call to pdistance is a 4D array, indexed by two IDs, quantile level and time.
	## Select the IDs we need already, then draw a ribbon bounded by the 5% and 95% quantiles
	q.distance <- qdistance(c(0.05, 0.95), trd, distance.times)[1,2,,]
	distance.plot <- distance.plot + geom_ribbon(aes(ymin=q.distance["0.05",],
		ymax=q.distance["0.95",]), fill=grey(0.8))
	distance.plot <- distance.plot + geom_line() # re-add the line, it was overwritten by the ribbon
	distance.plot + theme_bw()

	## Compute a speed distribution on a time scale of 5 minutes
	sp.d  <- speedDistribution(tr, timestepSize=60, time.scale= 300, grid.dim=300, grid.pad=0.5)
		
	## TODO: Clip the results to the 99% home range, apparently my existing implementation is broken
	plot(sp.d)
	
	## Other working functions:
	## Compute statistics about speed at given times
	# speed    : linear model speed
	# mu.speed : mean speed in BBMM
	# dspeed   : density of speed distribution
	# pspeed   : CDF of speed distribution
	# qspeed   : quantile of speed distribution
	
	## direction, ddirection, pdirection (quantiles don't make much sense here)
	## distance,  ddistance,  pdistance,  qdistance  -- ddistance computes probability of encounter
	## position: gives mean position and variance
	## velocity: gives mean velocity and variance
	
	
## TODO: add a few features that are still missing
##   - Encounter distributions? Probably needs improvements first

tr <- split(tr)[[1]]
	
ts <- sort(unique(tr@timestamps))
dcs <- cbind(ts, seq_len(length(ts)))
dcs	
diffusion(tr) <- list(dcs)
	
sapply(tr@diffusion, function(d) {
	d(tr@timestamps)	
})
