#speedDistribution <- function(tr, grid=NULL, timestepSize=60, time.scale=timestepSize,
#		xc=NULL, yc=NULL, grid.dim=100, grid.pad=0.2) {
#	tr <- na.omit(tr) # Filter out missing measurements, these break the algorithm
#	
#	for (i in 1:length(tr)) {
#		burst <- tr[[i]]
#		burst$date <- as.double(burst$date)
#		tr[[i]] <- burst
#	}
#		
#	if (inherits(grid,"asc")) {
#		# extract the grid lines from the provided grid
#		xc <- seq(from=attr(grid, "xll"), by=attr(grid, "cellsize"), length.out=attr(grid, "dim")[1])
#		yc <- seq(from=attr(grid, "yll"), by=attr(grid, "cellsize"), length.out=attr(grid, "dim")[2])
#	} else if (!is.null(grid)) {
#		stop("grid must be an instance of 'asc' if set")
#	}
#	
#	if (is.null(xc) || is.null(yc)) {
#		g <- .defaultGrid(tr, grid.dim, grid.pad)
#		xc <- g$x
#		yc <- g$y
#	}
#	
#	speed.avg <- matrix(NA, nrow=length(xc), ncol=length(yc), dimnames=list(xc, yc))
#	for (x in xc) {
#		for (y in yc) {
#			print(paste('(', x, ',', y, ')'))
#			sum <- weight <- 0
#			for (burst in tr) {
#				ts <- range(as.double(burst$date)) + (0.5 * c(time.scale, -time.scale))
#				ts <- seq(ts[1], ts[2], by=timestepSize)
#				
#				pos <- position(tr, ts)[1,,]
#				# Compute the probability of being at (x,y) at every time step
#				pr <- apply(pos, 2, function(p) {
#					exp(((x-p['x'])^2 + (y-p['y'])^2) / (-2 * p['var'])) / (2 * pi * p['var'])
#				})
#				i <- 1 # s.t. burst$date[i] <= t < burst$date[i+1]
#				for (t in ts) {
#					while (t >= burst$date[i+1]) {
#						i <- i+1
#					}

##					if ((x-pos['x',as.character(t)])^2 + (y-pos['y',as.character(t)])^2 < 9.210314522 * pos['var',as.character(t)]) {
##						next # only look inside the 99th percentile circle
##					}
#				
#					# Determine position distributions at t - dt/2 and t + dt/2
#					t.int <- t + c(-0.5, 0.5) * time.scale
#					
#					if (t.int[1] < burst$date[i]) {
#						# the start of the time interval is in the previous bridge
#						alpha <- (burst$date[i+1] - t) / burst$dt[i]
#						beta <- (burst$date[i] - t.int[1]) / burst$dt[i-1]
#						
#						mu.s <- burst[i-1,c('x','y')]
#						mu.i <- burst[i,c('x','y')]
#						mu.e <- burst[i+1,c('x','y')]
#						var.x <- ((1-alpha)/alpha)^2 * burst$loc.var[i+1] + 
#							burst$dt[i] * (1-alpha) / alpha * burst$diff.coeff[i]
#							
#						mu <- beta * mu.s + (1-beta) * (
#								  var.x * mu.i
#								+ burst$loc.var[i] * (c(x,y) - (1-alpha)*mu.e) / alpha
#							) / (var.x + burst$loc.var[i])

#						pos.s <- c(unlist(mu),
#							var= beta^2 * burst$loc.var[i-1]
#								+ burst$dt[i-1] * beta * (1-beta) * burst$diff.coeff[i-1]
#								+ (1-beta)^2 * var.x * burst$loc.var[i] / (var.x + burst$loc.var[i])
#						)
#					} else {
#						# The start of the time interval is in the same bridge as t
#						alpha <- (t - burst$date[i]) / burst$dt[i]
#						beta <- (t.int[1] - burst$date[i]) / (t - burst$date[i])
#						
#						mu.s <- burst[i,c('x','y')]
#						mu.e <- burst[i+1,c('x','y')]
#						var.x <- (burst$dt[i] * alpha * (1-alpha) * burst$diff.coeff[i]
#								+ (alpha/(1-alpha))^2 * burst$loc.var[i+1])
#						
#						mu <- beta * c(x,y) + (1-beta) * (
#							  mu.s * var.x
#							+ (c(x,y) - alpha * mu.e) / (1-alpha) * burst$loc.var[i]
#						) / (var.x + burst$loc.var[i])
#						
#						pos.s <- unlist(c(mu,
#							var=(t-burst$date[i]) * beta * (1-beta) * burst$diff.coeff[i] +
#								(1-beta)^2 * var.x * burst$loc.var[i] / (var.x + burst$loc.var[i])
#						))
#					}
#					
#					if (t.int[2] > burst$date[i+1]) {
#						# the end of the time interval is in the next bridge
#						alpha <- (t - burst$date[i]) / burst$dt[i]
#						beta <- (t.int[2] - burst$date[i+1]) / burst$dt[i+1]
#						
#						mu.s <- burst[i,c('x','y')]
#						mu.i <- burst[i+1,c('x','y')]
#						mu.e <- burst[i+2,c('x','y')]
#						var.x <- ((1-alpha)/alpha)^2 * burst$loc.var[i] + 
#							burst$dt[i] * (1-alpha) / alpha * burst$diff.coeff[i]
#							
#						mu <- beta * mu.e + (1-beta) * (
#								  var.x * mu.i
#								+ burst$loc.var[i+1] * (c(x,y) - (1-alpha)*mu.s) / alpha
#							) / (var.x + burst$loc.var[i+1])

#						pos.s <- c(unlist(mu),
#							var= beta^2 * burst$loc.var[i+2]
#								+ burst$dt[i+1] * beta * (1-beta) * burst$diff.coeff[i+1]
#								+ (1-beta)^2 * var.x * burst$loc.var[i+1] / (var.x + burst$loc.var[i+1])
#						)

#					} else {
#						# The end of the time interval is in the same bridge as t
#						alpha <- (burst$date[i+1] - t) / burst$dt[i]
#						beta <- (burst$date[i+1] - t.int[2]) / (burst$date[i+1] - t)
#						
#						mu.s <- burst[i,c('x','y')]
#						mu.e <- burst[i+1,c('x','y')]
#						var.x <- (burst$dt[i] * alpha * (1-alpha) * burst$diff.coeff[i]
#								+ (alpha/(1-alpha))^2 * burst$loc.var[i])
#						
#						mu <- beta * c(x,y) + (1-beta) * (
#							  mu.e * var.x
#							+ (c(x,y) - alpha * mu.s) / (1-alpha) * burst$loc.var[i+1]
#						) / (var.x + burst$loc.var[i+1])
#						
#						pos.e <- unlist(c(mu,
#							var=(burst$date[i+1] - t) * beta * (1-beta) * burst$diff.coeff[i] +
#								(1-beta)^2 * var.x * burst$loc.var[i+1] / (var.x + burst$loc.var[i+1])
#						))
#					}
#					
#					# pos.s and pos.e are now independent, as the fixed position is between them in time
#					# Get the Rice distribution parameters
#					nu <- sqrt(sum((pos.e[c('x','y')] - pos.s[c('x','y')])^2)) / time.scale
#					sigma <- sqrt(pos.s['var'] + pos.e['var']) / time.scale
#					
#					avg.speed <- sigma * sqrt(pi / 2) * lmomco::LaguerreHalf(-0.5 * (nu/sigma)^2)
#					sum <- sum + pr[as.character(t)] * avg.speed
##					sum <- sum + pr[as.character(t)] * nu
#					weight <- weight + pr[as.character(t)]
#				}
#			}
#			# The average speed is the weighted average over time (weighted by probability of being at (x,y)
#			speed.avg[as.character(x),as.character(y)] <- sum / weight
#		}
#	}
#	speed.avg
#}

