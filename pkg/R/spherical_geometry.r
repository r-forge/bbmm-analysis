## Some functions related to cartesion and spherical geometry
## functions with .ll are lat/long, with .ct are cartesian

ll2xyz <- function(ll) {
  ## Convert lat/long coordinates (in radians) to 3d Cartesian (on unit sphere)
  c( cos(ll[1])*c(cos(ll[2]),sin(ll[2])) , sin(ll[1]) )
}

xyz2ll <- function(xyz) {
  ## Convert 3d Cartesian coordinates to lat/long (in radians)
  ## Also computes angles for coordinates not on the unit sphere
  c( atan2(xyz[3], sqrt(sum(xyz[1:2]^2))), atan2(xyz[2],xyz[1]) )
}

dist.ct <- function(p,q) {
  sqrt(sum((p-q)^2))
}

dist.ll <- function(p,q) {
  ## Great-circle distance between p,q in lat/long coordinates (in radians)
  ## Normalized to unit sphere
  
  # Calculates the geodesic distance using the Haversine formula
  d <- q-p
  a <- sin(d[1]/2)^2 + cos(p[1]) * cos(q[1]) * sin(d[2]/2)^2
  2 * asin(min(1,sqrt(a)))
}

midpoint.ct <- function(p,q) {
  ## Midpoint of p and q in Cartesian coordinates (arbitrary d)
  (p+q)*0.5
}

midpoint.ll <- function(p,q) {
  ## Midpoint of shortest great-circle arc between p,q (lat/lon in radians)
  ## Undefined for antipodal points
  
  ## Convert to Cartesian, compute midpoint, convert back
  xyz2ll( midpoint.ct( ll2xyz(p), ll2xyz(q) ) )
}

circumcircle.ct <- function(p,q,r) {
  ## Compute circumcircle of triangle embedded in higher dimensional space
  ## (p,q,r in Cartesian coordinates, d arbitrary)
  ## Based on https://en.wikipedia.org/wiki/Circumscribed_circle#Higher_dimensions
  l  <- function(x) { sqrt(sum(x^2)) } # shortcut for vector length
  l2 <- function(x) { sum(x^2) } # squared vector length
  
  ## Translate p to origin
  q <- q-p
  r <- r-p
  
  # Compute translated circumcenter
  qxr2 <- l2(q)*l2(r) - sum(q*r)^2 # ||q cross r||^2
  d <- (l2(q)*r - l2(r)*q)
  m <- (sum(d*r)*q - sum(d*q)*r) / (2*qxr2)
  # Compute radius
  r <- l(q) * l(r) * l(q-r) / (2*sqrt(qxr2))
  
  # translate back circumcenter and return
  c(m+p, r)
}

circumcircle.ll <- function(p, q, r) {
  ## Midpoint of smallest circumcircle of p,q,r (lat/long in radians)
  ## Undefined for p,q,r and origin coplanar
  
  ## Convert to Cartesian, compute circumcircle, convert back
  m <- xyz2ll(circumcircle.ct( ll2xyz(p), ll2xyz(q), ll2xyz(r) )[1:3])
  c(m, dist.ll(m,p)) # Recompute radius along great circle
}

