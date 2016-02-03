# inspired by http://www.petercollingridge.co.uk/science-simulation/gel-electrophoresis-mathematics

mu2mu0 <- function(Nbp, V_cm, b) {
  # assuming 1% agarose gel
  # Nbp is length in base pairs of dsDNA fragment
  # V_cm is electric force in V/cm
  # b is pore size in nm - 1% agarose estimates range from 200-500 nm
  
  # Nk is the number of Kuhn lengths in a sequence and the estimate
  # is about 300 bp to 1 Kuhn length
  Nk = Nbp / 300
  
  # convert to V/m
  E = V_cm / 100
  
  # l is Kuhn length of dsDNA - supposedly around 100 nm
  l = 100
  
  # epsilon_k is empirically estimated (E in V/m)
  epsilon_k = 1.4e-4 * E
  
  # alpha is experimentally estimated
  alpha = 3
  
  first_term = (((b/l)^2)/(3*Nk))^2
  
  second_term_numerator = 2*epsilon_k*((b/l)^2)
  second_term_denominator = 5 + 2*alpha*epsilon_k*((b/l)^2)
  second_term = (second_term_numerator/second_term_denominator)^2
  
  sqrt(first_term + second_term)
}

velocity <- function(Nbp, V_cm, gel_percentage) {
  # mu_0 is maximum mobility of free dsDNA in buffer
  mu_0 = 3.4e-8
  
  # bit of a hack
  pore = 350 - 100*gel_percentage
  
  Ncrit = pore * 4
  
  if (Nbp >= Ncrit) {
    # proportion of maximum_mobility
    mobility = mu2mu0(Nbp, V_cm, pore) * mu_0
  } else {
    rcrit = mu2mu0(Ncrit, V_cm, pore)
    # maximum speed in gel relative to maximum speed in buffer
    rmax = rcrit + 2*pore/Ncrit      
    rNbp = ((rcrit - rmax)/Ncrit)*Nbp + rmax
    
    mobility = rNbp * mu_0
  }
  
  # res is in m^2 per V per second
  # convert to cm^2 per V per minute
  mobility = mobility*(100^2)*60
  
  # convert to velocity
  mobility * V_cm
}

make_lane = function(positions, intensities, width, height, multiplier=100) {
  p2 = as.integer(positions * multiplier)
  lane = matrix(0,ncol=width,nrow=height)
  for (i in 1:length(p2)) {
    dens = dgamma(seq(0,20), shape=2, rate=0.4)
    dens = dens * (intensities[i]/max(dens))
    dens = rev(dens)
    # only place bands that are still in gel
    if (p2[i] <= height) {
      lane[(p2[i]-length(dens)+1):p2[i],] = lane[(p2[i]-length(dens)+1):p2[i],] + dens
    }
  }
  lane
}

gel <- function(percentage, width, height, wellsize, spacing, margin) {
  nWells = (width - 2*margin) / (wellsize + spacing)
  if ((nWells - floor(nWells)) > wellsize) {
    nWells = floor(nWells) + 1
  } else {
    nWells = floor(nWells)
  }
  g = new.env(parent=globalenv())
  g$percentage = percentage
  g$nWells=nWells
  g$width=width
  g$height=height
  g$wellsize=wellsize
  g$spacing=spacing
  g$margin=margin
  g$sizes=list()
  g$positions=list()
  g$concentrations=list()
  for (i in 1:nWells) {
    g$sizes[[i]] = NA
    g$positions[[i]] = NA
    g$concentrations[[i]] = NA
  }
  class(g) <- "gel"
  g
}

loadGel <- function(gel, lane, sizes, concentrations=NULL) {
  if (is.null(concentrations)) {
    concentrations = rep(100, length(sizes))
  }
  if ( (lane < 1) || (lane > gel$nWells)) {
    stop("lane number not in gel range")
  }
  gel$sizes[[lane]] = sizes
  gel$concentrations[[lane]] = concentrations
}

runGel <- function(gel, time, voltage) {
  for (i in 1:gel$nWells) {
    if (!is.na(gel$sizes[[i]][1])) {
      gel$positions[[i]] = sapply(gel$sizes[[i]], function(x) {velocity(x, voltage, gel$percentage) } ) * time
    }
  }
}

resetGel <- function(gel) {
  for (i in gel$nWells) {
    g$positions[[i]] = NA
  }
}

plotGel <- function(gel, labels=NULL, scale=1) {
  if (is.null(labels)) {
    labels = c("ladder", 1:(gel$nWells-1))
  }
  
  # work with widths in 1 mm increments
  margin_px = floor(gel$margin * 10)
  well_px = floor(gel$wellsize * 10)
  spacing_px = floor(gel$spacing * 10)
  
  # height of gel is 100 * height in cm
  height_multiplier = 100
  height_px = gel$height * height_multiplier
  
  # start with the left margin
  img = matrix(0, ncol=margin_px, nrow=height_px)
  
  # add lane data
  for (i in 1:gel$nWells) {
    if (!is.na(gel$positions[[i]][1])) {
      img = cbind(img, make_lane(gel$positions[[i]], gel$concentrations[[i]], well_px, height_px, height_multiplier))
    } else {
      img = cbind(img, matrix(0, ncol=well_px, nrow=height_px))
    }
    if (i != gel$nWells) {
      img = cbind(img, matrix(0, ncol=spacing_px, nrow=height_px))
    }
  }
  
  img = cbind(img, matrix(0, ncol=margin_px, nrow=height_px))
  
  # convert to inches and scale
  width_in = (gel$width * 0.393701) * scale
  height_in = (gel$height * 0.393701) * scale

  # reverse row order of pixel matrix
  img = img[nrow(img):1,]
  
  # transpose to make lanes vertical in the plot
  img = t(img)
  
  # need different options for on-screen and off-screen plotting
  if (names(dev.cur()) == getOption("device")) {
    old.par = par(pin=c(width_in,height_in))
  } else {
    old.par = par(pin=c(width_in,height_in), omi=c(0,0,0,0), mar=c(5.1, 4.1, 0.1, 0.1))
  }
  image(img, col=gray(20:100/100), xaxt='n', yaxt='n')
  # add a scale on the y axis
  # assume the first lane is the ladder
  ylabels = gel$sizes[[1]]
  ypositions = 1 - (gel$positions[[1]] * 100)/(height_px)
  axis(2, at=ypositions, labels=ylabels, tick=F, las=2, cex.axis=0.6, pos=0.025)
  xfirst = gel$margin + 0.5 * gel$wellsize
  xstep = gel$wellsize + gel$spacing
  xpositions = (seq(from=xfirst, by=xstep, length.out=gel$nWells)) / gel$width
  axis(1, at=xpositions[1:length(labels)], labels=labels, tick=F, las=2)
  par(old.par)
}

ladder = c(50, 100, 200, 300, 400, 500, 750, 1000, 1400, 1550, 2000)
ladconc= c(50,  50,  50,  50,  50, 100,  75,  125,   75,   75,  125)

testit <- function() {
  mygel = gel(2.0, 12, 12, 0.5, 0.2, 0.5)
  loadGel(mygel, 1, ladder, ladconc)
  for (i in 1:15) {
    loadGel(mygel, i+1, c(50*i, 100*i), c(100, 100))
  }
  runGel(mygel, 90, 80/12)
  #png(file="~/Downloads/test.png", width=3, height=3, res=300, units="in", pointsize=6)
  plotGel(mygel)
  #dev.off()
}

