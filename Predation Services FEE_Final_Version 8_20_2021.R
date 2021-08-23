#### Simple lotka-voltera

install.packages("deSolve")
library(deSolve)


# Simple predator-prey model to simulate abundance 
# prey is density-depnedent with K


K = 1000	# carrying capacity of prey


# Use this set of parameters to get stability: long-lived pred and prey
# values and ratio-dependent equation from Eberhardt 1997: http://www.nrcresearchpress.com/doi/pdf/10.1139/z97-184

r1 = 0.38				# intrinsic (max) rate of increase of prey
r2 = 0.48				# intrinsic (max) rate of increase of pred
a = 	0.029			# equilibrium ratio (pred per prey at equilibrium)
c = 	6.9				# predation rate (number of prey taken per wolf per year)

N.V.start = K/2						# number of prey at start
N.P.start = round(0.001*N.V.start)	# there are 1 predator for every 100 prey at start


tmax = 30
output = matrix(nrow=tmax, ncol=5)
P = N.P.start
V = N.V.start
V.alt = N.V.start

for(t in 1:tmax){								# Simulate predator prey dynamics
	
	V.alt = V.alt + V.alt*r1*(1-(V.alt/K)) 		# logistic growth, no pred
	V = V + r1*V*(1-V/K)-c*P					# Prey N, this time step
	P = P + r2*P*(1-P/(a*V))					# Pred N, this time step  
	if(V <=2){V = 2}
	if(P <=2){P = 2}
	
	output[t,1] = t									# save time step
	output[t,2] = V									# save prey pop size
	output[t,3] = P									# save predator pop size
	output[t,4] = c*V								# Prey killed by predators
	output[t,5] = V.alt								# V alt, -DD, but no predation at all
	}

output = data.frame(output)
colnames(output) = c("t", "V", "P", 					# col 1-3
					"V.P.kill", "V.alt")				# col 4-5

plot(1:tmax, output$V.alt, type = "l", ylim=c(0, max(output$V.alt)))
points(1:tmax, output$V, type = "l", col="black", pch=19)
points(1:tmax, output$P, type = "l", col="red")

##################################################################################
### Now experiment with different relationships between density and cost/benefits

# Figure 2. Basic functional relationships of direct costs/benefits

N = 1:1000							# abundance of an animal
Nmax=1000							# carrying capacity
N.c.max = 0.9						# max cost per capita
N.b.max = 1							# max benefit per capita
N.c.min = 0.1*N.c.max				# min cost per capita
N.b.min = 0.1*N.b.max				# min benefit per capita		
N.c.k = log(N.c.max/N.c.min)*(N/Nmax) # exponential growth rate of costs
N.b.k = log(N.b.min/N.b.max)*(N/Nmax)	# exponential decline rate of benefits

N.c.pc.con = N.c.max				# per-capita (marginal) cost is constant
N.b.pc.con = N.b.max				# per-capita (pc) benefit is constant
N.c.pc.lin = N.c.max*(N/Nmax)		# pc cost is larger as pop grows (linear)
N.b.pc.lin = N.b.max*(1-N/Nmax)		# pc benefit is smalle as pop grows (linear)
N.c.pc.exp = N.c.min*exp(N.c.k)	# pc cost grows faster as pop grows
N.b.pc.exp = N.b.max*exp(N.b.k)	# pc benefit declines faster as pop grows
#N.c.pc.exp = N.c.min*N.c.k*N		# pc cost grows faster as pop grows
#N.b.pc.exp = N.b.max*N.b.k*N		# pc benefit declines faster as pop grows

N.c.tc.con = N.c.pc.con*N			# total costs from constant marginals
N.b.tc.con = N.b.pc.con*N			# total benefits from constant marginals
N.c.tc.lin = N.c.pc.lin*N			# linearly increasing total costs 
N.b.tc.lin = N.b.pc.lin*N			# linearly decreasing total benefits
N.c.tc.exp = N.c.pc.exp*N			# costs increasing exponentially
N.b.tc.exp = N.b.pc.exp*N			# benefits decreasing exponentially



####### Now, calculate indirect effects, via killing and eating prey...

V.c.max = 0.9			# maximum cost of prey per capita (pc)
V.b.max = 1.0			# maximum benefit of prey pc
P.c.max = 0.9			# maximum cost of predators pc	
P.b.max = 1.0			# maximum benefit of predators pc

V.c.min = 0.001*V.c.max	# min cost per capita
V.b.min = 0.001*V.b.max	# min benefit per capita		
P.c.min = 0.001*P.c.max	# min cost per capita
P.b.min = 0.001*P.b.max	# min benefit per capita	


#### First: Exponential/Asymtotic relationships of marginals (most complex)
	
V.c.k = log(V.c.max/V.c.min)*(1/K) 	# exponential growth rate of costs
V.b.k = log(V.b.min/V.b.max)*(1/K)	# exponential decline rate of benefits
P.c.k = log(P.c.max/P.c.min)*(1/K) 	# exponential growth rate of costs
P.b.k = log(P.b.min/P.b.max)*(1/K)	# exponential decline rate of benefits

out.asym.both = output				# store a version of the output to work with...
V = output$V						# prey population size, with predators
V.alt = output$V.alt				# prey population size, no predators
P = output$P						# predator population size

# C/B of prey, WITH predators in system... multiply pc C/B by exp(populaton size)
out.asym.both$V.c.pc = V.c.min*exp(V.c.k*V)			# V pc cost grows faster as pop grows
out.asym.both$V.b.pc = V.b.max*exp(V.b.k*V)			# V pc benefit declines faster as pop grows
out.asym.both$V.c.tot  = V.c.min*exp(V.c.k*V)*V		# V cost tot
out.asym.both$V.b.tot = V.b.max*exp(V.b.k*V)*V		# V beni tot

# C/B of prey, WITHOUT predators in system
out.asym.both$V.c.pc.alt = V.c.min*exp(V.c.k*V.alt)	# V alt cost per cap
out.asym.both$V.b.pc.alt = V.b.max*exp(V.b.k*V.alt)	# V alt beni per cap
out.asym.both$V.c.tot.alt = V.c.min*exp(V.c.k*V.alt)*V.alt	# V alt cost tot, +DD, no pred
out.asym.both$V.b.tot.alt = V.b.max*exp(V.b.k*V.alt)*V.alt	# V alt beni tot, -DD, no pred.

# Direct C/B of predators
out.asym.both$P.c.pc.dir = P.c.min*exp(P.c.k*P)			# P pc cost grows faster as pop grows
out.asym.both$P.b.pc.dir = P.b.max*exp(P.b.k*P)			# P pc benefit declines faster as pop grows
out.asym.both$P.c.tot.dir  = P.c.min*exp(P.c.k*P)*P		# P cost tot
out.asym.both$P.b.tot.dir = P.b.max*exp(P.b.k*V)*P		# P beni tot

# Indirect C/B of predators, calculated as "prey that do not exist this year compared to no-pred pop"
out.asym.both$P.b.tot.indir = out.asym.both$V.c.tot.alt-out.asym.both$V.c.tot		# ben = avoided cost = how many more prey costs in no-pred world?
out.asym.both$P.c.tot.indir = out.asym.both$V.b.tot.alt-out.asym.both$V.b.tot		# cost = foregone beni = how many fewer prey ben in with-pred world?

# Totabl C/B of predators, calculated as "prey that do not exist this year compared to no-pred pop"
out.asym.both$P.b.tot = out.asym.both$P.b.tot.dir + out.asym.both$P.b.tot.indir 
out.asym.both$P.c.tot = out.asym.both$P.c.tot.dir + out.asym.both$P.c.tot.indir 


###############################################
#### Next: Linear relationships of marginals

out.lin.both = output

# C/B of prey, with predators in system
out.lin.both$V.c.pc = V.c.max*(V/K)			# V pc cost is larger as pop grows
out.lin.both$V.b.pc = V.b.max*(1-(V/K))		# V pc benefit is smaller as pop grows
out.lin.both$V.c.tot  = out.lin.both$V.c.pc*V	# V cost tot
out.lin.both$V.b.tot = out.lin.both$V.b.pc*V		# V beni tot

# C/B of prey, without predators in system
out.lin.both$V.c.pc.alt = V.c.max*(V.alt/K)		# V alt cost per cap
out.lin.both$V.b.pc.alt = V.b.max*(1-(V.alt/K))	# V alt beni per cap
out.lin.both$V.c.tot.alt = out.lin.both$V.c.pc.alt*V.alt	# V alt cost tot, +DD, no pred
out.lin.both$V.b.tot.alt = out.lin.both$V.b.pc.alt*V.alt	# V alt beni tot, -DD, no pred.

# Direct C/B of predators
out.lin.both$P.c.pc.dir = P.c.max*(P/(a*V))			# P pc cost 
out.lin.both$P.b.pc.dir = P.b.max*(1-P/(a*V))		# P pc benefit
out.lin.both$P.c.tot.dir  = out.lin.both$P.c.pc.dir*P	# P cost tot
out.lin.both$P.b.tot.dir = out.lin.both$P.b.pc.dir*P		# P beni tot

# Indirect C/B of predators, calculted from "prey eaten this year"

# Indirect C/B of predators, calculated as "prey that do not exist this year compared to no-pred pop"
out.lin.both$P.b.tot.indir = out.lin.both$V.c.tot.alt-out.lin.both$V.c.tot		# ben = avoided cost = how many more prey costs in no-pred world?
out.lin.both$P.c.tot.indir = out.lin.both$V.b.tot.alt-out.lin.both$V.b.tot		# cost = foregone beni = how many fewer prey ben in with-pred world?

# Total C/B of predators, calculated as "prey that do not exist this year compared to no-pred pop"
out.lin.both$P.b.tot = out.lin.both$P.b.tot.dir + out.lin.both$P.b.tot.indir 
out.lin.both$P.c.tot = out.lin.both$P.c.tot.dir + out.lin.both$P.c.tot.indir 


###############################################
#### Constant marginals

out.con.both = output

# C/B of prey, with predators in system
out.con.both$V.c.pc = V.c.max			# V pc cost is constant
out.con.both$V.b.pc = V.b.max			# V pc benefit is constant
out.con.both$V.c.tot  = out.con.both$V.c.pc*V	# V cost tot
out.con.both$V.b.tot = out.con.both$V.b.pc*V		# V beni tot

# C/B of prey, without predators in system
out.con.both$V.c.pc.alt = V.c.max		# V alt cost per cap
out.con.both$V.b.pc.alt = V.b.max	# V alt beni per cap
out.con.both$V.c.tot.alt = out.con.both$V.c.pc.alt*V.alt	# V alt cost tot, +DD, no pred
out.con.both$V.b.tot.alt = out.con.both$V.b.pc.alt*V.alt	# V alt beni tot, -DD, no pred.

# Direct C/B of predators
out.con.both$P.c.pc.dir = P.c.max			# P pc cost 
out.con.both$P.b.pc.dir = P.b.max		# P pc benefit
out.con.both$P.c.tot.dir  = out.con.both$P.c.pc.dir*P	# P cost tot
out.con.both$P.b.tot.dir = out.con.both$P.b.pc.dir*P		# P beni tot

# Indirect C/B of predators, calculted from "prey eaten this year"

# Indirect C/B of predators, calculated as "prey that do not exist this year compared to no-pred pop"
out.con.both$P.b.tot.indir = out.con.both$V.c.tot.alt-out.con.both$V.c.tot		# ben = avoided cost = how many more prey costs in no-pred world?
out.con.both$P.c.tot.indir = out.con.both$V.b.tot.alt-out.con.both$V.b.tot		# cost = foregone beni = how many fewer prey ben in with-pred world?

# Total C/B of predators, calculated as "prey that do not exist this year compared to no-pred pop"
out.con.both$P.b.tot = out.con.both$P.b.tot.dir + out.con.both$P.b.tot.indir 
out.con.both$P.c.tot = out.con.both$P.c.tot.dir + out.con.both$P.c.tot.indir 


###################### Figures ###############################################

### Figure 2. Basic c/b marginal and total relationships
##### Marginal and total costs and benefits of a species, direct effects only #########
par(mfcol=c(3,2), las=1)
col.cost = "red"
col.ben = "dodgerblue"
lwd=1.5

### Marginal (per capita) values in this row
# constant marginal values
plot(N, rep(N.b.max, length(N)), type="l", col=col.ben, xlab="N", ylab="", ylim=c(0,1), lwd=lwd)
points(N,rep(N.c.max, length(N)), type="l", col=col.cost, lwd=lwd, lty=1)
mtext("", at = -20, line = 2)
legend(100, 0.5, legend=c("Cost", "Benefit"), col=c(col.cost, col.ben), lwd=lwd, lty=1, box.col="transparent")

# linearly increasinhg/decreasing marginal values
plot(N, N.b.pc.lin, type="l", col=col.ben, xlab="N", ylab="", ylim=c(0,1), lwd=lwd)
points(N,N.c.pc.lin, type="l", col=col.cost, lwd=lwd, lty=1)
mtext("", at = -20, line = 2)

# exponentially increasing/decreasing marginal values
plot(N, N.b.pc.exp, type="l", col=col.ben, xlab="N", ylab="", ylim=c(0,1), lwd=lwd)
points(N, N.c.pc.exp, type="l", col=col.cost, lwd=lwd, lty=1)
mtext("", at = -20, line = 2)


### total values in this row
plot(N, N.b.tc.con, type = "l", col=col.ben, xlab="N", ylab="", ylim=c(0,1000), lwd=lwd)
points(N, N.c.tc.con, type="l", col=col.cost, lwd=lwd, lty=1)
mtext("", at = -20, line = 2)

plot(N, N.b.tc.lin, col=col.ben, type = "l", xlab="N", ylab="", ylim=c(0,1000), lwd=lwd)
points(N, N.c.tc.lin, type="l", col=col.cost, lwd=lwd, lty=1)
mtext("", at = -20, line = 2)

plot(N, N.b.tc.exp, col=col.ben, type = "l", xlab="N", ylab="", ylim=c(0, 1000), lwd=lwd)
points(N, N.c.tc.exp, type="l", col=col.cost, lwd=lwd, lty=1)
mtext("", at = -20, line = 2)



############ Figure 3 Abundance ##########

par(mfrow=c(1,1), las=1)
ylim=c(0,1000)
xaxp=c(0,30,3)

dat = out.con.both
### Abundance of prey w/ and w/out pred, and abundance of pred. ### 
plot(dat$t, dat$V.alt, type = "l",				# V abundance, no pred
	xlab="Time", ylab="Abundance (N)", main = "",
	col = "grey50",
	ylim=ylim, lty=2, lwd=lwd, xaxp=xaxp)
lines(dat$t, dat$V, 								# V abundance, w/pred
	type = "l", xlab="", ylab="", col="grey50", lwd=lwd)
lines(dat$t, dat$P, 								# P abundance
	type = "l", xlab="", ylab="", col="grey50", lwd=lwd, lty=4)
#abline(h = 0, lty = 3, col = "grey70", lwd=lwd)
legend(0, 400, legend= c("Prey, no predation", "Prey, with predation", "Predators"),
	col = c("grey50", "grey50", "grey50"), 
	lty = c(2, 1, 4), lwd=lwd, box.col = "transparent", cex=0.8)



#####Figure 4: indirect pred beni = prey C/B with and without predation, 
#all 3 marginal val relationships

par(mfrow=c(3,2), las=1, mar = c(3,4,4,2))

ylim=c(-500, 1000)
lwd=1.5
xaxp=c(0, 30, 3)
col.prey.c = "coral3"
col.prey.b = "aquamarine4"
col.pred.c = "red"
col.pred.b = "dodgerblue"
#yaxp=c(y1=0,y2=1000, n=4)


# First row: constant marginals
dat = out.con.both

## indirect benefit of pred = prevented prey costs
plot(dat$t, dat$V.c.tot.alt, type = "l", 		# V cost, no pred
	xlab="", ylab="", 
	col = col.prey.c, lty = 2,
	ylim=ylim, lwd=lwd, xaxp=xaxp)
lines(dat$t, dat$V.c.tot, 						# V cost, pred 
	type = "l", col=col.prey.c, lwd=lwd)
lines(dat$t, dat$P.b.tot.indir, 				# P indir beni 
	type = "l", col=col.pred.b, lwd=lwd, lty=4)
abline(h = 0, lty = 3, col = "grey80", lwd=lwd)
legend(0, 0, legend= c("Prey cost w/out pred", "Prey cost w/pred", "Pred indirect benefit", "(avoided prey cost)"),
	col = c(col.prey.c, col.prey.c, col.pred.b, "transparent"), lty=c(2, 1,4, 1), lwd=lwd,box.col = "transparent", cex=0.8)
mtext("a", at = 0, line = 2)


## indirect cost of pred = lost prey benefits
plot(dat$t, dat$V.b.tot.alt, type = "l", 		# V ben, no pred
	xlab="", ylab="", 
	col = col.prey.b, lty = 2,
	ylim=ylim, lwd=lwd, xaxp=xaxp)
lines(dat$t, dat$V.b.tot, 						# V ben, pred
	type = "l", col=col.prey.b, lwd=lwd)
lines(dat$t, dat$P.c.tot.indir, 				# P indir beni
	type = "l", col=col.pred.c, lwd=lwd, lty=4)
abline(h = 0, lty = 3, col = "grey80", lwd=lwd)
legend(0, 0, legend= c("Prey ben. w/out pred", "Prey ben. w/pred", "Pred indirect cost", "(lost prey benefit)"),
	col = c(col.prey.b, col.prey.b, col.pred.c, "transparent"), lty=c(2, 1, 4, 1), lwd=lwd,box.col = "transparent", cex=0.8)
mtext("b", at = -0, line = 2)

	
# Second row: linear marginals
dat = out.lin.both

## indirect benefit of pred = prevented prey costs
plot(dat$t, dat$V.c.tot.alt, type = "l", 		# V cost, no pred
	xlab="", ylab="", 
	col = col.prey.c, lty = 2,
	ylim=ylim, lwd=lwd, xaxp=xaxp)
abline(h = 0, lty = 3, col = "grey80", lwd=lwd)
lines(dat$t, dat$V.c.tot, 						# V cost, pred 
	type = "l", col=col.prey.c, lwd=lwd)
lines(dat$t, dat$P.b.tot.indir, 				# P indir beni 
	type = "l", col=col.pred.b, lwd=lwd, lty=4)
mtext("c", at = 0, line = 2)


## indirect cost of pred = lost prey benefits
plot(dat$t, dat$V.b.tot.alt, type = "l", 		# V ben, no pred
	xlab="", ylab="", 
	col = col.prey.b, lty = 2,
	ylim=ylim, lwd=lwd, xaxp=xaxp)
lines(dat$t, dat$V.b.tot, 						# V ben, pred
	type = "l", col=col.prey.b, lwd=lwd)
lines(dat$t, dat$P.c.tot.indir, 				# P indir beni
	type = "l", col=col.pred.c, lwd=lwd, lty=4)
abline(h = 0, lty = 3, col = "grey80", lwd=lwd)
mtext("d", at = 0, line = 2)


# # Third row: asymtotic marginals
dat = out.asym.both

## indirect benefit of pred = prevented prey costs
plot(dat$t, dat$V.c.tot.alt, type = "l", 		# V cost, no pred
	xlab="", ylab="Value ($)", 
	col = col.prey.c, lty = 2,
	ylim=ylim, lwd=lwd, xaxp=xaxp)
lines(dat$t, dat$V.c.tot, 						# V cost, pred 
	type = "l", col=col.prey.c, lwd=lwd)
lines(dat$t, dat$P.b.tot.indir, 				# P indir beni 
	type = "l", col=col.pred.b, lwd=lwd, lty=4)
abline(h = 0, lty = 3, col = "grey80", lwd=lwd)
mtext("e", at = 0, line = 2)


## indirect cost of pred = lost prey benefits
plot(dat$t, dat$V.b.tot.alt, type = "l", 		# V ben, no pred
	xlab="", ylab="", 
	col = col.prey.b, lty = 2,
	ylim=ylim, lwd=lwd, xaxp=xaxp)
lines(dat$t, dat$V.b.tot, 						# V ben, pred
	type = "l", col=col.prey.b, lwd=lwd)
lines(dat$t, dat$P.c.tot.indir, 				# P indir beni
	type = "l", col=col.pred.c, lwd=lwd, lty=4)
abline(h = 0, lty = 3, col = "grey80", lwd=lwd)
mtext("f", at = 0, line = 2)



############ Figure 5 ####################
# Predator C/B figure: how do we get to predator net costs/benefits? Example

par(mfcol=c(2,2), las=1, mar = c(3,4,4,2))

#legend(0, 400, legend= c("Predator net benefit"), col = "black", lwd = lwd, lty = c(1.5), box.col = "transparent", cex=0.8)

dat=out.con.both

# Direct c/b of predators only- differet y scale
plot(dat$t, dat$P.b.tot.dir, type = "l", 		# Pred indirect benefit
	xlab="", ylab="", 
	col = col.pred.b, lty = 1,
	ylim=c(-10,30), lwd=lwd, xaxp=xaxp)
abline(h = 0, lty = 3, col = "grey80", lwd=lwd)
lines(dat$t, dat$P.c.tot.dir, 						# V cost, pred 
	type = "l", col=col.pred.c, lwd=lwd)
legend(0, 30, legend= c("Predator direct cost", "Predator direct benefit"), col = c(col.pred.c, col.pred.b), lty=1, lwd=lwd, box.col = "transparent", cex=0.8)
mtext("a", at = 0, line = 2)


### Direct and indirect c/b of predators
plot(dat$t, dat$P.b.tot.dir*2, 						# P benefi over t
	type = "l", xlab="", ylab="", 
	col = "dodgerblue", lwd = lwd,
	ylim=c(-250, 1000), xaxp=xaxp)
abline(h = 0, lty = 3, col = "grey80", lwd=lwd)
lines(dat$t, dat$P.b.tot.indir, 					# P indir beni over t
	type = "l", col="dodgerblue", lwd=lwd, lty=4)
lines(dat$t, dat$P.c.tot.dir, 					# P cost over t
	type = "l", col="red", lwd=lwd)
lines(dat$t, dat$P.c.tot.indir, 					# P indire cost over t
	type = "l", col="red", lwd=lwd, lty=4)
lines(dat$t, dat$P.b.tot-dat$P.c.tot, 			# P net effect
	type = "l", col="black", lwd=lwd, lty=1)
legend(0, 1000, legend= c("Predator indirect cost", "Predator indirect benefit", "Predator net value"), col = c(col.pred.c, col.pred.b, "black"), lty=c(3,3,1), lwd=lwd, box.col = "transparent", cex=0.8)
mtext("c", at = 0, line = 2)


# linear marginals
dat=out.lin.both

# Direct c/b of predators only- differet y scale
plot(dat$t, dat$P.b.tot.dir, type = "l", 		# Pred indirect benefit
	xlab="", ylab="", 
	col = col.pred.b, lty = 1,
	ylim=c(-10,30), lwd=lwd, xaxp=xaxp)
abline(h = 0, lty = 3, col = "grey80", lwd=lwd)
lines(dat$t, dat$P.c.tot.dir, 						# V cost, pred 
	type = "l", col=col.pred.c, lwd=lwd)
mtext("b", at = 0, line = 2)


### Direct and indirect c/b of predators
plot(dat$t, dat$P.b.tot.dir*2, 						# P benefi over t
	type = "l", xlab="", ylab="", 
	col = "dodgerblue", lwd = lwd,
	ylim=c(-250, 1000), xaxp=xaxp)
abline(h = 0, lty = 3, col = "grey80", lwd=lwd)
lines(dat$t, dat$P.b.tot.indir, 					# P indir beni over t
	type = "l", col="dodgerblue", lwd=lwd, lty=4)
lines(dat$t, dat$P.c.tot.dir, 					# P cost over t
	type = "l", col="red", lwd=lwd)
lines(dat$t, dat$P.c.tot.indir, 					# P indire cost over t
	type = "l", col="red", lwd=lwd, lty=4)
lines(dat$t, dat$P.b.tot-dat$P.c.tot, 			# P net effect
	type = "l", col="black", lwd=lwd, lty=1)
mtext("d", at = 0, line = 2)



######### SI Figures



#Fig SI 2

par(mfrow=c(1,2))
dat = out.asym.both

## indirect benefit of pred = prevented prey costs
plot(dat$t, dat$V.c.tot.alt, type = "l", 		# V cost, no pred
	xlab="", ylab="Value ($)", 
	col = col.prey.c, lty = 2,
	ylim=ylim, lwd=lwd, xaxp=xaxp)
lines(dat$t, dat$V.c.tot, 						# V cost, pred 
	type = "l", col=col.prey.c, lwd=lwd)
lines(dat$t, dat$P.b.tot.indir, 				# P indir beni 
	type = "l", col=col.pred.b, lwd=lwd, lty=4)
abline(h = 0, lty = 3, col = "grey80", lwd=lwd)
mtext("a", at = 0, line = 2)
legend(0, 0, legend= c("Prey cost w/out pred", "Prey cost w/pred", "Pred indirect benefit", "(avoided prey cost)"),
	col = c(col.prey.c, col.prey.c, col.pred.b, "transparent"), lty=c(2, 1,4, 1), lwd=lwd,box.col = "transparent", cex=0.8)
	mtext("a", at = 0, line = 2)


## indirect cost of pred = lost prey benefits
plot(dat$t, dat$V.b.tot.alt, type = "l", 		# V ben, no pred
	xlab="", ylab="", 
	col = col.prey.b, lty = 2,
	ylim=ylim, lwd=lwd, xaxp=xaxp)
lines(dat$t, dat$V.b.tot, 						# V ben, pred
	type = "l", col=col.prey.b, lwd=lwd)
lines(dat$t, dat$P.c.tot.indir, 				# P indir beni
	type = "l", col=col.pred.c, lwd=lwd, lty=4)
abline(h = 0, lty = 3, col = "grey80", lwd=lwd)
legend(0, 0, legend= c("Prey ben. w/out pred", "Prey ben. w/pred", "Pred indirect cost", "(lost prey benefit)"),
	col = c(col.prey.b, col.prey.b, col.pred.c, "transparent"), lty=c(2, 1, 4, 1), lwd=lwd,box.col = "transparent", cex=0.8)

mtext("b", at = 0, line = 2)



#### Asymtotic marginals

# dat=out.asym.both

# # Direct c/b of predators only- differet y scale
# plot(dat$t, dat$P.b.tot.dir, type = "l", 		# Pred indirect benefit
	# xlab="", ylab="", 
	# col = col.pred.b, lty = 1,
	# ylim=c(-10,30), lwd=lwd, xaxp=xaxp)
# abline(h = 0, lty = 3, col = "grey80", lwd=lwd)
# lines(dat$t, dat$P.c.tot.dir, 						# V cost, pred 
	# type = "l", col=col.pred.c, lwd=lwd)
# mtext("c", at = 0, line = 2)


# ### Direct and indirect c/b of predators
# plot(dat$t, dat$P.b.tot.dir*2, 						# P benefi over t
	# type = "l", xlab="", ylab="", 
	# col = "dodgerblue", lwd = lwd,
	# ylim=c(-250, 1000), xaxp=xaxp)
# abline(h = 0, lty = 3, col = "grey80", lwd=lwd)
# lines(dat$t, dat$P.b.tot.indir, 					# P indir beni over t
	# type = "l", col="dodgerblue", lwd=lwd, lty=4)
# lines(dat$t, dat$P.c.tot.dir, 					# P cost over t
	# type = "l", col="red", lwd=lwd)
# lines(dat$t, dat$P.c.tot.indir, 					# P indire cost over t
	# type = "l", col="red", lwd=lwd, lty=4)
# lines(dat$t, dat$P.b.tot-dat$P.c.tot, 			# P net effect
	# type = "l", col="black", lwd=lwd, lty=1)
# mtext("f", at = 0, line = 2)


