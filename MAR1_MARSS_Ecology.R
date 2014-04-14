######## R script described in "Applying MAR(1) and MARSS to an ecological dataset" ##########
###### Supplementary material for ##########
###### "Quantifying effects of abiotic and biotic drivers on community dynamics with multivariate autoregressive (MAR) models"########
#######Hampton, Holmes, Scheef, Scheuerell, Katz, Pendleton, Ward - submitted to Ecology 18 July 2013######
##### A PDF of the exercise "Applying MAR(1) and MARSS to an ecological dataset" is available, to place this R script in context, and can be used as a self-tutorial #############

######## MAR1 R Package ########

# Install MAR1 package

install.packages("MAR1", repos="http://cran.rstudio.com", dependencies=TRUE)

# Load package

library(MAR1)

######## Prepare the dataset ########

# Load and view L4 dataset included in MAR1 package

data(L4.AllDates)
head(L4.AllDates)
summary(L4.AllDates)

?L4.AllDates

# Begin transformations – average into monthly time steps

L4.byMonth <- prepare.data(data=L4.AllDates, increment="month")

summary(L4.byMonth)
dim(L4.byMonth)

# See all arguments and help file for prepare.data function

formals(prepare.data)
?prepare.data

# Fill gaps with linear interpolation, perform log and z-score transformations

L4.mar1 <- prepare.data(data=L4.AllDates, increment="month", fill.gap=1, replace.0s="rand.half", log=T, z.method="standard")

# Look at summary of transformed dataset

summary(L4.mar1)

# Create another dataset using the de-seasoning z-score method

L4.mar2 <- prepare.data(data=L4.AllDates, increment="month", fill.gap=1, replace.0s="rand.half", log=T, z.method="deseason")

summary(L4.mar2)

######## Build the MAR(1) model ########

# Look at the arguments available in run.mar and view help file

formals(run.mar)
?run.mar

# Initiate MAR(1) analysis using pop-up windows to select variables and set interaction restrictions
	# Set cnidarian, chaetognath, calanoid.lg, calanoid.sm, diatom, and dino as variates and surface.temp as a covariate
	# Do not set any restrictions

run1 <- run.mar(data=L4.mar1)

# Same model can be built without using the pop-up windows

myvar <- c(0,0,1,0,1,0,0,0,0,1,1,0,0,0,1,1,0,0,0,0,2)
myres <- matrix(0.5, nrow=length(which(myvar==1)), ncol=length(which(myvar!=0)))

run1 <- run.mar(data=L4.mar1, variables=myvar, restrictions=myres)

# Summarize and plot main results

run1
summary(run1)
plot(run1)

# Build another model using the de-seasoned dataset

run2<-run.mar(data=L4.mar2, variables=run1, restrictions=run1)

# Compare results of 1st and 2nd model in a plot

plot(run1, run2, legend=T)

# Run a 3rd model – this time you can toggle off one of the biologically implausible interactions (cnidarians on chaetognaths)

run3 <- run.mar(data=L4.mar2, variables=run2)

# Same model can be built without using the pop-up window

myres[2,1]<-0

run3 <- run.mar(data=L4.mar2, variables=run2, restrictions=myres)

# Compare results of all three models in a plot

plot(run1, run2, run3, legend=T)

# For top model (lowest AIC) from 3rd model, compare to other top models in a plot

dev.new()
plot(run3$top.bestfit)

# Compare AIC values of top models using a histogram

hist(run3$top.bestfit)


######## Building a state-space MAR model (MARSS) ########

# Install and load the MARSS package

install.packages("MARSS", repos="http://cran.rstudio.com", dependencies=TRUE)

library(MARSS)

# Compare non-state-space MAR to a state-space MAR

ss.run3 <- ss.mar1(L4.mar2, MAR.obj=run3)

ss.run3$B
run3$bestfit$B

ss.run3.plot <- list(restrictions.set=run1$restrictions.set, bestfit=list(B=ss.run3$B,C=ss.run3$C), bootstrap=NULL)
class(ss.run3.plot) <- "MAR"

dev.new()
plot(run3, ss.run3.plot, legend=T)

# Use a state-space model with L4 data aggregated by month with gaps not filled

L4.mar3 <- prepare.data(data=L4.AllDates, increment="month", fill.gap=0, replace.0s="rand.half", log=TRUE, z.method="deseason")

ss.run3.no.interp <- ss.mar1(L4.mar3, MAR.obj=run3)

ss.run3.no.interp$B

# Use a state-space model with a fixed observation variance

R <- diag(c(.1,.1,.1,.1,.05,.05))

ss.run3.fixed.R <- ss.mar1(L4.mar3, MAR.obj=run3, model=list(R=R))

ss.run3.fixed.R$B

# Compare non-state-space and all three state-space MAR models in a plot

ss.run3.no.interp.plot <- list(restrictions.set=run1$restrictions.set, bestfit=list(B=ss.run3.no.interp$B,C=ss.run3.no.interp$C), bootstrap=NULL)
class(ss.run3.no.interp.plot) <- "MAR"

ss.run3.fixed.R.plot <- list(restrictions.set=run1$restrictions.set, bestfit=list(B=ss.run3.fixed.R$B,C=ss.run3.fixed.R$C), bootstrap=NULL)
class(ss.run3.fixed.R.plot) <- "MAR"

dev.new()
plot(run3, ss.run3.plot, ss.run3.no.interp.plot, ss.run3.fixed.R.plot, legend=T)