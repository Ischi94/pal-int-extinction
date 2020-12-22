# load libraries
library(here)
library(divDyn)

# load isotope data
isot <- read.csv(file = here("data/Prokoph & Veizer 2015.csv"), header = TRUE)

# and stages
data(stages)

# This omits all (ant)arctic, and southern temperate for PF
isot <- subset(isot, isot$climate %in% c("trop", "temp")) 

# Detrending based on Veizer & Prokoph 2015
base <- -0.00003 * isot$gts2012 ^ 2 + 0.0046 * isot$gts2012
coriso <- isot[, 3] - base

# Time binning the data based on median, MAD and n of samples, for tropical and temperate separate and together
isostage <- numeric()
isomad <- numeric()
isostagetrop <- numeric()
isostagetemp <- numeric()
isostagetropn <- numeric()
isostagetempn <- numeric()
for (i in 1:length(stages[, 1])) {
  isostage[i] <-
    median(coriso[isot$gts2012 < stages$bottom[i] &
                    isot$gts2012 > stages$top[i]], na.rm = TRUE)
  isostagetrop[i] <-
    median(coriso[isot$gts2012 < stages$bottom[i] &
                    isot$gts2012 > stages$top[i] &
                    isot$climate == "trop"], na.rm = TRUE)
  isostagetropn[i] <-
    sum(
      isot$gts2012 < stages$bottom[i] &
        isot$gts2012 > stages$top[i] & isot$climate == "trop",
      na.rm = TRUE
    )
  isomad[i] <-
    mad(coriso[isot$gts2012 < stages$bottom[i] &
                 isot$gts2012 > stages$top[i]], na.rm = TRUE)
  isostagetemp[i] <-
    median(coriso[isot$gts2012 < stages$bottom[i] &
                    isot$gts2012 > stages$top[i] &
                    isot$climate == "temp"], na.rm = TRUE)
  isostagetempn[i] <-
    sum(
      isot$gts2012 < stages$bottom[i] &
        isot$gts2012 > stages$top[i] & isot$climate == "temp",
      na.rm = TRUE
    )
}


N.exp.mig2 <- isostage[14:94] # Original, with temperate
N.exp.mig3 <- N.exp.mig2 # Keep original
N.exp.mig3[39] <-
  (N.exp.mig3[38] + N.exp.mig3[40]) / 2  # Interpolating the missing Induan value
N.exp.mig3 <-
  16.9 - (4 * (N.exp.mig3 - 0.27))  # Adjusting the d18O to approximate degrees celsius based on Visser, Thunell, and Stott (2003; recommended in Veizer & Prokoph, 2015)

# Use original to calculate the tropical-only degrees celcius
N.exp.mig2 <- isostagetrop[14:94] # *4
N.exp.mig2[39] <- (N.exp.mig2[38] + N.exp.mig2[40]) / 2
N.exp.migtrop <- 16.9 - (4 * (N.exp.mig2 - 0.27))

# There is a large interval (all of the Jurassic, I think) when there were few tropical data, so  there we stitch in the tropical&temperate values. See Reddin et al 2018 for further justification
N.exp.mig3 <-
  c(N.exp.migtrop[1:45], N.exp.mig3[46:60], N.exp.migtrop[61:81])

# we use N.exp.mig3 for further processing, i.e. it contains the same information as TimeSeriesUsed.csv

