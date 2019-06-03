# read in the data

setwd("C:/Users/Antoine/Documents/Antoine- Massey/Thesis/C. jejuni/poultry_counts")

dat <- read.csv("poultry_counts.csv") ###encoding="utf8")
library(survival)

# convert sample column to an equivalent volume
#    50uL is 0.050 mL
#  1000uL is 1.000 mL
#  p100uL is 4.000 mL as the pellet contains all bugs in 200mL which is then buffered into 5mL, 
#then a 100uL sample taken (so equivalent to sampling 4/200 mL) (5/200=0.025, 0.1/0.025=4)

dat$Volume <- 0

levels(dat$Sample) -> samples

dat$Volume[dat$Sample == samples[1]] <- 1    # 1000 uL
dat$Volume[dat$Sample == samples[2]] <- 0.05 #   50 uL
dat$Volume[dat$Sample == samples[3]] <- 4    # 4000 uL


# generate month column

dat$Date <- as.character(dat$Date)
dat$Date <- as.Date(dat$Date, format="%d %B %Y")

dat$season <- factor(dat$Season, levels = c("Spring","Summer","Autumn", "Winter"))
levels(dat$season)

#dat$Weight <- dat$Weight.in.grams
#dat$Weight.in.grams <- NULL

# convert count to numbers, removing TNTC, UC, Dried

# TODO: Potentially TNTC is 'lots' so maybe think about how we can treat this as right censored or some such?

dat$Count <- as.numeric(as.character(dat$Count))

# convert weight to numbers (get rid of commas)
#dat$Weight <- as.numeric(gsub(",", as.character(dat$Weight), replacement="")) / 1000

# zero-inflated poisson model
library("pscl", lib.loc="~/R/win-library/3.1") 

model.9 <- zeroinfl (Count ~ Intervention + offset(log(Volume)) | Intervention, data=dat)
summary(model.9)

model.11<-zeroinfl (Count ~ Intervention*Source + offset(log(Volume)) | Intervention* Source, data=dat)
summary(model.11)

