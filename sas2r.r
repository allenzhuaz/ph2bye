# After fun the 'sas2r.sas' file 
library(Hmisc)
mydata <- sasxport.get("c:/mydata.xpt")

# An easier way using 'haven' package
library(haven)
mydata <- read_sas("<path to your SAS file>")
