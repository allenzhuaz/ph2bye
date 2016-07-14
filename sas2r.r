# After fun the 'sas2r.sas' file 
library(Hmisc)
mydata <- sasxport.get("c:/mydata.xpt")

# An easier way using 'haven' package, the import data type is dataframe
# so we can use 'dplyr' to manipulate the data.
library(haven)
mydata <- read_sas("c://mydata.sas7bdat")

