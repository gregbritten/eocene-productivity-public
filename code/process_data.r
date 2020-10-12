accum <- read.csv('data/DSDP_596_Fish_Accumulation_siteid_1_132.csv')
size  <- read.csv('data/596_size_structure.csv',skip=1)

##--ATTACH AGE TO SIZE DATA--#####
ids <- unique(size$SiteID)
size$age=size$accum <- NA
for(i in 1:length(ids)){
  size$age[size$SiteID==ids[i]] <- accum$age[accum$SiteId==ids[i]]
  size$accum[size$SiteID==ids[i]] <- accum$ich_accum[accum$SiteId==ids[i]]
}

