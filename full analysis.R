################################
# Analyses for Dale et al. Veronica sexual system and biome shift paper
#phylogeny from Thomas et al. 2023 
# section 3 code written by Luke Liddell
# all other sections written by Esther Dale
#################################

# Section 1 reading in and preparing data
## 1.1 biome data
## 1.2 phylogeny
## 1.3 trait data

# Section 2 biogeographic analysis - biome shifts
# 2.1 model selection
# 2.2 biogeographic stochastic mapping

# section 3 trait-independent and dependent biogeographic analysis - biome effect on trait evolution

# section 4 analyses
##

library(rexpokit)
library(cladoRcpp)
library(devtools)
#devtools::install_github(repo='nmatzke/BioGeoBEARS')
library(optimx)  
library(FD)      
library(snow)    
library(parallel)
library(geiger)
library(BioGeoBEARS)
library(evaluate) 
library(dplyr)
library(tidyr)
library(ape)
library(phangorn)

extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
extdata_dir
list.files(extdata_dir)

####read in phylogeny ####
tr = read.tree("Veronica consensus.tree")

#### prepaing biome data ####
clade <- "Veronica"
#alt_lat_bands.csv can be downloaded via https://github.com/annethomas/veronica-biogeo/raw/refs/heads/main/data/distribution_data/alt_lat_bands.csv
biomes.df <- read.csv("alt_lat_bands.csv")
biomes.df$species <- gsub("Veronica_", "", biomes.df$species) #remove genus from species name
biomes.df <- biomes.df[,1:4]#subset to biome data

## add in island species missing from Anne's biome dataset ##
biomes.df[nrow(biomes.df)+1,] <- c("barkeri", 0, 0, 1) # barkeri chathams, coastal [lowland]
biomes.df[nrow(biomes.df)+1,] <- c("benthamii", 0, 1, 1) #benthamii montane, cambell and Aucklands nzpcn.org.nz [lowland and mountain]
biomes.df[nrow(biomes.df)+1,] <- c("breviracemosa", 0, 0, 1)#breviracemosa kermadecs [lowland]
biomes.df[nrow(biomes.df)+1,] <- c("chathamica", 0, 0, 1) # chathamica coastal chatthams [lowland]
biomes.df[nrow(biomes.df)+1,] <- c("dieffenbachii", 0, 0, 1) #dieffenbachii coastal or inland Chathams[lowland]
biomes.df[nrow(biomes.df)+1,] <- c("insularis", 0, 0, 1) #insularis coastal cliffs [lowland]

row.names(biomes.df) <- biomes.df$species #make row names same as tree tip labels
name.check(tr,biomes.df)

#remove species not in the phylogeny
biomes.df <- biomes.df[biomes.df$species %in% tr$tip.label,]

#add in mountain column (combine subalpine and alpine)
biomes.df$mountain <- paste(biomes.df$alpine, biomes.df$subalpine, sep="")
biomes.df$mountain <- gsub("00", 0, biomes.df$mountain)
biomes.df$mountain[as.numeric(biomes.df$mountain)>0] <-1 

#get data in correct format for biogeobears
# adding in biome info for spp not in database ie exotic species
Species <- biomes.df$species # sp names
l <- biomes.df$lowland #lowland biome
m <- biomes.df$mountain #mountain (subalpine and alpine)

lm <- paste(l,m, sep="")
biomes=as.data.frame(cbind(Species,lm)) #make dataframe of biome info for each species
header <- paste(nrow(biomes), "2", "(L M)")#first row has number of spp and no. biomes

write.table(header,file="geog_file.txt",quote=F,row.names=F,col.names = F)
write.table(biomes,file="geog_file.txt",append=T,sep=' ',quote=F,row.names=F,col.names = F) #write dataframe to text file. Don't repeat cos will be repeated in text file

################ BiogeobBEARS model selection model selection #################
trfn <- "Veronica consensus.tree"
geogfn = file="geog_file.txt" #specify geog file
moref(geogfn)

tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)

#max number of areas a spp can occupy
max_range_size = 2
numstates_from_numareas(numareas=3, maxareas=3, include_null_range=TRUE)

###############################
####### DEC ANALYSIS #########
##############################

BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$force_sparse=FALSE    # sparse=FALSE causes pathology & isn't much faster at this scale
BioGeoBEARS_run_object$speedup=TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = TRUE
BioGeoBEARS_run_object$calc_ancprobs=TRUE    # get ancestral states from optim run

#### add in time stratification in future to account for age of nz alpine #####
# Set up a time-stratified analysis
# (un-comment to use; see example files in extdata_dir,
#  and BioGeoBEARS google group posts for further hints)
BioGeoBEARS_run_object$timesfn = "time_periods.txt"
BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"

# Input the maximum range size
BioGeoBEARS_run_object$max_range_size = max_range_size

# Multicore processing if desired
BioGeoBEARS_run_object$num_cores_to_use=1

BioGeoBEARS_run_object$force_sparse=FALSE

# Give BioGeoBEARS the location of the geography text file
BioGeoBEARS_run_object$geogfn = geogfn

# Give BioGeoBEARS the location of the phylogeny Newick file
BioGeoBEARS_run_object$trfn = trfn

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
# The stratified tree is described in this table:
#BioGeoBEARS_run_object$master_table

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE

#BioGeoBEARS_run_object$BioGeoBEARS_model_object #model object
#BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table # parameters of the model

# Run this to check inputs
check_BioGeoBEARS_run(BioGeoBEARS_run_object)
# For a slow analysis, run once, then set runslow=FALSE to just
# load the saved result.
runslow = TRUE
resfn = "DEC_M0_4_v1.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res
  
  save(res, file=resfn)
  resDEC = res
} else {
  # Loads to "res"
  load(resfn)
  resDEC = res
}

###############################
##### DEC + J ANALYSIS #######
##############################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size

# Set up the stratified part
BioGeoBEARS_run_object$timesfn = "time_periods.txt"
BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"


BioGeoBEARS_run_object$speedup=TRUE          # shorcut to speed ML search
BioGeoBEARS_run_object$use_optimx = TRUE
BioGeoBEARS_run_object$num_cores_to_use=1
BioGeoBEARS_run_object$force_sparse=FALSE   
BioGeoBEARS_run_object$calc_ancprobs=TRUE    # get ancestral states from optim run

BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
# The stratified tree is described in this table:
#BioGeoBEARS_run_object$master_table

# kept as defaults
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE

# Set up DEC+J model
# Get the ML parameter values from the 2-parameter nested model
# (this will ensure that the 3-parameter model always does at least as good)
dstart = resDEC$outputs@params_table["d","est"]
estart = resDEC$outputs@params_table["e","est"]
jstart = 0.0001

# Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# Add j as a free parameter
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn = "DEC+J_M0_4_v1.Rdata"
runslow = TRUE
if (runslow)
{
  #sourceall("/Dropbox/_njm/__packages/BioGeoBEARS_setup/")
  
  res = bears_optim_run(BioGeoBEARS_run_object)
  res
  
  save(res, file=resfn)
  
  resDECj = res
} else {
  # Loads to "res"
  load(resfn)
  resDECj = res
}
#ignore warning about unused control argument and parameter bounds

########################## DEC vs DEC+J plots ####################################
pdffn = "timestratified_DEC_vs_DEC+J_4_v1.pdf"
pdf(pdffn, width=5, height=20)

#### DEC ####
analysis_titletxt =paste(clade, "BioGeoBEARS DEC NZ", sep=" ")

# Setup
results_object = resDEC
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.4, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.4, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

#####DEC + J ######
analysis_titletxt =paste(clade, "BioGeoBEARS DEC+J NZ", sep=" ")

# Setup
results_object = resDECj
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.4, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.4, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

dev.off()

##############################
##### DIVALIKE ANALYSIS ######
##############################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size

# Set up the stratified part
BioGeoBEARS_run_object$timesfn = "time_periods.txt"
BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"


BioGeoBEARS_run_object$use_optimx = TRUE
BioGeoBEARS_run_object$num_cores_to_use=1

BioGeoBEARS_run_object$force_sparse=FALSE    # sparse=FALSE causes pathology & isn't much faster at this scale
BioGeoBEARS_run_object$speedup=TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$calc_ancprobs=TRUE    # get ancestral states from optim run

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
# The stratified tree is described in this table:
#BioGeoBEARS_run_object$master_table

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE

# Set up DIVALIKE model
# Remove subset-sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

# Allow classic, widespread vicariance; all events equiprobable
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

runslow = TRUE
resfn = "DIVALIKE_M0_4_v1.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res
  
  save(res, file=resfn)
  resDIVALIKE = res
} else {
  # Loads to "res"
  load(resfn)
  resDIVALIKE = res
}
#ignore warning about unused control arguments

##############################
#### DIVALIKE+J ANALYSIS #####
##############################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size

# Set up the stratified part

BioGeoBEARS_run_object$timesfn = "time_periods.txt"
BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"

BioGeoBEARS_run_object$use_optimx = TRUE
BioGeoBEARS_run_object$num_cores_to_use=1

BioGeoBEARS_run_object$force_sparse=FALSE    # sparse=FALSE causes pathology & isn't much faster at this scale
BioGeoBEARS_run_object$speedup=TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$calc_ancprobs=TRUE    # get ancestral states from optim run

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
# The stratified tree is described in this table:
#BioGeoBEARS_run_object$master_table

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE

# Set up DIVALIKE+J model
# Get the ML parameter values from the 2-parameter nested model
# (this will ensure that the 3-parameter model always does at least as good)
dstart = resDIVALIKE$outputs@params_table["d","est"]
estart = resDIVALIKE$outputs@params_table["e","est"]
jstart = 0.0001

# Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# Remove subset-sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

# Allow classic, widespread vicariance; all events equiprobable
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

# Add jump dispersal/founder-event speciation
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

# Under DIVALIKE+J, the max of "j" should be 2, not 3 (as is default in DEC+J)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 1.99999

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn = "DIVALIKE+J_M0_4_v1.Rdata"
runslow = TRUE
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res
  
  save(res, file=resfn)
  
  resDIVALIKEj = res
} else {
  # Loads to "res"
  load(resfn)
  resDIVALIKEj = res
}
####################### plot DIVALIKE and DIVALIKE_J ##########################
pdffn = "_timestratified_DIVALIKE_vs_DIVALIKE+J_4_v1.pdf"
pdf(pdffn, width=5, height=20)

#### DIVALIKE #######
analysis_titletxt =paste(clade, "BioGeoBEARS DIVALIKE", sep=" ")

# Setup
results_object = resDIVALIKE
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.4, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.4, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

##### DIVALIKE+J###################
analysis_titletxt =paste(clade, "BioGeoBEARS DIVALIKE+J", sep=" ")

# Setup
results_object = resDIVALIKEj
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.4, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.4, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

dev.off()


####################################
####### BAYAREALIKE ANALYSIS #######
####################################

BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size

# Set up the stratified part
BioGeoBEARS_run_object$timesfn = "time_periods.txt"
BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"

# See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

BioGeoBEARS_run_object$use_optimx = TRUE
BioGeoBEARS_run_object$num_cores_to_use=1

BioGeoBEARS_run_object$force_sparse=FALSE    # sparse=FALSE causes pathology & isn't much faster at this scale
BioGeoBEARS_run_object$speedup=TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$calc_ancprobs=TRUE    # get ancestral states from optim run

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
# The stratified tree is described in this table:
#BioGeoBEARS_run_object$master_table

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE

# Set up BAYAREALIKE model
# No subset sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

# No vicariance
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

# Adjust linkage between parameters
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

# Only sympatric/range-copying (y) events allowed, and with
# exact copying (both descendants always the same size as the ancestor)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

# Check the inputs
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

runslow = TRUE
resfn = "BAYAREALIKE_M0_4_v1.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res
  
  save(res, file=resfn)
  resBAYAREALIKE = res
} else {
  # Loads to "res"
  load(resfn)
  resBAYAREALIKE = res
}
# ignore warning about unused control arguments

#######################################################
############## BAYAREALIKE+J ANALYSIS ################
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size

# Set up the stratified part
BioGeoBEARS_run_object$timesfn = "time_periods.txt"
BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
# See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

BioGeoBEARS_run_object$use_optimx = TRUE
BioGeoBEARS_run_object$num_cores_to_use=1

BioGeoBEARS_run_object$force_sparse=FALSE    # sparse=FALSE causes pathology & isn't much faster at this scale
BioGeoBEARS_run_object$speedup=TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$calc_ancprobs=TRUE    # get ancestral states from optim run

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
# The stratified tree is described in this table:
#BioGeoBEARS_run_object$master_table

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE

# Set up BAYAREALIKE+J model
# Get the ML parameter values from the 2-parameter nested model
# (this will ensure that the 3-parameter model always does at least as good)
dstart = resBAYAREALIKE$outputs@params_table["d","est"]
estart = resBAYAREALIKE$outputs@params_table["e","est"]
jstart = 0.0001

# Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# No subset sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

# No vicariance
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

# *DO* allow jump dispersal/founder-event speciation (set the starting value close to 0)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

# Under BAYAREALIKE+J, the max of "j" should be 1, not 3 (as is default in DEC+J) or 2 (as in DIVALIKE+J)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999

# Adjust linkage between parameters
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

# Only sympatric/range-copying (y) events allowed, and with
# exact copying (both descendants always the same size as the ancestor)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","min"] = 0.0000001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","max"] = 4.9999999

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"] = 0.0000001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","max"] = 4.9999999

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999

check_BioGeoBEARS_run(BioGeoBEARS_run_object)
error.check <- evaluate("check_BioGeoBEARS_run(BioGeoBEARS_run_object)")

if(grepl("STOP ERROR", error.check[[2]][1])==1){ #if stop error. searches for "STOP ERROR" in error.check object
  print("Fix params")
  BioGeoBEARS_run_object <- fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object = BioGeoBEARS_run_object) #fixes up parameter issue
}else{print("Params fine")}

resfn = "BAYAREALIKE+J_M0_v1.Rdata"
runslow = TRUE
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res
  
  save(res, file=resfn)
  
  resBAYAREALIKEj = res
} else {
  # Loads to "res"
  load(resfn)
  resBAYAREALIKEj = res
}
# ignore error about unused control arguments

#################### plot BAYAREALIKE and BAYAREALIKE+J#####################
pdffn = "timestratified_BAYAREALIKE_vs_BAYAREALIKE+J_4_v1.pdf"
pdf(pdffn, width=5, height=20)

######## BAYAREALIKE #######
analysis_titletxt =paste(clade, "BioGeoBEARS BAYAREALIKE NZ", sep=" ")

# Setup
results_object = resBAYAREALIKE
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.4, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.4, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

###### BAYAREALIKE+J ########
analysis_titletxt =paste(clade, "BioGeoBEARS BAYAREALIKE+J on NZ M0_unconstrained", sep=" ")

# Setup
results_object = resBAYAREALIKEj
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.4, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.4, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

dev.off()
#cmdstr = paste("open ", pdffn, sep="")
#system(cmdstr)


#################### Model comparison #######################################
# Set up empty tables to hold the statistical results
restable = NULL
teststable = NULL

############################ DEC vs. DEC+J###############################################
# We have to extract the log-likelihood differently, depending on the
# version of optim/optimx
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDEC)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDECj)

numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
stats

# DEC, null model for Likelihood Ratio Test (LRT)
res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resDEC, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
# DEC+J, alternative model for Likelihood Ratio Test (LRT)
res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resDECj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)

rbind(res2, res1)
tmp_tests = conditional_format_table(stats)

restable = rbind(restable, res2, res1)
teststable = rbind(teststable, tmp_tests)

################################# DIVALIKE vs. DIVALIKE+J############################################
# We have to extract the log-likelihood differently, depending on the
# version of optim/optimx
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKE)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKEj)

numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
stats

# DIVALIKE, null model for Likelihood Ratio Test (LRT)
res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKE, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
# DIVALIKE+J, alternative model for Likelihood Ratio Test (LRT)
res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKEj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)

rbind(res2, res1)
conditional_format_table(stats)

tmp_tests = conditional_format_table(stats)

restable = rbind(restable, res2, res1)
teststable = rbind(teststable, tmp_tests)

################################### BAYAREALIKE vs. BAYAREALIKE+J##################################
# We have to extract the log-likelihood differently, depending on the
# version of optim/optimx
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKE)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKEj)

numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
stats

# BAYAREALIKE, null model for Likelihood Ratio Test (LRT)
res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKE, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
# BAYAREALIKE+J, alternative model for Likelihood Ratio Test (LRT)
res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKEj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)

rbind(res2, res1)
conditional_format_table(stats)

tmp_tests = conditional_format_table(stats)

restable = rbind(restable, res2, res1)
teststable = rbind(teststable, tmp_tests)

########## RESULTS: DEC, DEC+J, DIVALIKE, DIVALIKE+J, BAYAREALIKE, BAYAREALIKE+J######################
teststable$alt = c("DEC+J", "DIVALIKE+J", "BAYAREALIKE+J")
teststable$null = c("DEC", "DIVALIKE", "BAYAREALIKE")
row.names(restable) = c("DEC", "DEC+J", "DIVALIKE", "DIVALIKE+J", "BAYAREALIKE", "BAYAREALIKE+J")

# Look at the results!!
restable
teststable

## save results ##
saveRDS(restable,file="restable_Veronica_4.Rdata")
saveRDS(teststable,file="teststable_Veronica_4.Rdata")

##################################################################
#################### BSM #########################################
##################################################################
# BSMs simplified code (so identical to code in null simulations)
# to calculate transition rates between biomes lowland (L) and mountain (M)
#stratified 100-4my lowland only, 1.9-0my lowland and mountain

clade <- "Veronica"
model_name="DEC+J" #preferred model

tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
max_range_size = 3

####### all models set up ######
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$force_sparse=FALSE    # sparse=FALSE causes pathology & isn't much faster at this scale
BioGeoBEARS_run_object$speedup=TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = TRUE

#### add in time stratification in future to acount for age of nz alpine #####
# Set up a time-stratified analysis
BioGeoBEARS_run_object$timesfn = "time_periods.txt"
BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"



# Input the maximum range size
BioGeoBEARS_run_object$max_range_size = max_range_size

# Multicore processing if desired
BioGeoBEARS_run_object$num_cores_to_use=1

# Give BioGeoBEARS the location of the geography text file
BioGeoBEARS_run_object$geogfn = geogfn

# Give BioGeoBEARS the location of the phylogeny Newick file
BioGeoBEARS_run_object$trfn = trfn

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE

#### DEC model as default - change options for other models ####

#### DIVALIKE & +J ####
if(model_name %in% c("DIVALIKE", "DIVALIKE+J")){
  # Remove subset-sympatry
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
  
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"
  
  # Allow classic, widespread vicariance; all events equiprobable
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5
}# end of DIVALIKE if

#### BAYAREALIKE & +J ####
if(model_name %in% c("BAYAREALIKE", "BAYAREALIKE+J")){
  # No subset sympatry
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
  
  # No vicariance
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0
  
  # Adjust linkage between parameters
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"
  
  # Only sympatric/range-copying (y) events allowed, and with
  # exact copying (both descendants always the same size as the ancestor)
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999
} #end of BAYREALIKE if

#### set up J part of model  if needed ####
if(model_name%in% c("DEC+J", "DIVALIKE+J", "BAYAREALIKE+J")){
  
  # Add j as a free parameter
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.001
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.001
  
  #set j min and max for DIVALIKE+J and BAYREALIKE+J
  if(model_name=="DIVALIKE+J"){
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 1.99999
  } #  end of DIVALIKE+J if
  
  if(model_name=="BAYAREALIKE+J"){
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999
    # uncomment these bits if it crashes, a BAYREALIKE+J Windows issue
    # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","min"] = 0.0000001
    # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","max"] = 4.9999999
    
    # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"] = 0.0000001
    # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","max"] = 4.9999999
    
  } #  end of BAYAREALIKE+J if
  
}#end of +J if


# Run this to check inputs
try(check_BioGeoBEARS_run(BioGeoBEARS_run_object))

error.check <- evaluate("check_BioGeoBEARS_run(BioGeoBEARS_run_object)")

if(grepl("STOP ERROR", error.check[[2]][1])==1){ #if stop error. searches for "STOP ERROR" in error.check object
  print("Fix params")
  BioGeoBEARS_run_object <- fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object = BioGeoBEARS_run_object) #fixes up parameter issue
}else{print("Params fine")}

# run model #
res = bears_optim_run(BioGeoBEARS_run_object)
resfn = paste(clade,"_", model_name,"_4_v1.Rdata", sep="")
#save(res, file=resfn)


#######################################################
# Stochastic mapping
#######################################################
clado_events_tables = NULL
ana_events_tables = NULL
lnum = 0

#######################################################
# Get the inputs for Biogeographical Stochastic Mapping
# Note: this can be slow for large state spaces and trees, since 
# the independent likelihoods for each branch are being pre-calculated
# E.g., for 10 areas, this requires calculation of a 1024x1024 matrix
# for each branch.  On a tree with ~800 tips and thus ~1600 branches, this was about 1.6 gigs
# for storage of "BSM_inputs_file.Rdata".
# Update: 2015-09-23 -- now, if you used multicore functionality for the ML analysis,
# the same settings will be used for get_inputs_for_stochastic_mapping().
#######################################################

BSM_inputs_fn = "BSM_inputs_file.Rdata"

stochastic_mapping_inputs_list = get_inputs_for_stochastic_mapping(res=res)
save(stochastic_mapping_inputs_list, file=BSM_inputs_fn)

# Check inputs (doesn't work the same on unconstr)
names(stochastic_mapping_inputs_list)
stochastic_mapping_inputs_list$phy2
stochastic_mapping_inputs_list$COO_weights_columnar
stochastic_mapping_inputs_list$unconstr
#set.seed(seed=as.numeric(Sys.time()))

# does 100 runs #
runBSMslow = TRUE
BSM_output = runBSM(res, stochastic_mapping_inputs_list=stochastic_mapping_inputs_list, maxnum_maps_to_try=150, nummaps_goal=1000, maxtries_per_branch=40000, save_after_every_try=TRUE, savedir=getwd(), seedval=12345, wait_before_save=0.01)

RES_clado_events_tables = BSM_output$RES_clado_events_tables
RES_ana_events_tables = BSM_output$RES_ana_events_tables

ana_events_tables = BSM_output$RES_ana_events_tables
# Extract BSM output

# save
resfn.clado = paste("Clado",clade, model_name,"multi.Rdata", sep="_")
resfn.ana = paste("Ana",clade, model_name, "multi.Rdata", sep="_")
saveRDS(RES_clado_events_tables, file=resfn.clado)
saveRDS(RES_ana_events_tables, file=resfn.ana)


#######################################################
# Plot all stochastic maps to PDF
#######################################################
# Setup
include_null_range = TRUE
areanames = names(tipranges@df)
areas = areanames
max_range_size = 2

states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=include_null_range)
colors_list_for_states = get_colors_for_states_list_0based(areanames=areanames, states_list_0based=states_list_0based, max_range_size=max_range_size, plot_null_range=TRUE)
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
stratified = TRUE

# Loop through the maps and plot to PDF
pdffn = paste0("lowland mountain/pdfs/",clade, model_name, "_", length(RES_clado_events_tables), "BSMs_v1.pdf")
pdf(file=pdffn, width=10, height=20)

nummaps_goal = 1000
for (i in 1:nummaps_goal)
{
  clado_events_table = RES_clado_events_tables[[i]]
  analysis_titletxt = paste0(model_name, " - Stochastic Map #", i, "/", nummaps_goal)
  plot_BSM(results_object=res, clado_events_table=clado_events_table, stratified=stratified, analysis_titletxt=analysis_titletxt, addl_params=list("j"), label.offset=0.5, plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, colors_list_for_states=colors_list_for_states, show.tip.label=TRUE, include_null_range=include_null_range)
} # END for (i in 1:nummaps_goal)

dev.off()

##################################################################
########### generate biome shift table for Veronica ##############
##################################################################

RES_ana_events_tables <- readRDS(paste("Ana",clade, model_name,"multi.Rdata", sep="_"))
RES_clado_events_tables <- readRDS(paste("Clado",clade, model_name,"multi.Rdata", sep="_"))

###### generating summary tables ###########
# add run column to all 100 tables in RES ana and clado tables list
for (i in 1:1000){
  RES_clado_events_tables[[i]]$run <- rep(i, nrow(RES_clado_events_tables[[i]]))
  RES_ana_events_tables[[i]]$run <- rep(i, nrow(RES_ana_events_tables[[i]]))
  RES_ana_events_tables[[i]] <- RES_ana_events_tables[[i]][,-c(45:46)] #the 45th columns throws an error for some clades, so will remove here
}#for loop 100

# combine the list of 100 dfs into one master df
master_RES_clado <-as.data.frame(bind_rows(RES_clado_events_tables))
master_RES_ana <-as.data.frame(bind_rows(RES_ana_events_tables))

#########################
### anagenetic shifts####
#rename columns
colnames(master_RES_ana)[44] <- "shift"
colnames(master_RES_ana)[41] <- "shift.timing"

############ Cladogenetic events ######################
#rename columns
colnames(master_RES_clado)[26] <- "shift"
colnames(master_RES_clado)[7] <- "shift.timing"
colnames(master_RES_clado)[25] <- "event_type"

#subset to rows with clado shifts
master_RES_clado <- master_RES_clado[!master_RES_clado$shift=="",]
master_RES_clado <- master_RES_clado[!is.na(master_RES_clado$event_type),]

################### combine ANA and CLADO #############################
master_RES <- rbind(master_RES_ana[,c("node", "run", "event_type", "shift", "shift.timing")], master_RES_clado[, c("node", "run", "event_type", "shift", "shift.timing")] ) #subset to useful columns of ana and clado and combine to new df: node, clade, run and shift

# save
saveRDS(master_RES, "Biome shift table with modes of speciation Veronica.Rdata")

################################################
###### simplify shift types ################
#################################################

# remove sympatric speciation
sympatric <- c("L->L,L", "M->M,M", "LM->LM,LM")
master_RES <- master_RES[!master_RES$shift %in% sympatric,]

#remove extinction and subset so only shifts involving new biomes included
master_RES <- master_RES[master_RES$event_type %in% c("d", "founder (j)"),]

##rename shifts to merge shifts with same outcome##
master_RES$shift[master_RES$shift %in% c("L->M,L")] <- "L->L,M"
master_RES$shift[master_RES$shift %in% c("M->M,L")] <- "M->L,M"

#change shift timing to numeric
master_RES$shift.timing <- as.numeric(as.character(master_RES$shift.timing))

#remove event type column cos subsetted sympatric speciation, so it wouldn't cover all speciation events
master_RES <- master_RES[,c("node", "run", "shift", "shift.timing")]

### save ###
saveRDS(master_RES, "Biome shift table expansion shifts Veronica.Rdata")

## table of biome shift types ##
all.shifts <- readRDS("Biome shift table with modes of speciation Veronica.Rdata")
all.shifts <- all.shifts[!all.shifts$shift %in% sympatric,]
all.shifts$shift[all.shifts$shift %in% c("L->M,L")] <- "L->L,M"
all.shifts$shift[all.shifts$shift %in% c("M->M,L")] <- "M->L,M"
all.shifts$shift[all.shifts$shift %in% c("LM->LM,L")] <- "LM->L,LM"
all.shifts$shift[all.shifts$shift %in% c("LM->M,L")] <- "LM->L,M"
all.shifts$shift[all.shifts$shift %in% c("LM->M,LM")] <- "LM->LM,M"

totals <- table(all.shifts$shift)/1000
LtoM <- sum(totals['L->L,M'], totals['L->M'], na.rm = T)
LtoLM <- sum(totals['L->LM']  )
LMtoM <- sum(totals['LM->M'], 0.5*totals['LM->L,M'], totals['LM->LM,M'], na.rm = T)
LMtoL <- sum(totals['LM->L'], 0.5*totals['LM->L,M'], totals['LM->L,LM'], na.rm = T)
MtoL <- sum(totals['M->L,M'], totals['M->L'], na.rm = T)
MtoLM <- sum(totals['M->LM'], na.rm = T)

# arrow sizes are these sums*0.2 e.g. 20 shifts per run would be size 4

####################################################################
##### preparing sexual system data ###########
#read in csv
ver.df <- read.csv("Veronica sexual system data.csv")

#tidy up checked gender column
table(ver.df$Gender_literature)
ver.df$Gender_literature[ver.df$Gender_literature==""] <- ver.df$Gender[ver.df$Gender_literature==""] #fill in blanks with McGlone and Richardson data
ver.df$Gender_literature[ver.df$Gender_literature=="H or possibly GD"] <- "H"
ver.df$Gender_literature[ver.df$Gender_literature=="H or GD"] <- "H"
D_or_GD <- which(ver.df$Gender_literature=="D or GD") #species with multiple options
ver.df$Gender_literature[D_or_GD] <- "D" #will treat as D initially

#make binary gender column with gynodioecious and dioecious both as D for dimorphic
ver.df$gender.binary <- gsub("GD", "D", traits$Gender_literature)

saveRDS(ver.df, file="Veronica sexual system data.Rdata")

##########################################################################
##### Section 2 biome dependent and independent biogeographic models #####
##########################################################################

############################################################################

##########################################################################
############## Section Analyses ############
##########################################################################

####################################################################
###### modelling trait changes with biomes in Veronica #############
####################################################################
library(ape)
library(BioGeoBEARS)
library(dplyr)
library(ggplot2) #for heatmap table
library(pGLS)
library(phytools)
library(Rphylopars)
library(plotrix)
library(psychometric)
library(tidyr)

tr <- read.tree("Veronica consensus.tree")
clade <- "Veronica"
traits <- read.csv(file="Trait data from literature.csv")
shifts <- readRDS("Biome shift table expansion shifts Veronica.Rdata")

#combine trait datasets
traits$species <- gsub("Veronica_", "", traits$species) #remove genus name to match tree
traits <- traits[traits$species %in% tr$tip.label,1:12]
#file sources from McGlone and Richardson 2022 data repository https://doi.org/10.7931/qzf7-9q70
traits2 <- read.csv("McGlone_Richardson_2022_data.csv")

traits2 <- traits2[traits2$Genus=="Veronica",]
traits2$SpeciesName <- gsub("Veronica ", "", traits2$SpeciesName)
traits <- left_join(traits, traits2[,c(3,16, 25)], by=c("species"="SpeciesName"))

#fill in missing values
traits.log <- traits
traits.log[,2:14] <-log(traits.log[,2:14]) 

BM2 <- phylopars(trait_data = traits.log, tree = tr, model = "BM") #log transformed version to avoid negative values
traits.complete <- as.data.frame(exp(BM2$anc_recon[1:111,]))
traits.complete$species <- as.character(rownames(traits.complete))
traits.complete <- traits.complete[,c(14,1:13)]

## add functional versions of traits ##
traits.complete$leaf.area.cm2 <- pi*traits.complete$leaf.length.mm*traits.complete$leaf.width.mm*0.01 #assume ellipse shape of leaf
traits.complete$flower.area.mm2 <- traits.complete$corolla.length.mm*traits.complete$corolla.width.mm/2 #assumes diamond or triangle shape of flower
traits.complete$seed.volume.mm3 <-4/3*pi*traits.complete$seed.length.mm*traits.complete$seed.width.mm^2  #assumes ellipsoid
traits.complete$seedSize.mm <- apply(traits.complete[,10:11],1,max)

#### estimate ancestral traits and fill in missing values with closest relatives ####
traits <- traits.complete
traits.log <- traits
traits.log[,2:18] <-log(traits.log[,2:18])  

#pairs(traits[,2:18], pch=19)
#pairs(traits[,c(2:3,6:7,13:14)], pch=19)

BM1 <-phylopars(trait_data = traits, tree = tr, model = "BM", phylo_correlated = F) #remove phylogenetic correlation cos was getting negative values

#log transformed avoids negative values
BM2 <- phylopars(trait_data = traits.log, tree = tr, model = "BM", phylo_correlated = F) #log transformed version to avoid negative values

AIC(BM1, BM2) #BM2 (log transformed) fits best

#back transform values
ancestral.ests <- exp(BM2$anc_recon)
ancestral.ests <- as.data.frame(ancestral.ests)

### add binary version of corolla colour (probability of being white) ####
flower.colour.binary <- gsub("1", "white", traits$corolla.colour)
flower.colour.binary <- gsub("2", "colourful", flower.colour.binary)
flower.colour.binary <- gsub("3", "colourful", flower.colour.binary)
names(flower.colour.binary) <- traits$species

simmap.colour <- make.simmap(tr, flower.colour.binary, pi=c(1,0), nsim=1000)
probability.colour <- describe.simmap(simmap.colour)
#saveRDS(probability.colour, file="flower colour simmap probabilities.Rdata")
#probability.colour <- readRDS("flower colour simmap probabilities.Rdata")
probability.colour.df <- as.data.frame(probability.colour$ace)
probability.colour.df$node <- rownames(probability.colour.df)
colnames(probability.colour.df)[1] <- "colour.prob"
ancestral.ests$node <- rownames(ancestral.ests)

# mean and std error of shifts #
count.col.means <- colSums(probability.colour$count)/1000
count.col.sd <- c(sd(probability.colour$count[,1]),sd(probability.colour$count[,2]),sd(probability.colour$count[,3]))
count.col.se <- count.col.sd/sqrt(1000)

## plot colour simmap ##
cols<-setNames(c("magenta4","white"),c("colourful", "white"))
obj<-densityMap(simmap.colour,states=c("colourful", "white"),plot=FALSE)
obj$cols[1:1001]<-colorRampPalette(c("magenta4","white"), space="Lab")(1001)
plot(obj)

## simmap sexual system ##
gender <- readRDS("Veronica sexual system data.Rdata")

gender <-gender[,c(4,7)] 
gender$gender.binary.numeric <- gsub("H", 0, gender$gender.binary)
gender$gender.binary.numeric <- gsub("D", 1, gender$gender.binary.numeric)
gender$gender.binary.numeric <- as.numeric(gender$gender.binary.numeric)

ver.gender <- gender$gender.binary #make vector of gender
names(ver.gender) <- gender$SpeciesName
cols<-setNames(c("#0072B2", "#F0E442"),levels(as.factor(ver.gender))) #set colours

simmap.gender <- make.simmap(tr, ver.gender, pi=c(0,1), nsim=1000)
gender.summary <- summary(simmap.gender)

probability.gender <- describe.simmap(simmap.gender)
#saveRDS(probability.gender, file="sexual system simmap probabilities.Rdata")
#probability.gender <- readRDS("sexual system simmap probabilities.Rdata")
probability.gender.df <- as.data.frame(probability.gender$ace)
probability.gender.df$node <- rownames(probability.gender.df)
colnames(probability.gender.df)[1] <- "dimorphic.prob"

# mean and std error of shifts #
count.gender.means <- colSums(probability.gender$count)/1000
count.gender.sd <- c(sd(probability.gender$count[,1]),sd(probability.gender$count[,2]),sd(probability.gender$count[,3]))
count.gender.se <- count.gender.sd/sqrt(1000)

##############biomes data#################
#biomes.df <- read.csv("C:/Users/esthe/My Drive/UNIL postdoc/NZ flora/analyses/alt_lat_bands_mid.csv")
biomes.df <- read.csv("alt_lat_bands_mid.csv")

biomes.df$species <- gsub("Veronica_", "", biomes.df$species) #remove genus from species name
biomes.df <- biomes.df[,1:4]#subset to biome data
## add in island species missing from Anne's biome dataset ##
biomes.df[nrow(biomes.df)+1,] <- c("barkeri", 0, 0, 1) # barkeri chathams, coastal [lowland]
biomes.df[nrow(biomes.df)+1,] <- c("benthamii", 0, 1, 1) #benthamii montane, cambell and Aucklands nzpcn.org.nz [lowland and mountain]
biomes.df[nrow(biomes.df)+1,] <- c("breviracemosa", 0, 0, 1)#breviracemosa kermadecs [lowland]
biomes.df[nrow(biomes.df)+1,] <- c("chathamica", 0, 0, 1) # chathamica coastal chatthams [lowland]
biomes.df[nrow(biomes.df)+1,] <- c("dieffenbachii", 0, 0, 1) #dieffenbachii coastal or inland Chathams[lowland]
biomes.df[nrow(biomes.df)+1,] <- c("insularis", 0, 0, 1) #insularis coastal cliffs [lowland]
row.names(biomes.df) <- biomes.df$species #make row names same as tree tip labels
#remove species not in the phylogeny
biomes.df <- biomes.df[biomes.df$species %in% tr$tip.label,]

#add in mountain column (combine subalpine and alpine)
biomes.df$mountain <- paste(biomes.df$alpine, biomes.df$subalpine, sep="")
biomes.df$mountain <- gsub("00", 0, biomes.df$mountain)
biomes.df$mountain[as.numeric(biomes.df$mountain)>0] <-1 
biomes.df$mountain.numeric <- as.numeric(biomes.df$mountain)
biomes.df$mountain.numeric[biomes.df$mountain==1 & biomes.df$lowland==1] <- 0.5
biomes.df$lowland.numeric <- as.numeric(biomes.df$lowland)
biomes.df$lowland.numeric[biomes.df$mountain==1 & biomes.df$lowland==1] <- 0.5
rownames(biomes.df) <- biomes.df$species
biomes.df<- biomes.df[match(tr$tip.label, biomes.df$species),]
biomes.df$biomes.factor <- biomes.df$lowland
biomes.df$biomes.factor <- gsub("1", "L", biomes.df$biomes.factor)
biomes.df$biomes.factor <- paste(biomes.df$biomes.factor, biomes.df$mountain, sep="")
biomes.df$biomes.factor <- gsub("1", "M", biomes.df$biomes.factor)
biomes.df$biomes.factor <- gsub("0", "", biomes.df$biomes.factor)
biomes.df$biomes.factor <- as.factor(biomes.df$biomes.factor)
biomes.df <- left_join(biomes.df, traits, by="species")

#figure 1 density map of flower colour with ring of biome and sexual system and clades
set.width <- rep(0.25, nrow(traits))
names(set.width) <- traits$species
set.width <- set.width[match(tr$tip.label, traits$species)] #order to match phylogeny

biome.cols <- biomes.df$biomes.factor
names(biome.cols) <- gender$SpeciesName 
biome.cols <- gsub("^LM$", "#0072B2", biome.cols)
biome.cols <- gsub("^L$", "#009E73", biome.cols)
biome.cols <- gsub("^M$", "#56B4E9", biome.cols)
names(biome.cols) <- biomes.df$species
biome.cols <- biome.cols[match(tr$tip.label, biomes.df$species)] #reorder to match phylo

gender.cols <- gender$gender.binary.numeric
gender.cols <- gsub("0", "grey", gender.cols)
gender.cols <- gsub("1", "black", gender.cols)
names(gender.cols) <- gender$species
gender.cols <- gender.cols[match(tr$tip.label, gender$SpeciesName)]

clade.names <- c("sun hebes", "snow hebes", "semi-whipcord hebes", "speedwell hebes", "core hebes")
clade.nodes <- c(217,212,210,197,115)

pdf(file="Fig 1 flower colour simmap with clade biome and clade.pdf", height=12, width=12)
plot(obj,lwd=4.5,outline=TRUE,fsize=c(0.7,0.9), type='fan',offset=7, ftype="off",xlim=c(-8,8), outline = T, legend=F)

draw.circle(x=0,y=0,radius=(max(nodeHeights(tr))-6),nv=100, border="black", lty=2, lwd=1) #6mya circle
draw.circle(x=0,y=0,radius=(max(nodeHeights(tr))-4),nv=100, border="black", lty=1, lwd=1) #4mya circle
draw.circle(x=0,y=0,radius=(max(nodeHeights(tr))-2),nv=100, border="black", lty=2, lwd=1) #2mya circle
#draw.circle(x=0,y=0,radius=(max(nodeHeights(tr))-0),nv=100, border="black", lty=2, lwd=1) #0mya circle

points(2.62,2.97, pch=16, cex=2) #point showing lowland transition in core hebes
ring(set.width, tr, col=gender.cols, offset=0.5, style="ring")
ring(set.width, tr, col=biome.cols, offset=0.25, style="ring")

for(i in c(1,2,4:length(clade.names))){
  arc.cladelabels(tree=NULL, text=clade.names[i], node=clade.nodes[i], mark.node=F, ln.offset = 1.15, lab.offset =1.2, lwd=5, cex=1.5)
  
}
arc.cladelabels(tree=NULL, text="semi-whipcord", node=clade.nodes[3], mark.node=F, ln.offset = 1.15, lab.offset =1.25, lwd=5, cex=1.5)
arc.cladelabels(tree=NULL, text="hebes", node=clade.nodes[3], mark.node=F, ln.offset = 1.15, lab.offset =1.2, lwd=5, cex=1.5)

add.color.bar(4,obj$cols,title=substitute(paste(bold('Flower colour'))), subtitle="(probability of white flowers)",lims=c(0,1),
              x=-8,y=-8,prompt=FALSE,fsize=1.25)

add.simmap.legend(leg=c("Hermaphroditic", "Dimorphic"), c("grey", "black"), prompt=F, vertical=TRUE, x=-8, y=-6.75, fsize=1.25)
add.simmap.legend(leg=c("Mountain", "Both", "Lowland"), c("#56B4E9", "#0072B2", "#009E73"), prompt=F, vertical=TRUE, x=-8, y=-5, fsize=1.25)
text(-8.25,-6.3, substitute(paste(bold('Sexual system'))), cex=1.25, pos=4)
text(-8.25,-4.5, substitute(paste(bold('Biome'))), cex=1.25, pos=4)

dev.off()

#################
#save.image(file="Workspace after phylogeny.Rdata")

#load("Workspace after phylogeny.Rdata")

#add colour and gender to ancestral traits df
probability.gender.df <- left_join(probability.gender.df, probability.colour.df, by="node")
ancestral.ests <- left_join(ancestral.ests, probability.gender.df[,c(1,3,4)], by="node")
ancestral.ests <- ancestral.ests[,c(1:17,19:20)] #remove node column and white column

traits.imputed.full <- ancestral.ests #the tip values of all traits i.e. with imputed values
traits.imputed.full$node <- rownames(traits.imputed.full)
traits.imputed <- ancestral.ests[1:111,] #the tip values of all traits i.e. with imputed values

#phylo structure information
prt.ver <- prt(tr)

### ancestral state estimation dataframe ####
vars <- colnames(ancestral.ests)
ancestral.ests$node <- rownames(ancestral.ests)

ancestral.ests <- pivot_longer(ancestral.ests, cols=vars, names_to = c("trait"), values_to = "value")

saveRDS(ancestral.ests, file="ancestral trait estimation Veronica all vars.Rdata")

#### read in biome node data ###
load("BAYAREALIKE+J_M0_v1.Rdata")
node.biomes.df <- as.data.frame(res$ML_marginal_prob_each_state_at_branch_top_AT_node)
node.biomes.df$node <- row.names(node.biomes.df)
node.biomes.df <- node.biomes.df[,c(5,2:4)]
colnames(node.biomes.df)[2:4] <- c("L", "M", "LM")
node.biomes.df$most.likely <- colnames(node.biomes.df[,2:4])[apply(node.biomes.df[,2:4],1,which.max)]

#################################################################
########### trait changes by branch #############################
#################################################################
# plot the phylogeny with  node numbers and axis
pdf(file="Veronica phylogeny with node numbers.pdf", height=30, width=10)
plot(tr, label.offset = 0.3)
nodelabels()
tiplabels()
#axis(1, at=seq(from=0, by=0.25, to=6.25))
axis(1, at=seq(from=max(prt.ver$time_bp), by=-0.5, length=14), labels = seq(from=0, by=0.5, length=14))
mtext("Time (Ma)", side=1, line=3, at=max(prt.ver$time_bp)/2)

abline(v=seq(from=max(prt.ver$time_bp)-1, by=-1, length=6), lty=3, col="black")
abline(v=seq(from=max(prt.ver$time_bp)-0.5, by=-1, length=6), lty=3, col="grey")
dev.off()

#change node names to numbers
for(i in 1:length(tr$tip.label)){
  tip.i <- paste("^",tr$tip.label[i], "$", sep="")
  tip.rows <- grep(tip.i, ancestral.ests$node)
  tip.no <- prt.ver$node[grep(tip.i,prt.ver$label)]
  ancestral.ests$node[tip.rows] <- tip.no
  
  tip.rows.2 <-  grep(tip.i, traits.imputed.full$node)
  traits.imputed.full$node[tip.rows.2] <- tip.no
}

#change to numeric
ancestral.ests <- ancestral.ests[,c(2,1,3)]
ancestral.ests[,2:ncol(ancestral.ests)] <- lapply(ancestral.ests[,2:ncol(ancestral.ests)], function(x) as.numeric(as.character(x)) )

traits.imputed.full$node <- as.numeric(as.character(traits.imputed.full$node))

# set up table to fill with data - each row is a branch section
branches.df<- data.frame("node.older"=numeric(), "node.younger"=numeric(),"age.max"=numeric(), "age.min"=numeric(), "age.median"=numeric(), "biome.older"=character(), "biome.younger"=character())

new.df2 <- as.data.frame(matrix(ncol=7, nrow=length(tr$edge.length))) #blank df with as many rows as branch sections
colnames(new.df2) <- colnames(branches.df)
new.df2[,1:2] <- tr$edge #node number for ends of each branch length
for (j in 1:nrow(new.df2)){
  new.df2$age.max[j] <- prt.ver$time_bp[grep(paste("^", new.df2$node.older[j], "$", sep=""), prt.ver$node)]  #age of oldest node of branch length  
  new.df2$age.min[j] <- prt.ver$time_bp[grep(paste("^", new.df2$node.younger[j], "$", sep=""), prt.ver$node)]  #age of oldest node of branch length  
  
}
new.df2$age.median <- (new.df2$age.max + new.df2$age.min) /2 #age of middle of branch section
branches.df <- rbind(branches.df, new.df2) #add to master df
branches.df <- branches.df[order(branches.df$node.younger),]#order by node number so will match vcv later
rownames(branches.df) <- branches.df$node.younger
# add counts of each shift type for each branch section
shift.types <- data.frame(matrix(ncol=length(table(shifts$shift)), nrow=nrow(branches.df)))
colnames(shift.types) <- names(table(shifts$shift))
branches.df <- cbind(branches.df, shift.types)
shifts$shift.factor <- as.factor(shifts$shift) #make factor version of shifts so missing shift types will be 0s

#get average shift rates for each shift type for each branch
for(k in 1:nrow(branches.df)){
  node.young.k <- branches.df$node.younger[k] # younger node of kth branch length
  branches.df[k,8:11]<- table(shifts$shift.factor[shifts$node==node.young.k])/1000 # table of shift types for kth branch length (shifts labelled as the younger end of branch)
}

branches.df$total.shifts <- rowSums(branches.df[,8:11]) #total shifts for each branch length

#add sums for each biome shifted into
branches.df$'->L' <- branches.df$`M->L,M`+branches.df$`M->LM`
branches.df$'->M' <- branches.df$`L->L,M`+branches.df$`L->LM`
branches.df$'net ->M' <- branches.df$'->M' - branches.df$'->L' 

#data frame for actual trait values for each trait by branch for older adn younger nodes (as opposed to trait changes)
branches.actual.younger.node.df <- branches.df
branches.actual.younger.node.df <- cbind(branches.actual.younger.node.df, data.frame(matrix(ncol = length(vars), nrow = nrow(branches.df))))
colnames(branches.actual.younger.node.df)[16:34] <- vars

#add columns for probability of being in a biome
branches.actual.younger.node.df$L.prob.older <- rep(NA, nrow(branches.actual.younger.node.df))
branches.actual.younger.node.df$M.prob.older <- rep(NA, nrow(branches.actual.younger.node.df))
branches.actual.younger.node.df$LM.prob.older <- rep(NA, nrow(branches.actual.younger.node.df))
branches.actual.younger.node.df$L.prob.younger <- rep(NA, nrow(branches.actual.younger.node.df))
branches.actual.younger.node.df$M.prob.younger <- rep(NA, nrow(branches.actual.younger.node.df))
branches.actual.younger.node.df$LM.prob.younger <- rep(NA, nrow(branches.actual.younger.node.df))

branches.actual.older.node.df <- branches.actual.younger.node.df

# data frame for trait shifts for each branch section for each trait
trait.shifts <- data.frame(matrix(ncol = length(vars), nrow = nrow(branches.df)))
colnames(trait.shifts) <- vars
branches.df <- cbind(branches.df, trait.shifts) #add trait shift cols to branches.df

for(i in 1:nrow(branches.df)){
  node.old.i <- branches.df$node.older[i]
  node.young.i <- branches.df$node.younger[i]
  
  branches.df$biome.older[i] <- node.biomes.df$most.likely[node.biomes.df$node==node.old.i]
  branches.df$biome.younger[i] <- node.biomes.df$most.likely[node.biomes.df$node==node.young.i]
  
  branches.actual.older.node.df$L.prob.older[i] <- node.biomes.df$L[node.biomes.df$node==node.old.i]
  branches.actual.older.node.df$L.prob.younger[i] <- node.biomes.df$L[node.biomes.df$node==node.young.i]
  
  branches.actual.older.node.df$M.prob.older[i] <- node.biomes.df$M[node.biomes.df$node==node.old.i]
  branches.actual.older.node.df$M.prob.younger[i] <- node.biomes.df$M[node.biomes.df$node==node.young.i]
  
  branches.actual.older.node.df$LM.prob.older[i] <- node.biomes.df$LM[node.biomes.df$node==node.old.i]
  branches.actual.older.node.df$LM.prob.younger[i] <- node.biomes.df$LM[node.biomes.df$node==node.young.i]
  
  
  #loop through vars
  for(j in 1:length(vars)){
    var.col <- grep(vars[j], colnames(branches.df)) #column number of jth trait in branches.df
    trait.older.j <- ancestral.ests$value[ancestral.ests$node==node.old.i & ancestral.ests$trait==vars[j]] # value of ancestral estimate for right node, clade and variable combo for older node
    trait.younger.j <- ancestral.ests$value[ancestral.ests$node==node.young.i & ancestral.ests$trait==vars[j]] #same but younger node of branch
    trait.shift.j <- trait.younger.j - trait.older.j
    try(branches.df[i,var.col] <- trait.shift.j)
    branches.actual.younger.node.df[i, var.col] <- trait.younger.j
    branches.actual.older.node.df[i, var.col] <- trait.older.j
    #branches.df[i,var.col] <- trait.shift.j
  }
}

branches.actual.younger.node.df[,35:40] <-branches.actual.older.node.df[,35:40]  #copy over biome probability columns
branches.actual.younger.node.df[,6:7] <- branches.df[,6:7] #copy over biome occupancy columns
branches.actual.older.node.df[,6:7] <- branches.df[,6:7] #copy over biome occupancy columns

saveRDS(branches.df, "Trait and biome shifts by branch section Veronica full traits.Rdata")

##################################
######## modelling ###############
##################################

#variance covariance matrix of nodes and tips
vcv <- vcvPhylo(tr, anc.nodes = T) #generate variance covariance matrix for phylo, including internal nodes

#check skewness of biome shift rates
var.list <- colnames(branches.df)[8:ncol(branches.df)] #list of shift types and traits

pdf("Histograms shift types and traits.pdf", height=10, width=6)
par(mfrow=c(5,2))
for(i in 1:length(var.list)){
  hist(branches.df[,which(colnames(branches.df)==var.list[i])], 20, xlab=var.list[i], main=var.list[i])
}
plot.new()
#dev.off()
for(i in 1:8){
  hist(log(branches.df[,which(colnames(branches.df)==var.list[i])]+1), 20, xlab=paste("Log", var.list[i], sep=" "), main=var.list[i])
}

dev.off()

###### trait change vs biome shift rates ######
##### loop through plotting and fitting models #####
shift.types <- var.list[1:8]
trait.list <- var.list[9:length(var.list)]

shift.names <- c("L to L and M shifts","L to LM shifts", "M to L and M shifts", 
                 "M to LM shifts", "Total biome shifts", "Shifts into lowland", 
                 "Shifts into mountains", "Net shifts into mountains")
trait.names <-c("Change in height (m)", "Change in internode length (mm)", 
                "Change in leaf length (mm)", "Change in leaf width (mm)", 
                "Change in number of flowers", "Change in inflorescence length (cm)", 
                "Change in corolla length (mm)", "Change in corolla width (mm)", 
                "Change in seed length (mm)", "Change in seed width (mm)", 
                "Change in corolla colour", "Change in flower size (mm)", "Change in height (m) McGlone",
                "Change in leaf area (cm2)", "Change in flower area (mm2)", "Change in seed volume (mm3)", 
                "Change in seed size (mm)", "Change in sexual system (probability of being dimorphic)","Change in flower colour (probability of being colourful)") 

res.trait.change.df <- data.frame("shift.type"=character(), "trait"=character(), 
                                  "coefficient.estimate"=numeric(), "coefficient.se"=numeric(),
                                  "R.squared"=numeric(), "R.sq.CI.lower"=numeric(),"R.sq.CI.upper"=numeric(),"p.value"=numeric())

pdf(file="Shifts vs trait changes Veronica.pdf", height=15, width=8)
par(mfrow=c(4,3))
for(i in 1:7){
  data.i <- branches.df[,c(which(colnames(branches.df)==shift.types[i]), 14:ncol(branches.df))]
  colnames(data.i)[1] <- "shift.i"
  for(j in 1:length(trait.list)){
    data.j <- data.i[,c(1, grep(trait.list[j], colnames(data.i)))]
    colnames(data.j)[2] <- "trait.j"  
    plot(data.j$trait.j~jitter(log(data.j$shift.i+1),100), pch=19, xlab= shift.names[i], ylab= trait.names[j])
    abline(h=0, lty=3, col="lightgrey")
    model.j <- pGLS(trait.j~log(shift.i+1), data=data.j, covarmatrix = vcv)
    if(model.j$coefficients[2,4]<0.05){abline(a=model.j$coefficients[1,1], b=model.j$coefficients[2,1])}
    
    R.sq.ci <- CI.Rsq(model.j$`R-Sq`, n=nrow(branches.df), k=1, level=0.95)
    
    new.line <- c(shift.names[i], trait.list[j], model.j$coefficients[2,1:2], model.j$`R-Sq`, R.sq.ci[3:4], model.j$coefficients[2,4])
    names(new.line) <- colnames(res.trait.change.df)
    res.trait.change.df <- rbind(res.trait.change.df, new.line)
    
    
  }#j loop
}#i loop
for(i in 8){ #non-logged version for net shifts
  data.i <- branches.df[,c(which(colnames(branches.df)==shift.types[i]), 14:ncol(branches.df))]
  colnames(data.i)[1] <- "shift.i"
  for(j in 1:length(trait.list)){
    data.j <- data.i[,c(1, grep(trait.list[j], colnames(data.i)))]
    colnames(data.j)[2] <- "trait.j"  
    plot(data.j$trait.j~jitter(data.j$shift.i,100), pch=19, xlab= shift.names[i], ylab= trait.names[j])
    abline(h=0, lty=3, col="lightgrey")
    model.j <- pGLS(trait.j~shift.i, data=data.j, covarmatrix = vcv)
    if(model.j$coefficients[2,4]<0.05){abline(a=model.j$coefficients[1,1], b=model.j$coefficients[2,1])}
    
    R.sq.ci <- CI.Rsq(model.j$`R-Sq`, n=nrow(branches.df), k=1, level=0.95)
    
    new.line <- c(shift.names[i], trait.list[j], model.j$coefficients[2,1:2], model.j$`R-Sq`, R.sq.ci[3:4], model.j$coefficients[2,4])
    names(new.line) <- colnames(res.trait.change.df)
    
    res.trait.change.df <- rbind(res.trait.change.df, new.line)
    
    
  }#j loop
} #i loop

dev.off()

res.trait.change.df$cohens.f2 <- res.trait.change.df$R.squared/(1-res.trait.change.df$R.squared)
res.trait.change.df$cohens.f2.lower <- res.trait.change.df$R.sq.CI.lower/(1-res.trait.change.df$R.sq.CI.lower)
res.trait.change.df$cohens.f2.upper <- res.trait.change.df$R.sq.CI.upper/(1-res.trait.change.df$R.sq.CI.upper)

res.trait.change.df$p.value.adjusted <- p.adjust(res.trait.change.df$p.value, "BH")

res.trait.change.df.trimmed <- res.trait.change.df[res.trait.change.df$shift.type %in% c("Net shifts into mountains", "Shifts into lowland", "Shifts into mountains")& res.trait.change.df$trait %in% c("Height.m", "FlowerSize", "leaf.area.cm2","no.Flowers", "inflorescence.cm","internode.mm", "seedSize.mm"),]
res.trait.change.df.trimmed$p.BH <- p.adjust(res.trait.change.df.trimmed$p.value, "BH")

## heatmap table of results ##
mean.traits <- c(apply(traits[2:18], 2, mean, na.rm = T),mean(probability.gender.df$dimorphic.prob), mean(probability.colour.df$colour.prob)) #get mean value of each trait so can normalise coefficients
names(mean.traits)[c(18,19)] <- c("dimorphic.prob","colour.prob")

included.traits <- trait.list

table.data <- res.trait.change.df[res.trait.change.df$shift.type %in% c("Shifts into lowland", "Shifts into mountains", "Net shifts into mountains")& res.trait.change.df$trait %in% included.traits,c(1:3, 8,10)]
table.data$mean.trait <- rep(mean.traits,3)
table.data$coefficient.scaled <- table.data$coefficient.estimate/table.data$mean.trait

# pdf(file="Model results heatmap cohensf2.pdf", height=10, width=8)
# ggplot(table.data, aes(x = factor(shift.type, levels = c("Shifts into lowland", 
#                                                          "Shifts into mountains", 
#                                                          "Net shifts into mountains")), 
#                        y = factor(trait, levels = rev(c("height.m", "Height.m","internode.mm", 
#                                                         "leaf.length.mm", "leaf.width.mm", "leaf.area.cm2",
#                                                        "no.Flowers","inflorescence.cm", 
#                                                         "corolla.length.mm", "corolla.width.mm", "FlowerSize", "flower.area.mm2",
#                                                         "corolla.colour", "dimorphic.prob", "colour.prob","seed.length.mm", 
#                                                         "seed.width.mm", "seed.volume.mm3", "seedSize.mm"))))) +
#   
#   #make heatmap with geom_tile
#   geom_tile(aes(fill = coefficient.scaled)) + #fill colour divided by mean trait values to scale it
#   
#   #add text of the coefficient estimate to the tile
#   geom_text(aes(label = round(coefficient.estimate, digits = 2)), color = "white") + #values rounded to 2 dp
#   #formatting
#    scale_x_discrete(position = "top", labels = c("into L\n(total)","into M\n(total)", "into M\n(net)")) + #labels at top
#   scale_y_discrete(labels = rev(c("Height (m)","Height (m) McGlone", "Internode length (mm)", "Leaf length (mm)","Leaf width (mm)", "Leaf area (cm2)",
#                                 "Number of flowers","Inflorescence length (cm)","Corolla length (mm)", "Corolla width (mm)", "Flower size (mm)", "Flower area (mm2)",
#                                   "Corolla colour", "Sexual system (probability of being dimorphic)","Corolla colour (probability of being colourful)", "Seed length (mm)", "Seed width (mm)", "Seed volume (mm3)", "Seed size (mm)"))) + #labels
#   scale_fill_gradient2(high = "#D55E00", mid="white", low = "#0072B2") +
#                    theme(legend.position = "bottom", legend.title.align = 0.5, legend.text = element_text(size=10),
#                            panel.grid = element_blank(), 
#                             panel.background = element_rect(fill = "white"),
#                            axis.ticks = element_blank()) +
#   geom_tile(data = table.data[table.data$cohens.f2.lower>0,], fill = NA, color = "black", linewidth = 0.75, height=1, width=1) + #outlines significant tiles in black
#                        labs(y = "Traits", x="Biome shift type", fill="strength of\nassociation")
# 
# dev.off()


table.data2 <- table.data[table.data$trait %in% res.trait.change.df.trimmed$trait,]
table.data2$p.adjusted <-p.adjust(table.data2$p.value, "BH") 

######## heat map binary traits #############
table.binary <- table.data[table.data$trait %in% c("colour.prob", "dimorphic.prob"),]
pdf(file="Model results heatmap reproductive traits.pdf", height=3, width=5)
ggplot(table.binary, aes(x = factor(shift.type, levels = c("Shifts into lowland", 
                                                           "Shifts into mountains", 
                                                           "Net shifts into mountains")), 
                         y = factor(trait, levels = rev(c("dimorphic.prob", "colour.prob"))))) +
  
  #make heatmap with geom_tile
  geom_tile(aes(fill = coefficient.estimate)) + #fill colour divided by mean trait values to scale it
  
  #add text of the coefficient estimate to the tile
  geom_text(aes(label = round(coefficient.estimate, digits = 2)), color = "black") + #values rounded to 2 dp
  #formatting
  scale_x_discrete(position = "top", labels = c("into L\n(total)","into M\n(total)", "into M\n(net)")) + #labels at top
  scale_y_discrete(labels = rev(c("Sexual system\n(probability of being dimorphic)","Corolla colour\n(probability of being pigmented)"))) + #labels
  scale_fill_gradient2(high = "#D55E00", mid="white", low = "#0072B2") +
  theme(legend.position = "bottom", legend.title.align = 0.5, legend.text = element_text(size=10),
        panel.grid = element_blank(), 
        panel.background = element_rect(fill = "white"),
        axis.ticks = element_blank()) +
  geom_tile(data = table.binary[table.binary$p.value<0.05,], fill = NA, color = "black", linewidth = 0.75, height=1, width=1) + #outlines significant tiles in black
  labs(y = "", x="Biome shift type", fill="strength of\nassociation")

dev.off()

pdf(file="Model results heatmap continuous traits BH p adjustment.pdf", height=6, width=6)
ggplot(table.data2, aes(x = factor(shift.type, levels = c("Shifts into lowland", 
                                                          "Shifts into mountains", 
                                                          "Net shifts into mountains")), 
                        y = factor(trait, levels = rev(c("Height.m", "internode.mm",
                                                         "leaf.area.cm2", "no.Flowers","inflorescence.cm", 
                                                         "FlowerSize", "seedSize.mm"))))) +
  
  #make heatmap with geom_tile
  geom_tile(aes(fill = coefficient.scaled)) + #fill colour divided by mean trait values to scale it
  
  #add text of the coefficient estimate to the tile
  geom_text(aes(label = round(coefficient.estimate, digits = 2)), color = "white") + #values rounded to 2 dp
  #formatting
  scale_x_discrete(position = "top", labels = c("into L\n(total)","into M\n(total)", "into M\n(net)")) + #labels at top
  scale_y_discrete(labels = rev(c("Height (m)", "Internode length (mm)", "Leaf area (cm2)",
                                  "Number of flowers","Inflorescence length (cm)", "Flower size (mm)",
                                  "Seed size (mm)"))) + #labels
  scale_fill_gradient2(high = "#D55E00", mid="white", low = "#0072B2") +
  theme(legend.position = "bottom", legend.title.align = 0.5, legend.text = element_text(size=10),
        panel.grid = element_blank(), 
        panel.background = element_rect(fill = "white"),
        axis.ticks = element_blank()) +
  geom_tile(data = table.data2[table.data2$p.adjusted<0.05,], fill = NA, color = "black", linewidth = 0.75, height=1, width=1) + #outlines significant tiles in black
  labs(y = "Traits", x="Biome shift type", fill="strength of\nassociation")

dev.off()

pdf(file="Model results heatmap continuous traits no p adjustment.pdf", height=6, width=6)
ggplot(table.data2, aes(x = factor(shift.type, levels = c("Shifts into lowland", 
                                                          "Shifts into mountains", 
                                                          "Net shifts into mountains")), 
                        y = factor(trait, levels = rev(c("Height.m", "internode.mm",
                                                         "leaf.area.cm2", "no.Flowers","inflorescence.cm", 
                                                         "FlowerSize", "seedSize.mm"))))) +
  
  #make heatmap with geom_tile
  geom_tile(aes(fill = coefficient.scaled)) + #fill colour divided by mean trait values to scale it
  
  #add text of the coefficient estimate to the tile
  geom_text(aes(label = round(coefficient.estimate, digits = 2)), color = "white") + #values rounded to 2 dp
  #formatting
  scale_x_discrete(position = "top", labels = c("into L\n(total)","into M\n(total)", "into M\n(net)")) + #labels at top
  scale_y_discrete(labels = rev(c("Height (m)", "Internode length (mm)", "Leaf area (cm2)",
                                  "Number of flowers","Inflorescence length (cm)", "Flower size (mm)",
                                  "Seed size (mm)"))) + #labels
  scale_fill_gradient2(high = "#D55E00", mid="white", low = "#0072B2") +
  theme(legend.position = "bottom", legend.title.align = 0.5, legend.text = element_text(size=10),
        panel.grid = element_blank(), 
        panel.background = element_rect(fill = "white"),
        axis.ticks = element_blank()) +
  geom_tile(data = table.data2[table.data2$p.value<0.05,], fill = NA, color = "black", linewidth = 0.75, height=1, width=1) + #outlines significant tiles in black
  labs(y = "Traits", x="Biome shift type", fill="strength of\nassociation")

dev.off()

# ########### trait value influence on biome shifts #################
# ### trait value and end of branch vs. probability of being in the biome at the start
# res.trait.value.younger.node.df <- data.frame("biome.prob.type"=character(), "trait"=character(), 
#                          "coefficient.estimate"=numeric(), "coefficient.se"=numeric(),
#                          "R.squared"=numeric(), "R.sq.CI.lower"=numeric(),"R.sq.CI.upper"=numeric(),"p.value"=numeric())
# trait.names2 <- gsub("Change in ", "", trait.names)
# 
# biome.prob.young <- colnames(branches.actual.younger.node.df)[38:40]
# biome.prob.old <- colnames(branches.actual.younger.node.df)[35:37]
# biome.prob.names <- c("Probability of being L", "Probability of being M", "Probability of being LM")
# 
# pdf(file="Biome probabilities older vs trait actual values younger Veronica lm.pdf", height=15, width=8)
# par(mfrow=c(4,3))
# for(i in 1:length(biome.prob.old)){
#   data.i <- branches.actual.younger.node.df[,c(which(colnames(branches.actual.younger.node.df)==biome.prob.old[i]), 16:34)]
# 
#   colnames(data.i)[1] <- "biome.prob.i"
#   for(j in 1:length(trait.list)){
#     data.j <- data.i[,c(1, grep(trait.list[j], colnames(data.i)))]
#     colnames(data.j)[2] <- "trait.j"  
#     data.j$trait.j <- data.j$trait.j/max(data.j$trait.j) #scale trait so are comparable
#     if(trait.list[j] %in% c("seed.Size.mm", "colour.prob", "dimorphic.prob")){ #no-logged version for non skewed vars
#       plot(jitter(data.j$trait.j,100)~jitter(data.j$biome.prob.i, 100), pch=19, xlab= biome.prob.names[i], ylab= trait.names2[j])
#       model.j <- lm(trait.j~biome.prob.i, data=data.j)
#       
#     }else{ #logged version
#       plot(jitter(log(data.j$trait.j),100)~jitter(data.j$biome.prob.i, 100), pch=19, xlab= biome.prob.names[i], ylab= trait.names2[j])
#       model.j <- lm(log(trait.j)~biome.prob.i, data=data.j)
#       
#     }
#    summary.j <- summary(model.j)
#     if(summary.j$coefficients[2,4]<0.05){abline(a=summary.j$coefficients[1,1], b=summary.j$coefficients[2,1])}
#     
#     R.sq.ci <- CI.Rsq(summary.j$r.squared, n=nrow(branches.actual.younger.node.df), k=1, level=0.95)
#     
#     new.line <- c(biome.prob.old[i], trait.list[j], summary.j$coefficients[2,1:2], summary.j$r.squared, R.sq.ci[3:4], summary.j$coefficients[2,4])
#     names(new.line) <- colnames(res.trait.value.younger.node.df)
#     res.trait.value.younger.node.df <- rbind(res.trait.value.younger.node.df, new.line)
#     
#     
#   }#j loop
# }#i loop
# 
# dev.off()
# 
# res.trait.value.younger.node.df$cohens.f2 <- res.trait.value.younger.node.df$R.squared/(1-res.trait.value.younger.node.df$R.squared)
# res.trait.value.younger.node.df$cohens.f2.lower <- res.trait.value.younger.node.df$R.sq.CI.lower/(1-res.trait.value.younger.node.df$R.sq.CI.lower)
# res.trait.value.younger.node.df$cohens.f2.upper <- res.trait.value.younger.node.df$R.sq.CI.upper/(1-res.trait.value.younger.node.df$R.sq.CI.upper)
# 
# res.trait.value.younger.node.df$p.value.adjusted <- p.adjust(res.trait.value.younger.node.df$p.value, "BH")
# 
# res.trait.value.younger.node.df$mean.trait <-  rep(mean.traits,3)
# 
# res.trait.value.younger.node.df.trimmed <- res.trait.value.younger.node.df[res.trait.value.younger.node.df$trait %in% c("Height.m", "FlowerSize", "leaf.area.cm2","no.Flowers", "inflorescence.cm","internode.mm", "seedSize.mm", "colour.prob"),]
# res.trait.value.younger.node.df.trimmed$p.value.adjusted<- p.adjust(res.trait.value.younger.node.df.trimmed$p.value, "BH")
# 
# pdf(file="Model results older biome probabilites and actual values younger traits heatmap.pdf", height=6, width=6)
# ggplot(res.trait.value.younger.node.df.trimmed, aes(x = factor(biome.prob.type, levels = c("L.prob.older", "LM.prob.older",  "M.prob.older")), 
#                                             y = factor(trait, levels = rev(c("Height.m","internode.mm", 
#                                                                              "leaf.area.cm2",
#                                                                              "no.Flowers","inflorescence.cm", 
#                                                                               "FlowerSize", "dimorphic.prob",
#                                                                              "colour.prob","seedSize.mm" ))))) +
#   
#   #make heatmap with geom_tile
#   geom_tile(aes(fill = coefficient.estimate)) + #fill colour divided by mean trait values to scale it
#   
#   #add text of the coefficient estimate to the tile
#   geom_text(aes(label = round(coefficient.estimate, digits = 2)), color = "white") + #values rounded to 2 dp
#   #formatting
#   scale_x_discrete(position = "top", labels = (c("L","LM", "M"))) + #labels at top
# 
#   scale_y_discrete(labels = rev(c("Height (m)","Internode length (mm)", "Leaf area (cm2)",
#                                   "Number of flowers","Inflorescence length (cm)","Flower size (mm)",
#                                   "Sexual system (probability of being dimorphic)",
#                                   "Corolla colour (probability of being colourful)","Seed size (mm)"))) + #labels
#    scale_fill_gradient2(high = "#D55E00", mid="white", low = "#0072B2", guide = "colourbar") +
#   theme(legend.position = "bottom", legend.title.align = 0.5, legend.text = element_text(size=8),
#         panel.grid = element_blank(), 
#         panel.background = element_rect(fill = "white"),
#         axis.ticks = element_blank()) +
#   geom_tile(data = res.trait.value.younger.node.df.trimmed[res.trait.value.younger.node.df.trimmed$p.value.adjusted<0.05,], fill = NA, color = "black", linewidth = 0.75, height=0.96, width=0.98) + #outlines significant tiles in black
#   labs(y = "Traits at end of branch", x="Biome (probability at start of branch)", fill="strength of\nassociation")
# 
# dev.off()
# 
# ### trait value at start of branch vs. probability of being in the biome at the end
# res.trait.value.older.node.df <- data.frame("biome.prob.type"=character(), "trait"=character(), 
#                                               "coefficient.estimate"=numeric(), "coefficient.se"=numeric(),
#                                             "p.value"=numeric())
#                                             
# 
# pdf(file="Biome probabilities younger vs trait actual values older Veronica glm.pdf", height=15, width=8)
# par(mfrow=c(4,3))
# for(i in 1:length(biome.prob.young)){
#   data.i <- branches.actual.older.node.df[,c(which(colnames(branches.actual.older.node.df)==biome.prob.young[i]), 16:34)]
#     colnames(data.i)[1] <- "biome.prob.i"
#     data.i$biome.prob.i.binary <- round(data.i$biome.prob.i, digits=0)
#   
#   for(j in 1:length(trait.list)){
#     data.j <- data.i[,c(1, 21,grep(trait.list[j], colnames(data.i)))]
#     colnames(data.j)[3] <- "trait.j"  
#     data.j$trait.j <- data.j$trait.j/max(data.j$trait.j) #scale trait
#     # plot(jitter(data.j$biome.prob.i, 200)~jitter(data.j$trait.j,200), pch=19, ylab= biome.prob.names[i], xlab= trait.names2[j], xlim=c(0,max(data.j$trait.j)))
#     # model.j <- pGLS(biome.prob.i~trait.j, data=data.j, covarmatrix = vcv)
#     # abline(model.j)
#     # if(model.j$coefficients[2,4]<0.05){abline(a=model.j$coefficients[1,1], b=model.j$coefficients[2,1], col="red", lwd=2)}
#     # 
#     if(trait.list[j] %in% c("seed.Size.mm", "colour.prob", "dimorphic.prob")){ #no-logged version for non skewed vars
#       plot(jitter(data.j$biome.prob.i.binary,0.25)~jitter(data.j$trait.j, 100), pch=19, ylab= biome.prob.names[i], xlab= trait.names2[j])
#       model.j <- glm(biome.prob.i.binary~trait.j, data=data.j, family="binomial")
#       newdata <- data.frame("trait.j" =seq(min(data.j$trait.j), max(data.j$trait.j), length=20))
#       
#       predict.j <- predict.glm(model.j, newdata, type = "response")
#       
#     }else{ #logged version
#       plot(jitter(data.j$biome.prob.i.binary,0.25)~jitter(log(data.j$trait.j),100), pch=19, ylab= biome.prob.names[i], xlab= paste("log", trait.names2[j], sep=" "))
#       data.j$trait.j.log <- log(data.j$trait.j)
#       model.j <- glm(biome.prob.i.binary~trait.j.log, data=data.j, family="binomial")
#       newdata <- data.frame("trait.j.log" =seq(min(data.j$trait.j.log), max(data.j$trait.j.log), length=20))
#       
#      predict.j <- predict.glm(model.j, newdata, type = "response")
#      
#     }
#     summary.j <- summary(model.j)
#     if(summary.j$coefficients[2,4]<0.05){lines(predict.j~newdata[,1])}
#     
#     
#   #    plot(jitter(data.j$trait.j,100)~jitter(data.j$biome.prob.i, 100), pch=19, xlab= biome.prob.names[i], ylab= trait.names2[j])
#   # model.j <- pGLS(trait.j~biome.prob.i, data=data.j, covarmatrix = vcv)
#   #if(model.j$coefficients[2,4]<0.05){abline(a=model.j$coefficients[1,1], b=model.j$coefficients[2,1])}
# 
#       
#   #  R.sq.ci <- CI.Rsq(summary.j$r.squared, n=nrow(branches.actual.older.node.df), k=1, level=0.95)
#     
#     new.line <- c(biome.prob.young[i], trait.list[j], summary.j$coefficients[2,1:2], summary.j$coefficients[2,4])
#     names(new.line) <- colnames(res.trait.value.older.node.df)
#     res.trait.value.older.node.df[nrow(res.trait.value.older.node.df)+1,] <- new.line
#     
#     
#   }#j loop
# }#i loop
# 
# dev.off()
# 
# # res.trait.value.older.node.df$cohens.f2 <- res.trait.value.older.node.df$R.squared/(1-res.trait.value.older.node.df$R.squared)
# # res.trait.value.older.node.df$cohens.f2.lower <- res.trait.value.older.node.df$R.sq.CI.lower/(1-res.trait.value.older.node.df$R.sq.CI.lower)
# # res.trait.value.older.node.df$cohens.f2.upper <- res.trait.value.older.node.df$R.sq.CI.upper/(1-res.trait.value.older.node.df$R.sq.CI.upper)
# 
# res.trait.value.older.node.df$p.value.adjusted <- p.adjust(res.trait.value.older.node.df$p.value, "BH")
# 
# res.trait.value.older.node.df$mean.trait <-  rep(mean.traits,3)
# 
# res.trait.value.older.node.df.trimmed <- res.trait.value.older.node.df[res.trait.value.older.node.df$trait %in% c("Height.m", "FlowerSize", "leaf.area.cm2","no.Flowers", "inflorescence.cm","internode.mm", "seedSize.mm", "colour.prob", "dimoprhic.prob"),]
# res.trait.value.older.node.df.trimmed$p.value.adjusted<- p.adjust(res.trait.value.older.node.df.trimmed$p.value, "BH")
# 
# res.trait.value.older.node.df.trimmed[,3:7] <- sapply(res.trait.value.older.node.df.trimmed[,3:7],as.numeric)
# 
# pdf(file="Model results younger biome probabilites and actual values older traits heatmap.pdf", height=6, width=6)
# ggplot(res.trait.value.older.node.df.trimmed, aes(x = factor(biome.prob.type, levels = c("L.prob.younger", "LM.prob.younger",  "M.prob.younger")), 
#                                                     y = factor(trait, levels = rev(c("Height.m","internode.mm", 
#                                                                                      "leaf.area.cm2",
#                                                                                      "no.Flowers","inflorescence.cm", 
#                                                                                      "FlowerSize", "dimorphic.prob",
#                                                                                      "colour.prob","seedSize.mm" ))))) +
#   
#   #make heatmap with geom_tile
#   geom_tile(aes(fill = coefficient.estimate)) + #fill colour divided by mean trait values to scale it
#   
#   #add text of the coefficient estimate to the tile
#   geom_text(aes(label = round(coefficient.estimate, digits = 2)), color = "white") + #values rounded to 2 dp
#   #formatting
#    scale_x_discrete(position = "top", labels = (c("L","LM", "M"))) + #labels at top
# 
#   scale_y_discrete(labels = rev(c("Height (m)","Internode length (mm)", "Leaf area (cm2)",
#                                   "Number of flowers","Inflorescence length (cm)",
#                                   "Flower size (mm)","Sexual system (probability of being dimorphic)",
#                                   "Corolla colour (probability of being colourful)", "Seed size (mm)"))) + #labels
#   scale_fill_gradient2(high = "#D55E00", mid="white", low = "#0072B2") +
#   theme(legend.position = "bottom", legend.title.align = 0.5, legend.text = element_text(size=10),
#         panel.grid = element_blank(), 
#         panel.background = element_rect(fill = "white"),
#         axis.ticks = element_blank()) +
#   geom_tile(data = res.trait.value.older.node.df.trimmed[res.trait.value.older.node.df.trimmed$p.value.adjusted<0.05,], fill = NA, color = "black", linewidth = 0.75, height=0.96, width=0.98) + #outlines significant tiles in black
#   labs(y = "Traits at start of branch", x="Biome (probability at end of branch)", fill="strength of\nassociation")
# 
# dev.off()

############## trait values by biome probability at nodes ##############
node.biomes.df$node <- as.numeric(as.character(node.biomes.df$node))
traits.imputed.full <- left_join(traits.imputed.full, node.biomes.df, by="node")

res.trait.value.biome.df <- data.frame("biome"=character(), "trait"=character(), 
                                       "coefficient.estimate"=numeric(), "coefficient.se"=numeric(),
                                       "R.squared"=numeric(), "R.sq.CI.lower"=numeric(),"R.sq.CI.upper"=numeric(),"p.value"=numeric())



biomes <- c("L", "LM", "M")


### trait values by biome ###
#complete trait values (with imputed)
trait.names2 <- gsub("Change in ", "", trait.names)
#y.lims <- list(c(0,15),c(0,50), c(0,120), c(0,50), c(0,))
trait.aov.results <- data.frame("trait" =character(), "F value"=numeric(), "p-value"=numeric(), "contrast.L.LM"=numeric(), "contrast.L.M"=numeric(), "contrast.LM.M"=numeric())
trait.aov.results.log <- trait.aov.results

contrast.letters.full <- list(c("a", "a", "a"), c("a", "a", "b"), c("a", "a", "b"), c("a", "b", "c"), c("a", "a", "b"), c("a", "a", "b"), 
                              c("a", "a", "a"), c("a", "a", "a"), c("a", "a", "a"), c("a", "a", "a"), 
                              c("a", "a", "b"), c("a", "a", "a"), c("a", "a", "a"), c("a", "b", "c"), c("a", "a", "a"), c("a", "a", "a") )
contrast.letters <- list(c("a", "a", "a"), c("a", "a", "b"), c("a", "b", "c"), c("a", "a", "b"), c("a", "a", "b"), 
                         c("a", "a", "a"), c("a", "a", "a") )


trait.list.selected <- trait.list[c(13,2,14,5:6,12,17)]
trait.names.selected <- trait.names2[c(1,2,14,5:6,12,17)]


pdf(file="Traits by biomes boxplots selected.pdf", height=8, width=6)
par(mfrow=c(3,3), mar=c(4.5,4.5,0.5,0.5))

for(i in 1:length(trait.list.selected)){
  data.i <- biomes.df[,c(1:8, grep(trait.list.selected[i], colnames(biomes.df)))]
  # rownames(data.i) <- biomes.df$species
  
  colnames(data.i)[9] <- "trait.i"
  aov.1 <- phylANOVA(tr, data.i$biomes.factor, data.i$trait.i, posthoc = T, p.adj = "BH")
  aov.res.i <- c(trait.list.selected[i], aov.1$F, aov.1$Pf, aov.1$Pt[2,1], aov.1$Pt[3,1],aov.1$Pt[3,2])
  trait.aov.results[i,] <- aov.res.i
  
  aov.2 <- phylANOVA(tr, data.i$biomes.factor, log(data.i$trait.i), posthoc = T, p.adj = "BY")
  trait.aov.results.log[i,] <- c(trait.list.selected[i], aov.2$F, aov.2$Pf, aov.2$Pt[2,1], aov.2$Pt[3,1],aov.2$Pt[3,2])
  
  plot1 <- boxplot((data.i$trait.i)~data.i$biomes.factor, col=c("#009E73","#0072B2", "#56B4E9"), ylab = trait.names.selected[i], xlab="biome", log="y")
  if(aov.2$Pf<0.05){text(1:3, plot1$stats[3,], unlist(contrast.letters[i]), col="white")}
  #rownames(biomes.df) <- biomes.df$species
  text(0.6,max(data.i$trait.i), paste(letters[i], ")", sep=""))
}
dev.off()


#### threshold model biomes and binary traits ####
library(coda)

thresh.df <- biomes.df[,c(1,4:5, 19)]
thresh.df$corolla.binary <- thresh.df$corolla.colour #make a binary version of flower colour
thresh.df$corolla.binary <- gsub("1", "white", thresh.df$corolla.colour)
thresh.df$corolla.binary[!thresh.df$corolla.binary=="white"] <- "colourful"
thresh.df <- left_join(thresh.df, gender[,c(1,4)], by=c("species"="SpeciesName"))
rownames(thresh.df) <- thresh.df$species
thresh.df <- thresh.df[,c(2:3,5:6)]
thresh.df$gender.binary.numeric <- as.factor(thresh.df$gender.binary.numeric)
thresh.df$corolla.binary <- as.factor(thresh.df$corolla.binary)
thresh.df$lowland <- as.factor(thresh.df$lowland)
thresh.df$mountain <- as.factor(thresh.df$mountain)

thresh.m.gender <- thresh.df[,c(2:3)]
thresh.l.gender <- thresh.df[,c(1,3)]
thresh.m.colour <- thresh.df[,c(2,4)]
thresh.l.colour <- thresh.df[,c(1,4)]
thresh.list <- list(thresh.l.colour, thresh.l.gender, thresh.m.colour, thresh.m.gender)
ngen=10e6

thresh.res.df <- data.frame("variable1"=character(), "variable2"=character(), 
                            "mean.r"=numeric(), "mean.r.no.burnin"=numeric(),
                            "hpd.lower"=numeric(), "hpd.upper"=numeric())

for(i in 1:length(thresh.list)){
  df.i <- thresh.list[[i]]
  mcmc.i <- threshBayes(tr, df.i, type=c("disc", "disc"), ngen=ngen, plot=F, control = list(sample = 1000))
  pdf(file=paste(colnames(df.i)[1],colnames(df.i)[2],"mcmc_plots.pdf", sep="_"))
  plot(density(mcmc.i))
  plot(mcmc.i)
  dev.off()
  print(mcmc.i)
  
  r.mcmc <- tail(mcmc.i$par$r, 0.8*nrow(mcmc.i$par))
  r.i <- mean(mcmc.i$par$r)
  r.i.no.burnin <- mean(r.mcmc)
  class(r.mcmc) <- "mcmc"
  hpd.i <- HPDinterval(r.mcmc)
  
  thresh.res.df[i,] <- c(colnames(df.i), r.i,r.i.no.burnin, hpd.i[1,])
  
  
}

saveRDS(thresh.res.df, file="Threshbayes results2.Rdata")


#### phylogenetic glm ###
library(phylolm)
thresh.df$colour.binary <- gsub("colourful", 1,thresh.df$corolla.binary)
thresh.df$colour.binary <- as.factor(gsub("white", 0,thresh.df$colour.binary))

thresh.df$dimorphic.binary <- gsub("D", 1, thresh.df$gender.binary)
thresh.df$dimorphic.binary <- as.factor(gsub("H", 0, thresh.df$dimorphic.binary))


glm1 <- phyloglm(dimorphic.binary~lowland, data=thresh.df, phy=tr, method = "logistic_MPLE", btol=10, boot = 100)
glm2 <- phyloglm(dimorphic.binary~mountain, data=thresh.df, phy=tr, method = "logistic_MPLE", btol=22, boot = 1000)
glm2.5 <- phyloglm(dimorphic.binary~lowland+mountain, data=thresh.df, phy=tr, method = "logistic_MPLE", btol=20, boot = 100)

glm3 <- phyloglm(colour.binary~lowland, data=thresh.df, phy=tr, method = "logistic_MPLE", btol=10, boot = 100) 
glm4 <- phyloglm(colour.binary~mountain, data=thresh.df, phy=tr, method = "logistic_MPLE", btol=10, boot = 100) 
glm4.5 <- phyloglm(colour.binary~mountain+lowland, data=thresh.df, phy=tr, method = "logistic_MPLE", btol=10, boot = 100) 

summary(glm2.5)

# biome as factor with three states
thresh.df$biome.factor <- biomes.df$biomes.factor #rows already in right order

glm5 <- phyloglm(dimorphic.binary~biome.factor, data=thresh.df, phy=tr, method = "logistic_MPLE", btol=20, boot = 1000)
summary(glm5)
glm6 <- phyloglm(colour.binary~biome.factor, data=thresh.df, phy=tr, method = "logistic_MPLE", btol=10, boot = 1000)

saveRDS(glm5, file="pglm sexual system and biomes.Rdata")
saveRDS(glm6, file="pglm flower colour and biomes.Rdata")

