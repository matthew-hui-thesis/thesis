required.packages <- c("data.table", "magrittr", "readr")
packages2beInstalled <- required.packages[!(required.packages %in% unname(installed.packages()[,1]))]
if(length(packages2beInstalled)) install.packages(packages2beInstalled)

library(data.table)
library(magrittr)
library(readr)

options(nwarnings = 10000)
if (Sys.info()["sysname"] == "Windows") {tryCatch({memory.limit(size = 16000)})}

set.seed(runif(1, max = .Machine$integer.max) %>% ceiling)

displayStr <- "   "

total_start_Time <- Sys.time()

#+++++++++++++++++++++++++++#
#### +++++Function+++++ #####
#+++++++++++++++++++++++++++#

diffTime_to_display <- function(dTime) {
  Ndays = as.numeric(dTime, units="days") %>% floor
  Nhours = as.numeric(dTime, units="hours") %>% floor
  Nmins = as.numeric(dTime, units="mins") %>% floor
  Nsecs = as.numeric(dTime, units="secs") %>% floor
  displayStr = paste0(Ndays, " days ", Nhours - Ndays * 24, " hours ", Nmins - Nhours * 60, " mins ", Nsecs - Nmins * 60, " secs")
  return(displayStr)
}

list_for_dput <- function(arrVar) {
  listVar = lapply(arrVar, function(cVar) {get(cVar)})
  names(listVar) = arrVar
  return(listVar)
}

dgetList_to_var <- function(lstVar) {invisible(
  for (i in seq(length(lstVar))) {assign(names(lstVar)[i], lstVar[[i]], envir = .GlobalEnv)}
)}

cat("\n\n")
cat("*****************************\n")
cat("******** Preparation ********\n")
cat("*****************************\n")

#+++++++++++++++++++++++++++#
#### +++++Settings+++++ #####
#+++++++++++++++++++++++++++#

baseWd <- "C:/Users/Matthew/Documents/NONMEM database/Mixture Model Evaluation Database/2M59049/D04"
projName <- "2M59049_D04"

# Skip run part
skip_generating_data_temp <- F
skip_generating_sim_stream <- F
skip_generating_fit_mix_stream <- F
skip_generating_fit_nomix_stream <- F
skip_generating_fit_mixTV_stream <- F
skip_sim <- F
skip_fit_mix <- F
skip_fit_nomix <- F
skip_fit_mixTV <- F

# Test run indicator (specify the range of dataset# to try) (NULL to perform full)
testRun <- NULL
if (!is.null(testRun)) {
  cat("\n!!!!---- Non-full run has been specified ----!!!!\n")
}

cat("\nReading settings...")

#### Check for previous setting ####
if ((regexpr("SETTING.R", dir()) %>% equals(-1) %>% not %>% any %>% not) ||
    !dir.exists(baseWd)) {
  #### Create directory if not exist yet ####
  dir.create(baseWd, showWarnings = F, recursive = T)
  setwd(baseWd)
  
  #### Dose and directory Settings ####
  dose <- 500
  subDirSize <- 100 #Number of dataset to contain within a subdirectory
  
  #### PsN Settings ####
  
  #Simulation PsN options
  sim_nmfe_options <- "-prdefault"
  sim_clean <- 3
  sim_nm_output <- ""
  sim_threads <- if (Sys.info()["sysname"] == "Windows") {4} else {16}
  
  #Fit mixture PsN options
  fitMix_nmfe_options <- "-prdefault"
  fitMix_clean <- 3
  fitMix_nm_output <- "coi,cor,cov,ext,grd,phi,phm,shk"
  fitMix_threads <- if (Sys.info()["sysname"] == "Windows") {4} else {16}
  fitMix_min_retries <- 0
  fitMix_nodes <- ""
  fitMix_picky <- F
  fitMix_retries <- 4
  fitMix_tweak_inits <- T
  
  #Fit no mixture PsN options
  fitNoMix_nmfe_options <- "-prdefault"
  fitNoMix_clean <- 3
  fitNoMix_nm_output <- "coi,cor,cov,ext,grd,phi,phm,shk"
  fitNoMix_threads <- if (Sys.info()["sysname"] == "Windows") {4} else {16}
  fitNoMix_min_retries <- 0
  fitNoMix_nodes <- ""
  fitNoMix_picky <- F
  fitNoMix_retries <- 4
  fitNoMix_tweak_inits <- T
  
  #Fit mixture (true value) PsN options
  fitMixTV_nmfe_options <- "-prdefault"
  fitMixTV_clean <- 3
  fitMixTV_nm_output <- "coi,cor,cov,ext,grd,phi,phm,shk"
  fitMixTV_threads <- if (Sys.info()["sysname"] == "Windows") {4} else {16}
  fitMixTV_min_retries <- 0
  fitMixTV_nodes <- ""
  
  #Initial estimates deviations (in fraction) from true values
  IE_fix_minDev <- 0.05 # Both direction
  IE_fix_maxDev <- 0.2 # Both direction
  IE_LOGIT_MIXP_minDev <- 0.1 # in absolute value
  IE_LOGIT_MIXP_maxDev <- 1 # in absolute value
  IE_rnd_minDev <- 0.05 # Always larger
  IE_rnd_maxDev <- 0.5 # Always larger
  
  #### Sampling time setting ####
  # Timings to take samples depend on the pharmacokinetic profiles
  # Provide array of factors of tmax and t12 (half-life) to use in determining the sampling times
  # Two arrays must be of same length
  
  tmaxFactor <- c(0.25, 0.5, 0.75, 1,    1,   1,    1, 1,   1, 1, 1, 1)
  t12Factor <-  c(   0,   0,    0, 0, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4)
  
  #### Defining parameter space ####
  
  paramNames <- c("TVCL", "RCL", "CVCL", "TVVD", "CVVD", "TVKA", "CVKA", "MIXP", "W_PROP", "SampleSize")
  paramRngs <- list(
    list(c(1, 3.7), c(3.7, 6.3), c(6.3, 9)),
    list(c(1, 2), c(2, 3), c(3, 4)),
    list(c(0.1, 0.23), c(0.23, 0.37), c(0.37, 0.5)),
    list(c(50, 367), c(367, 683), c(683, 1000)),
    list(c(0.1, 0.23), c(0.23, 0.37), c(0.37, 0.5)),
    list(c(0.1, 0.37), c(0.37, 0.63), c(0.63, 0.9)),
    list(c(0.1, 0.23), c(0.23, 0.37), c(0.37, 0.5)),
    list(c(0.05, 0.35), c(0.35, 0.65), c(0.65, 0.95)),
    list(c(0.1, 0.17), c(0.17, 0.23), c(0.23, 0.3)),
    list(c(20, 50), c(50, 100), c(100, 300))
  ) %>% set_names(paste0("rng_", paramNames))
  
  #### Number of dataset to simulate per scenario
  
  NdatPerScn <- 1
  
  #### Defining control stream format ####
  
  simStrTemp <- paste(sep = "\n",
                      "$PROBLEM REPLACE_PROB_NAME",
                      "",
                      "$INPUT ID AMT TIME DV",
                      "",
                      "$DATA REPLACE_DATA_TEMPLATE_NAME IGNORE=@",
                      "",
                      "$SIM (REPLACE_SEED NEW) ONLYSIM",
                      "",
                      "$SUBROUTINE ADVAN2 TRANS2",
                      "",
                      "$PK",
                      ";THETA",
                      "\tTVCL=THETA(1)",
                      "\tRCL=THETA(2)",
                      "\tCVCL=THETA(3)",
                      "\tTVVD=THETA(4)",
                      "\tCVVD=THETA(5)",
                      "\tTVKA=THETA(6)",
                      "\tCVKA=THETA(7)",
                      "\t;MIXP: SEE $MIX BLOCK",
                      "\tW_PROP=THETA(9)",
                      ";ETA",
                      "\tETA_CL=ETA(1)",
                      "\tETA_VD=ETA(2)",
                      "\tETA_KA=ETA(3)",
                      ";TRANSFORMATION",
                      "\tSD_ETA_CL=SQRT(LOG(CVCL**2+1))",
                      "\tSD_ETA_VD=SQRT(LOG(CVVD**2+1))",
                      "\tSD_ETA_KA=SQRT(LOG(CVKA**2+1))",
                      ";CLEARANCE",
                      "\tIF (MIXNUM.EQ.1) POPCL=TVCL",
                      "\tIF (MIXNUM.EQ.1) SIMMIX=1",
                      "\tIF (MIXNUM.EQ.2) POPCL=TVCL*RCL",
                      "\tIF (MIXNUM.EQ.2) SIMMIX=2",
                      "\tCL=POPCL*EXP(SD_ETA_CL*ETA_CL)",
                      "\tSIM_CL=CL",
                      ";VOLUME_OF_DISTRIBUTION",
                      "\tPOPVD=TVVD",
                      "\tV=POPVD*EXP(SD_ETA_VD*ETA_VD)",
                      "\tS2=V",
                      "\tSIM_VD=V",
                      ";ABSORPTION_RATE_CONSTANT",
                      "\tPOPKA=TVKA",
                      "\tKA=POPKA*EXP(SD_ETA_KA*ETA_KA)",
                      "\tSIM_KA=KA",
                      "",
                      "$MIX",
                      "MIXP1=THETA(8)",
                      "NSPOP=2",
                      "P(1)=MIXP1",
                      "P(2)=1-P(1)",
                      "",
                      "$THETA",
                      "REPLACE_TVCL ;TVCL",
                      "REPLACE_RCL ;RCL",
                      "REPLACE_CVCL ;CVCL",
                      "REPLACE_TVVD ;TVVD",
                      "REPLACE_CVVD ;CVVD",
                      "REPLACE_TVKA ;TVKA",
                      "REPLACE_CVKA ;CVKA",
                      "REPLACE_MIXP ;MIXP",
                      "REPLACE_W_PROP ;W_PROP",
                      "",
                      "$ERROR",
                      "IF (ICALL.EQ.4.AND.NEWIND.NE.2) THEN",
                      "\tDO WHILE (W_PROP*EPS(1).LE.-1)",
                      "\t\tCALL SIMEPS(EPS)",
                      "\tENDDO",
                      "ENDIF",
                      "Y=F*(1+W_PROP*EPS(1))",
                      "IPRED=F",
                      "",
                      "$OMEGA",
                      "1 FIXED ;ETA_CL",
                      "1 FIXED ;ETA_VD",
                      "1 FIXED ;ETA_KA",
                      "",
                      "$SIGMA",
                      "1 FIXED ;EPS_PROP",
                      "",
                      "$TABLE ID AMT TIME DV SIMMIX POPCL SIM_CL POPVD SIM_VD POPKA SIM_KA IPRED NOPRINT NOAPPEND ONEHEADER FILE=REPLACE_DATA_NAME"
  )
  
  #### Defining Fiting mixture control stream format ####
  
  fitMixStrTemp <- paste(sep = "\n",
                         "$PROBLEM REPLACE_PROB_NAME",
                         "",
                         "$INPUT ID AMT TIME DV SIMMIX POPCL=DROP SIM_CL POPVD=DROP SIM_VD POPKA=DROP SIM_KA SIM_IPRED",
                         "",
                         "$DATA REPLACE_DATA_NAME IGNORE=@",
                         "",
                         "$SUBROUTINE ADVAN2 TRANS2",
                         "",
                         "$PK",
                         ";THETA",
                         "\tTVCL=THETA(1)",
                         "\tRCL=THETA(2)",
                         "\tCVCL=THETA(3)",
                         "\tTVVD=THETA(4)",
                         "\tCVVD=THETA(5)",
                         "\tTVKA=THETA(6)",
                         "\tCVKA=THETA(7)",
                         "\t;MIXP: SEE $MIX BLOCK",
                         "\tW_PROP=THETA(9)",
                         ";ETA",
                         "\tETA_CL=ETA(1)",
                         "\tETA_VD=ETA(2)",
                         "\tETA_KA=ETA(3)",
                         ";TRANSFORMATION",
                         "\tSD_ETA_CL=SQRT(LOG(CVCL**2+1))",
                         "\tSD_ETA_VD=SQRT(LOG(CVVD**2+1))",
                         "\tSD_ETA_KA=SQRT(LOG(CVKA**2+1))",
                         ";CLEARANCE",
                         "\tIF (MIXNUM.EQ.1) POPCL=TVCL",
                         "\tIF (MIXNUM.EQ.1) MIX=1",
                         "\tIF (MIXNUM.EQ.2) POPCL=TVCL*RCL",
                         "\tIF (MIXNUM.EQ.2) MIX=2",
                         "\tCL=POPCL*EXP(SD_ETA_CL*ETA_CL)",
                         ";VOLUME_OF_DISTRIBUTION",
                         "\tPOPVD=TVVD",
                         "\tV=POPVD*EXP(SD_ETA_VD*ETA_VD)",
                         "\tS2=V",
                         "\tVD=V",
                         ";ABSORPTION_RATE_CONSTANT",
                         "\tPOPKA=TVKA",
                         "\tKA=POPKA*EXP(SD_ETA_KA*ETA_KA)",
                         "",
                         "$MIX",
                         "LOGIT_MIXP=THETA(8)",
                         "MIXP1=EXP(LOGIT_MIXP)/(1+EXP(LOGIT_MIXP))",
                         "NSPOP=2",
                         "P(1)=MIXP1",
                         "P(2)=1-P(1)",
                         "",
                         "$THETA",
                         "(0,REPLACE_IE_TVCL) ;TVCL",
                         "(1,REPLACE_IE_RCL) ;RCL",
                         "(0,REPLACE_IE_CVCL) ;CVCL",
                         "(0,REPLACE_IE_TVVD) ;TVVD",
                         "(0,REPLACE_IE_CVVD) ;CVVD",
                         "(0,REPLACE_IE_TVKA) ;TVKA",
                         "(0,REPLACE_IE_CVKA) ;CVKA",
                         "(-INF,REPLACE_IE_LOGIT_MIXP,INF) ;LOGIT_MIXP",
                         "(0,REPLACE_IE_W_PROP) ;W_PROP",
                         "",
                         "$ERROR",
                         "Y=F*(1+W_PROP*EPS(1))",
                         "IPRED=F",
                         "",
                         "$OMEGA",
                         "1 FIXED ;ETA_CL",
                         "1 FIXED ;ETA_VD",
                         "1 FIXED ;ETA_KA",
                         "",
                         "$SIGMA",
                         "1 FIXED ;EPS_PROP",
                         "",
                         "$ESTIMATION METHOD=1 INTERACTION MAXEVAL=9999 PRINT=1 NOABORT",
                         "",
                         "$COVARIANCE PRINT=E",
                         "",
                         "$TABLE ID AMT TIME DV SIMMIX MIX SIM_CL CL SIM_VD VD SIM_KA KA SIM_IPRED IPRED CWRES ECWRES EWRES NWRES UNCONDITIONAL NOPRINT ONEHEADER FILE=REPLACE_OUTPUT_NAME"
  )
  
  #### Defining Fiting without mixture control stream format ####
  
  fitNoMixStrTemp <- paste(sep = "\n",
                           "$PROBLEM REPLACE_PROB_NAME",
                           "",
                           "$INPUT ID AMT TIME DV SIMMIX POPCL=DROP SIM_CL POPVD=DROP SIM_VD POPKA=DROP SIM_KA SIM_IPRED",
                           "",
                           "$DATA REPLACE_DATA_NAME IGNORE=@",
                           "",
                           "$SUBROUTINE ADVAN2 TRANS2",
                           "",
                           "$PK",
                           ";THETA",
                           "\tTVCL=THETA(1)",
                           "\tRCL=THETA(2)",
                           "\tCVCL=THETA(3)",
                           "\tTVVD=THETA(4)",
                           "\tCVVD=THETA(5)",
                           "\tTVKA=THETA(6)",
                           "\tCVKA=THETA(7)",
                           "\t;MIXP=THETA(8)",
                           "\tW_PROP=THETA(9)",
                           ";ETA",
                           "\tETA_CL=ETA(1)",
                           "\tETA_VD=ETA(2)",
                           "\tETA_KA=ETA(3)",
                           ";TRANSFORMATION",
                           "\tSD_ETA_CL=SQRT(LOG(CVCL**2+1))",
                           "\tSD_ETA_VD=SQRT(LOG(CVVD**2+1))",
                           "\tSD_ETA_KA=SQRT(LOG(CVKA**2+1))",
                           ";CLEARANCE",
                           "\tPOPCL=TVCL",
                           "\tCL=POPCL*EXP(SD_ETA_CL*ETA_CL)",
                           ";VOLUME_OF_DISTRIBUTION",
                           "\tPOPVD=TVVD",
                           "\tV=POPVD*EXP(SD_ETA_VD*ETA_VD)",
                           "\tS2=V",
                           "\tVD=V",
                           ";ABSORPTION_RATE_CONSTANT",
                           "\tPOPKA=TVKA",
                           "\tKA=POPKA*EXP(SD_ETA_KA*ETA_KA)",
                           "",
                           "$THETA",
                           "(0,REPLACE_IE_TVCL) ;TVCL",
                           "1 FIXED ;RCL",
                           "(0,REPLACE_IE_CVCL) ;CVCL",
                           "(0,REPLACE_IE_TVVD) ;TVVD",
                           "(0,REPLACE_IE_CVVD) ;CVVD",
                           "(0,REPLACE_IE_TVKA) ;TVKA",
                           "(0,REPLACE_IE_CVKA) ;CVKA",
                           "0 FIXED ;MIXP",
                           "(0,REPLACE_IE_W_PROP) ;W_PROP",
                           "",
                           "$ERROR",
                           "Y=F*(1+W_PROP*EPS(1))",
                           "IPRED=F",
                           "",
                           "$OMEGA",
                           "1 FIXED ;ETA_CL",
                           "1 FIXED ;ETA_VD",
                           "1 FIXED ;ETA_KA",
                           "",
                           "$SIGMA",
                           "1 FIXED ;EPS_PROP",
                           "",
                           "$ESTIMATION METHOD=1 INTERACTION MAXEVAL=9999 PRINT=1 NOABORT",
                           "",
                           "$COVARIANCE PRINT=E",
                           "",
                           "$TABLE ID AMT TIME DV SIMMIX SIM_CL CL SIM_VD VD SIM_KA KA SIM_IPRED IPRED CWRES ECWRES EWRES NWRES UNCONDITIONAL NOPRINT ONEHEADER FILE=REPLACE_OUTPUT_NAME"
  )
  
  #### Defining Fiting mixture (true value) control stream format ####
  
  fitMixTVStrTemp <- paste(sep = "\n",
                           "$PROBLEM REPLACE_PROB_NAME",
                           "",
                           "$INPUT ID AMT TIME DV SIMMIX POPCL=DROP SIM_CL POPVD=DROP SIM_VD POPKA=DROP SIM_KA SIM_IPRED",
                           "",
                           "$DATA REPLACE_DATA_NAME IGNORE=@",
                           "",
                           "$SUBROUTINE ADVAN2 TRANS2",
                           "",
                           "$PK",
                           ";THETA",
                           "\tTVCL=THETA(1)",
                           "\tRCL=THETA(2)",
                           "\tCVCL=THETA(3)",
                           "\tTVVD=THETA(4)",
                           "\tCVVD=THETA(5)",
                           "\tTVKA=THETA(6)",
                           "\tCVKA=THETA(7)",
                           "\t;MIXP: SEE $MIX BLOCK",
                           "\tW_PROP=THETA(9)",
                           ";ETA",
                           "\tETA_CL=ETA(1)",
                           "\tETA_VD=ETA(2)",
                           "\tETA_KA=ETA(3)",
                           ";TRANSFORMATION",
                           "\tSD_ETA_CL=SQRT(LOG(CVCL**2+1))",
                           "\tSD_ETA_VD=SQRT(LOG(CVVD**2+1))",
                           "\tSD_ETA_KA=SQRT(LOG(CVKA**2+1))",
                           ";CLEARANCE",
                           "\tIF (MIXNUM.EQ.1) POPCL=TVCL",
                           "\tIF (MIXNUM.EQ.1) MIX=1",
                           "\tIF (MIXNUM.EQ.2) POPCL=TVCL*RCL",
                           "\tIF (MIXNUM.EQ.2) MIX=2",
                           "\tCL=POPCL*EXP(SD_ETA_CL*ETA_CL)",
                           ";VOLUME_OF_DISTRIBUTION",
                           "\tPOPVD=TVVD",
                           "\tV=POPVD*EXP(SD_ETA_VD*ETA_VD)",
                           "\tS2=V",
                           "\tVD=V",
                           ";ABSORPTION_RATE_CONSTANT",
                           "\tPOPKA=TVKA",
                           "\tKA=POPKA*EXP(SD_ETA_KA*ETA_KA)",
                           "",
                           "$MIX",
                           "MIXP1=THETA(8)",
                           "NSPOP=2",
                           "P(1)=MIXP1",
                           "P(2)=1-P(1)",
                           "",
                           "$THETA",
                           "REPLACE_TV_TVCL ;TVCL",
                           "REPLACE_TV_RCL ;RCL",
                           "REPLACE_TV_CVCL ;CVCL",
                           "REPLACE_TV_TVVD ;TVVD",
                           "REPLACE_TV_CVVD ;CVVD",
                           "REPLACE_TV_TVKA ;TVKA",
                           "REPLACE_TV_CVKA ;CVKA",
                           "REPLACE_TV_MIXP ;MIXP",
                           "REPLACE_TV_W_PROP ;W_PROP",
                           "",
                           "$ERROR",
                           "Y=F*(1+W_PROP*EPS(1))",
                           "IPRED=F",
                           "",
                           "$OMEGA",
                           "1 FIXED ;ETA_CL",
                           "1 FIXED ;ETA_VD",
                           "1 FIXED ;ETA_KA",
                           "",
                           "$SIGMA",
                           "1 FIXED ;EPS_PROP",
                           "",
                           "$ESTIMATION METHOD=1 INTERACTION MAXEVAL=0 PRINT=1 NOABORT",
                           "",
                           "$COVARIANCE PRINT=E",
                           "",
                           "$TABLE ID AMT TIME DV SIMMIX MIX SIM_CL CL SIM_VD VD SIM_KA KA SIM_IPRED IPRED CWRES ECWRES EWRES NWRES UNCONDITIONAL NOPRINT ONEHEADER FILE=REPLACE_OUTPUT_NAME"
  )
  
  #### Backup settings ####
  dput(list_for_dput(c(
    "dose", "subDirSize",
    "sim_nmfe_options", "sim_clean", "sim_nm_output", "sim_threads",
    "fitMix_nmfe_options", "fitMix_clean", "fitMix_nm_output", "fitMix_threads",
    "fitMix_min_retries", "fitMix_nodes", "fitMix_picky", "fitMix_retries",
    "fitMix_tweak_inits",
    "fitNoMix_nmfe_options", "fitNoMix_clean", "fitNoMix_nm_output", "fitNoMix_threads",
    "fitNoMix_min_retries", "fitNoMix_nodes", "fitNoMix_picky", "fitNoMix_retries",
    "fitNoMix_tweak_inits",
    "fitMixTV_nmfe_options", "fitMixTV_clean", "fitMixTV_nm_output", "fitMixTV_threads",
    "fitMixTV_min_retries", "fitMixTV_nodes",
    "IE_fix_minDev", "IE_fix_maxDev", "IE_LOGIT_MIXP_minDev", "IE_LOGIT_MIXP_maxDev",
    "IE_rnd_minDev", "IE_rnd_maxDev",
    "tmaxFactor", "t12Factor",
    "paramNames", "paramRngs",
    "NdatPerScn",
    "simStrTemp", "fitMixStrTemp", "fitNoMixStrTemp", "fitMixTVStrTemp"
  )), paste0(projName, "_SETTING.R"))
  
  cat("\rReading settings... Done!")
} else {
  setwd(baseWd)
  cat("\n\tPrevious settings found! Previous settings will be read...")
  setting_FileName <- regexpr("SETTING.R", dir()) %>% equals(-1) %>% not %>% which %>% dir()[.]
  readProjName <- setting_FileName %>% substr(start = 1, stop = nchar(.) - 10)
  if (projName != readProjName) {
    cat("\n\tProject name specified does not match previous setting! Previous name will be used...")
    projName <- readProjName
  }
  dgetList_to_var(dget(setting_FileName))
  cat(" Done!")
}

#++++++++++++++++++++++++++++++#
#### +++++Summarizing+++++ #####
#++++++++++++++++++++++++++++++#

cat("\nSummarizing settings...")

Nparam <- length(paramRngs) #Number of parameters to sample
Nrngs <- sapply(paramRngs, length) %>% unname #Number of ranges for each parameter
Nscn <- prod(Nrngs) #Number of scenarios
Ndat <- Nscn * NdatPerScn #Number of datasets

toRun <- if (is.null(testRun)) {seq(Ndat)} else {testRun}

SerialFormat <- log10(Ndat) %>% ceiling %>% sprintf("%02d", .) %>% paste0("%", ., "d") #Formatting for serial
simProbNameArr <- paste0("SIM_", projName, "_", sprintf(SerialFormat, seq(Ndat))) #Array of simulation stream problem names
fitMixProbNameArr <- paste0("FIT_MIX_", projName, "_", sprintf(SerialFormat, seq(Ndat))) #Array of fit mixture stream problem names
fitNoMixProbNameArr <- paste0("FIT_NOMIX_", projName, "_", sprintf(SerialFormat, seq(Ndat))) #Array of fit no mixture stream problem names
fitMixTVProbNameArr <- paste0("FIT_MIXTV_", projName, "_", sprintf(SerialFormat, seq(Ndat))) #Array of fit mixture (true value) stream problem names
datTempNameArr <- paste0("DATATEMP_", projName, "_", sprintf(SerialFormat, seq(Ndat)), ".csv") #Array of dataset template names
datNameArr <- paste0("DATA_", projName, "_", sprintf(SerialFormat, seq(Ndat)), ".txt") #Array of dataset names
simStrNameArr <- paste0("SIM_", projName, "_", sprintf(SerialFormat, seq(Ndat)), ".mod") #Array of simulation stream names
simLstNameArr <- paste0("SIM_", projName, "_", sprintf(SerialFormat, seq(Ndat)), ".lst") #Array of simulation lst names
fitMixStrNameArr <- paste0("FIT_MIX_", projName, "_", sprintf(SerialFormat, seq(Ndat)), ".mod") #Array of fit mixture stream names
fitMixLstNameArr <- paste0("FIT_MIX_", projName, "_", sprintf(SerialFormat, seq(Ndat)), ".lst") #Array of fit mixture lst names
fitMixTxtNameArr <- paste0("FIT_MIX_", projName, "_", sprintf(SerialFormat, seq(Ndat)), ".txt") #Array of fit mixture txt names
fitMixExtNameArr <- paste0("FIT_MIX_", projName, "_", sprintf(SerialFormat, seq(Ndat)), ".ext") #Array of fit mixture ext names
fitMixShkNameArr <- paste0("FIT_MIX_", projName, "_", sprintf(SerialFormat, seq(Ndat)), ".shk") #Array of fit mixture shk names
fitMixPhiNameArr <- paste0("FIT_MIX_", projName, "_", sprintf(SerialFormat, seq(Ndat)), ".phi") #Array of fit mixture phi names
fitMixPhmNameArr <- paste0("FIT_MIX_", projName, "_", sprintf(SerialFormat, seq(Ndat)), ".phm") #Array of fit mixture phm names
fitNoMixStrNameArr <- paste0("FIT_NOMIX_", projName, "_", sprintf(SerialFormat, seq(Ndat)), ".mod") #Array of fit no mixture stream names
fitNoMixLstNameArr <- paste0("FIT_NOMIX_", projName, "_", sprintf(SerialFormat, seq(Ndat)), ".lst") #Array of fit no mixture lst names
fitNoMixTxtNameArr <- paste0("FIT_NOMIX_", projName, "_", sprintf(SerialFormat, seq(Ndat)), ".txt") #Array of fit no mixture txt names
fitNoMixExtNameArr <- paste0("FIT_NOMIX_", projName, "_", sprintf(SerialFormat, seq(Ndat)), ".ext") #Array of fit no mixture ext names
fitNoMixShkNameArr <- paste0("FIT_NOMIX_", projName, "_", sprintf(SerialFormat, seq(Ndat)), ".shk") #Array of fit no mixture shk names
fitNoMixPhiNameArr <- paste0("FIT_NOMIX_", projName, "_", sprintf(SerialFormat, seq(Ndat)), ".phi") #Array of fit no mixture phi names
fitMixTVStrNameArr <- paste0("FIT_MIXTV_", projName, "_", sprintf(SerialFormat, seq(Ndat)), ".mod") #Array of fit mixture (true value) stream names
fitMixTVLstNameArr <- paste0("FIT_MIXTV_", projName, "_", sprintf(SerialFormat, seq(Ndat)), ".lst") #Array of fit mixture (true value) lst names
fitMixTVTxtNameArr <- paste0("FIT_MIXTV_", projName, "_", sprintf(SerialFormat, seq(Ndat)), ".txt") #Array of fit mixture (true value) txt names
fitMixTVExtNameArr <- paste0("FIT_MIXTV_", projName, "_", sprintf(SerialFormat, seq(Ndat)), ".ext") #Array of fit mixture (true value) lst names
fitMixTVPhiNameArr <- paste0("FIT_MIXTV_", projName, "_", sprintf(SerialFormat, seq(Ndat)), ".phi") #Array of fit mixture (true value) phi names
fitMixTVPhmNameArr <- paste0("FIT_MIXTV_", projName, "_", sprintf(SerialFormat, seq(Ndat)), ".phm") #Array of fit mixture (true value) phm names

#Array of names of subdirectories
subDirLL <- seq(Ndat) %>% subtract(1) %>% divide_by(100) %>% floor %>% multiply_by(100) %>% add(1)
subDirUL <- seq(Ndat) %>% divide_by(100) %>% ceiling %>% multiply_by(100)
subDirUL[subDirUL > Ndat] <- Ndat
subDirArr <- paste0(sprintf(SerialFormat, subDirLL), "-", sprintf(SerialFormat, subDirUL))

#Paths of files
datTempPathArr <- paste0(subDirArr, "/", datTempNameArr)
datPathArr <- paste0(subDirArr, "/", datNameArr)
datPathArrBS <- paste0(subDirArr, "\\", datNameArr)
simStrPathArr <- paste0(subDirArr, "/", simStrNameArr)
simLstPathArr <- paste0(subDirArr, "/", simLstNameArr)
fitMixStrPathArr <- paste0(subDirArr, "/", fitMixStrNameArr)
fitMixLstPathArr <- paste0(subDirArr, "/", fitMixLstNameArr)
fitMixTxtPathArr <- paste0(subDirArr, "/", fitMixTxtNameArr)
fitMixExtPathArr <- paste0(subDirArr, "/", fitMixExtNameArr)
fitMixShkPathArr <- paste0(subDirArr, "/", fitMixShkNameArr)
fitMixPhiPathArr <- paste0(subDirArr, "/", fitMixPhiNameArr)
fitMixPhmPathArr <- paste0(subDirArr, "/", fitMixPhmNameArr)
fitNoMixStrPathArr <- paste0(subDirArr, "/", fitNoMixStrNameArr)
fitNoMixLstPathArr <- paste0(subDirArr, "/", fitNoMixLstNameArr)
fitNoMixTxtPathArr <- paste0(subDirArr, "/", fitNoMixTxtNameArr)
fitNoMixExtPathArr <- paste0(subDirArr, "/", fitNoMixExtNameArr)
fitNoMixShkPathArr <- paste0(subDirArr, "/", fitNoMixShkNameArr)
fitNoMixPhiPathArr <- paste0(subDirArr, "/", fitNoMixPhiNameArr)
fitMixTVStrPathArr <- paste0(subDirArr, "/", fitMixTVStrNameArr)
fitMixTVLstPathArr <- paste0(subDirArr, "/", fitMixTVLstNameArr)
fitMixTVTxtPathArr <- paste0(subDirArr, "/", fitMixTVTxtNameArr)
fitMixTVExtPathArr <- paste0(subDirArr, "/", fitMixTVExtNameArr)
fitMixTVPhiPathArr <- paste0(subDirArr, "/", fitMixTVPhiNameArr)
fitMixTVPhmPathArr <- paste0(subDirArr, "/", fitMixTVPhmNameArr)

#Spliting simulation stream template
simSplitLst <- c("REPLACE_PROB_NAME",
                 "REPLACE_DATA_TEMPLATE_NAME",
                 "REPLACE_SEED",
                 "REPLACE_TVCL",
                 "REPLACE_RCL",
                 "REPLACE_CVCL",
                 "REPLACE_TVVD",
                 "REPLACE_CVVD",
                 "REPLACE_TVKA",
                 "REPLACE_CVKA",
                 "REPLACE_MIXP",
                 "REPLACE_W_PROP",
                 "REPLACE_DATA_NAME")
for (i in seq(length(simSplitLst))) {
  if (i == 1) {Left <- simStrTemp; simStrTempSplit <- character(0)}
  splitted <- strsplit(Left, simSplitLst[i]) %>% unlist
  if (!is.na(splitted[1])) {simStrTempSplit %<>% c(., splitted[1])}
  Left <- splitted[2]
  if (i == length(simSplitLst) & (!is.na(splitted[2]))) {simStrTempSplit %<>% c(., splitted[2])}
}

#Spliting fit mixture stream template
fitMixSplitLst <- c("REPLACE_PROB_NAME",
                    "REPLACE_DATA_NAME",
                    "REPLACE_IE_TVCL",
                    "REPLACE_IE_RCL",
                    "REPLACE_IE_CVCL",
                    "REPLACE_IE_TVVD",
                    "REPLACE_IE_CVVD",
                    "REPLACE_IE_TVKA",
                    "REPLACE_IE_CVKA",
                    "REPLACE_IE_LOGIT_MIXP",
                    "REPLACE_IE_W_PROP",
                    "REPLACE_OUTPUT_NAME")
for (i in seq(length(fitMixSplitLst))) {
  if (i == 1) {Left <- fitMixStrTemp; fitMixStrTempSplit <- character(0)}
  splitted <- strsplit(Left, fitMixSplitLst[i]) %>% unlist
  if (!is.na(splitted[1])) {fitMixStrTempSplit %<>% c(., splitted[1])}
  Left <- splitted[2]
  if (i == length(fitMixSplitLst) & (!is.na(splitted[2]))) {fitMixStrTempSplit %<>% c(., splitted[2])}
}

#Spliting fit no mixture stream template
fitNoMixSplitLst <- c("REPLACE_PROB_NAME",
                      "REPLACE_DATA_NAME",
                      "REPLACE_IE_TVCL",
                      "REPLACE_IE_CVCL",
                      "REPLACE_IE_TVVD",
                      "REPLACE_IE_CVVD",
                      "REPLACE_IE_TVKA",
                      "REPLACE_IE_CVKA",
                      "REPLACE_IE_W_PROP",
                      "REPLACE_OUTPUT_NAME")
for (i in seq(length(fitNoMixSplitLst))) {
  if (i == 1) {Left <- fitNoMixStrTemp; fitNoMixStrTempSplit <- character(0)}
  splitted <- strsplit(Left, fitNoMixSplitLst[i]) %>% unlist
  if (!is.na(splitted[1])) {fitNoMixStrTempSplit %<>% c(., splitted[1])}
  Left <- splitted[2]
  if (i == length(fitNoMixSplitLst) & (!is.na(splitted[2]))) {fitNoMixStrTempSplit %<>% c(., splitted[2])}
}

#Spliting fit mixture stream template
fitMixTVSplitLst <- c("REPLACE_PROB_NAME",
                      "REPLACE_DATA_NAME",
                      "REPLACE_TV_TVCL",
                      "REPLACE_TV_RCL",
                      "REPLACE_TV_CVCL",
                      "REPLACE_TV_TVVD",
                      "REPLACE_TV_CVVD",
                      "REPLACE_TV_TVKA",
                      "REPLACE_TV_CVKA",
                      "REPLACE_TV_MIXP",
                      "REPLACE_TV_W_PROP",
                      "REPLACE_OUTPUT_NAME")
for (i in seq(length(fitMixTVSplitLst))) {
  if (i == 1) {Left <- fitMixTVStrTemp; fitMixTVStrTempSplit <- character(0)}
  splitted <- strsplit(Left, fitMixTVSplitLst[i]) %>% unlist
  if (!is.na(splitted[1])) {fitMixTVStrTempSplit %<>% c(., splitted[1])}
  Left <- splitted[2]
  if (i == length(fitMixTVSplitLst) & (!is.na(splitted[2]))) {fitMixTVStrTempSplit %<>% c(., splitted[2])}
}

#PsN prefixes
sim_psn_prefix <- paste0("execute -threads=", sim_threads, if (sim_nm_output == "") {""} else {paste0(" -nm_output=", sim_nm_output)},
                     " -clean=", sim_clean, if (sim_nmfe_options == "") {""} else {paste0(" -nmfe_options=", sim_nmfe_options)})
fitMix_psn_prefix <- paste0("execute -threads=", fitMix_threads,
                            if (fitMix_nm_output == "") {""} else {paste0(" -nm_output=", fitMix_nm_output)},
                            " -clean=", fitMix_clean,
                            if (fitMix_nmfe_options == "") {""} else {paste0(" -nmfe_options=", fitMix_nmfe_options)},
                            " -min_retries=", fitMix_min_retries,
                            if (fitMix_nodes == "") {""} else {paste0(" -nodes=", fitMix_nodes)},
                            if (fitMix_picky) {" -picky"} else {""},
                            " -retries=", fitMix_retries,
                            if (fitMix_tweak_inits) {" -tweak_inits"} else {""}
                            )
fitNoMix_psn_prefix <- paste0("execute -threads=", fitNoMix_threads,
                            if (fitNoMix_nm_output == "") {""} else {paste0(" -nm_output=", fitNoMix_nm_output)},
                            " -clean=", fitNoMix_clean,
                            if (fitNoMix_nmfe_options == "") {""} else {paste0(" -nmfe_options=", fitNoMix_nmfe_options)},
                            " -min_retries=", fitNoMix_min_retries,
                            if (fitNoMix_nodes == "") {""} else {paste0(" -nodes=", fitNoMix_nodes)},
                            if (fitNoMix_picky) {" -picky"} else {""},
                            " -retries=", fitNoMix_retries,
                            if (fitNoMix_tweak_inits) {" -tweak_inits"} else {""}
                            )
fitMixTV_psn_prefix <- paste0("execute -threads=", fitMixTV_threads,
                              if (fitMixTV_nm_output == "") {""} else {paste0(" -nm_output=", fitMixTV_nm_output)},
                              " -clean=", fitMixTV_clean,
                              if (fitMixTV_nmfe_options == "") {""} else {paste0(" -nmfe_options=", fitMixTV_nmfe_options)},
                              " -min_retries=", fitMixTV_min_retries,
                              if (fitMixTV_nodes == "") {""} else {paste0(" -nodes=", fitMixTV_nodes)}
)
cat(" Done!")

#++++++++++++++++++++++++++++++++++++++++#
#### +++++Simulating parameters+++++ #####
#++++++++++++++++++++++++++++++++++++++++#

cat("\nSimulating parameters...")

dir.create(baseWd, showWarnings = F, recursive = T)
setwd(baseWd)

#### Simulating population parameters ####

times <- c(1, Nrngs[-Nparam]) %>% cumprod
each <- c(1, Nrngs[Nparam:2]) %>% cumprod %>% sort(decreasing = T)
if (!file.exists(paste0(projName, "_SIMULATED_PARAMETERS.csv"))) {
  # Simulate only if not simulated before
  dt_param <- sapply(seq(Nparam), function(i) {
    dt_lim = paramRngs[[i]] %>% unlist %>% matrix(ncol = 2, byrow = T) %>% data.table %>% set_colnames(c("LL", "UL"))
    LL_arr = dt_lim[, LL][rep(Nrngs[i] %>% seq, times = times[i], each = each[i])] %>% rep(each = NdatPerScn)
    UL_arr = dt_lim[, UL][rep(Nrngs[i] %>% seq, times = times[i], each = each[i])] %>% rep(each = NdatPerScn)
    param_arr = (LL_arr + (UL_arr - LL_arr) * runif(Ndat)) %>% signif(digits = 4)
    if (paramNames[i] == "SampleSize") {param_arr %<>% round}
    return(param_arr)
  }) %>% as.data.table %>% set_colnames(paramNames)
  write_csv(dt_param, paste0(projName, "_SIMULATED_PARAMETERS.csv"))
} else {
  # If simulated before, read previous simulated values
  cat("\n\tPopulation parameters had been simulated! Reading simulated values...")
  dt_param <- fread(paste0(projName, "_SIMULATED_PARAMETERS.csv"), header = T)
}

#### Simulating seed numbers ####

if (!file.exists(paste0(projName, "_SIMULATED_SEED_NUMBERS.csv"))) {
  # Simulate only if not simulated before
  seedArr <- runif(Ndat, max = .Machine$integer.max) %>% ceiling
  write_csv(data.table(seedArr), paste0(projName, "_SIMULATED_SEED_NUMBERS.csv"))
} else {
  # If simulated before, read previous simulated values
  cat("\n\tSeed Numbers had been simulated! Reading simulated values...")
  seedArr <- fread(paste0(projName, "_SIMULATED_SEED_NUMBERS.csv"), header = T) %>% unlist %>% unname
}

#### Calculating sampling time according to simulated population parameters ####

if (!file.exists(paste0(projName, "_SAMPLING_TIMES.csv"))) {
  # Simulate only if not simulated before
  if (length(tmaxFactor) != length(t12Factor)) {stop("tmaxFactor and t12Factor have different lengths!")} else {
    Ntime <- length(tmaxFactor)  
  }
  tmaxArr <- dt_param[, log(TVKA/((TVCL*MIXP + TVCL*RCL*(1-MIXP))/TVVD))/(TVKA - (TVCL*MIXP + TVCL*RCL*(1-MIXP))/TVVD)]
  t12Arr <- dt_param[, log(2)/((TVCL*MIXP + TVCL*RCL*(1-MIXP))/TVVD)]
  matTime <- matrix(
    (rep(tmaxFactor, each = Ndat) * rep(tmaxArr, Ntime) + rep(t12Factor, each = Ndat) * rep(t12Arr, Ntime)) %>%
      round(1),
    ncol = Ntime
  )
  write_csv(as.data.table(matTime), paste0(projName, "_SAMPLING_TIMES.csv"))
} else {
  # If simulated before, read previous simulated values
  cat("\n\tSampling times had been calculated before! Reading calculated values...")
  Ntime <- length(tmaxFactor) 
  matTime <- fread(paste0(projName, "_SAMPLING_TIMES.csv"), header = T) %>% as.matrix
}
  
cat(" Done!")

cat("\n\n")
cat("********************************************\n")
cat("******** Generating necessary files ********\n")
cat("********************************************\n")

#+++++++++++++++++++++++++++++++++++#
#### +++++Generating files+++++ #####
#+++++++++++++++++++++++++++++++++++#

#### Generating dataset templates ####

cat("\nGenerating dataset templates...")
if (!skip_generating_data_temp) {
  dir.create(paste0(baseWd, "/DATA"), showWarnings = F, recursive = T)
  setwd(paste0(baseWd, "/DATA"))
  cat("\n\tChecking progress and deciding batch to run...")
  dataTempProg <- which(!file.exists(datTempPathArr[toRun]))
  if (dataTempProg %>% length %>% is_greater_than(0)) {
    cat("\n")
    dataTemp_start_Time <- Sys.time()
    j <- 0
    for (i in toRun[dataTempProg]) {
      j <- j + 1
      thisN = dt_param[i, SampleSize]
      # Generating dataset template table
      tab_out = data.table(ID = rep(seq(thisN), each = Ntime + 1),
                           AMT = rep(c(dose, rep(0, Ntime)), times = thisN),
                           TIME = rep(c(0, matTime[i,]), times = thisN),
                           DV = rep(0, thisN * (Ntime + 1)))
      # Creating non-existent subdirectory
      thisDir = paste0(baseWd, "/DATA/", subDirArr[i])
      if (!dir.exists(thisDir)) {dir.create(thisDir, recursive = T)}
      # Writing dataset template as csv
      write_csv(tab_out, datTempPathArr[i])
      # Progress report
      dataTemp_End_Time <- Sys.time()
      if (j%%100 == 0 || j == length(dataTempProg)) {
        cat(paste0("\r\t", strrep(" ", nchar(displayStr) - 2)))
        displayStr <- paste0("\r\tGenerating dataset templates... ", j, "/", length(dataTempProg),
                             if (j == length(dataTempProg)) {"  |  Done!"} else {""},
                             "  |  Elapsed Time: ",
                             diffTime_to_display(dataTemp_End_Time - dataTemp_start_Time))
        cat(displayStr)
      }
    }
  } else {
    cat(" Done previously!")
  }
} else {
  cat(" Process skipped!")
}

#### Generating simulation control streams ####

cat("\nGenerating simulation control streams...")

if (!skip_generating_sim_stream) {
  # Creating strings
  simStrArr <- paste0(
    if (is.na(simStrTempSplit[1])) {""} else {simStrTempSplit[1]}, simProbNameArr,
    simStrTempSplit[2], datTempNameArr, simStrTempSplit[3], seedArr,
    simStrTempSplit[4], dt_param[, TVCL], simStrTempSplit[5], dt_param[, RCL], simStrTempSplit[6], dt_param[, CVCL],
    simStrTempSplit[7], dt_param[, TVVD], simStrTempSplit[8], dt_param[, CVVD],
    simStrTempSplit[9], dt_param[, TVKA], simStrTempSplit[10], dt_param[, CVKA],
    simStrTempSplit[11], dt_param[, MIXP], simStrTempSplit[12], dt_param[, W_PROP],
    simStrTempSplit[13], datNameArr,
    if (is.na(simStrTempSplit[14])) {""} else {simStrTempSplit[14]}
  )
  # Generating files
  setwd(paste0(baseWd, "/DATA"))
  cat("\n\tChecking progress and deciding batch to run...")
  simStrProg <- which(!file.exists(simStrPathArr[toRun]))
  
  if (simStrProg %>% length %>% is_greater_than(0)) {
    cat("\n")
    simStr_start_Time <- Sys.time()
    j <- 0
    for (i in toRun[simStrProg]) {
      j <- j + 1
      write(simStrArr[i], simStrPathArr[i]) %>% invisible
      simStr_End_Time <- Sys.time()
      if (j%%100 == 0 || j == length(simStrProg)) {
        cat(paste0("\r\t", strrep(" ", nchar(displayStr) - 2)))
        displayStr <- paste0("\r\tGenerating simulation control streams... ", j, "/", length(simStrProg),
                             if (j == length(simStrProg)) {"  |  Done!"} else {""},
                             "  |  Elapsed Time: ",
                             diffTime_to_display(simStr_End_Time - simStr_start_Time))
        cat(displayStr)
      }
    }
  } else {
    cat(" Done previously!")
  }
} else {
  cat(" Process skipped!")
}


#### Generating fit mixture control streams ####

cat("\nGenerating fit mixture control streams...")

if (!skip_generating_fit_mix_stream) {
  # Generating initial estimates
  fitMix_IE_TVCL <- dt_param[, TVCL + sample(c(-1, 1), Ndat, replace = T) *
                               runif(Ndat, IE_fix_minDev, IE_fix_maxDev) * TVCL] %>% round(digits = 3)
  fitMix_IE_RCL <- dt_param[, RCL + sample(c(-1, 1), Ndat, replace = T) *
                              runif(Ndat, IE_fix_minDev, IE_fix_maxDev) * RCL] %>% round(digits = 3)
  fitMix_IE_RCL[fitMix_IE_RCL < 1.1] <- 1.1
  fitMix_IE_CVCL <- dt_param[, CVCL +
                               runif(Ndat, IE_rnd_minDev, IE_rnd_maxDev) * CVCL] %>% round(digits = 3)
  fitMix_IE_TVVD <- dt_param[, TVVD + sample(c(-1, 1), Ndat, replace = T) *
                               runif(Ndat, IE_fix_minDev, IE_fix_maxDev) * TVVD] %>% round(digits = 3)
  fitMix_IE_CVVD <- dt_param[, CVVD +
                               runif(Ndat, IE_rnd_minDev, IE_rnd_maxDev) * CVVD] %>% round(digits = 3)
  fitMix_IE_TVKA <- dt_param[, TVKA + sample(c(-1, 1), Ndat, replace = T) *
                               runif(Ndat, IE_fix_minDev, IE_fix_maxDev) * TVKA] %>% round(digits = 3)
  fitMix_IE_CVKA <- dt_param[, CVKA +
                               runif(Ndat, IE_rnd_minDev, IE_rnd_maxDev) * CVKA] %>% round(digits = 3)
  fitMix_IE_LOGIT_MIXP <- dt_param[, log(MIXP/(1-MIXP)) + sample(c(-1, 1), Ndat, replace = T) *
                                     runif(Ndat, IE_LOGIT_MIXP_minDev, IE_LOGIT_MIXP_maxDev)] %>% round(digits = 3)
  fitMix_IE_LOGIT_MIXP[abs(fitMix_IE_LOGIT_MIXP) < 0.01] <- 0.01
  fitMix_IE_W_PROP <- dt_param[, W_PROP +
                                 runif(Ndat, IE_rnd_minDev, IE_rnd_maxDev) * W_PROP] %>% round(digits = 3)
  
  # Creating strings
  fitMixStrArr <- paste0(
    if (is.na(fitMixStrTempSplit[1])) {""} else {fitMixStrTempSplit[1]}, fitMixProbNameArr,
    fitMixStrTempSplit[2],
    if (Sys.info()["sysname"] == "Windows") {
      paste0("..\\..\\DATA\\", datPathArrBS)
    } else {
      paste0("../../DATA/", datPathArr)
    },
    fitMixStrTempSplit[3], fitMix_IE_TVCL, fitMixStrTempSplit[4], fitMix_IE_RCL,
    fitMixStrTempSplit[5], fitMix_IE_CVCL,
    fitMixStrTempSplit[6], fitMix_IE_TVVD, fitMixStrTempSplit[7], fitMix_IE_CVVD,
    fitMixStrTempSplit[8], fitMix_IE_TVKA, fitMixStrTempSplit[9], fitMix_IE_CVKA,
    fitMixStrTempSplit[10], fitMix_IE_LOGIT_MIXP, fitMixStrTempSplit[11], fitMix_IE_W_PROP,
    fitMixStrTempSplit[12], fitMixTxtNameArr,
    if (is.na(fitMixStrTempSplit[13])) {""} else {fitMixStrTempSplit[13]}
  )
  # Generating files
  dir.create(paste0(baseWd, "/FIT_MIX"), showWarnings = F, recursive = T)
  setwd(paste0(baseWd, "/FIT_MIX"))
  cat("\n\tChecking progress and deciding batch to run...")
  fitMixStrProg <- which(!file.exists(fitMixStrPathArr[toRun]))
  
  if (fitMixStrProg %>% length %>% is_greater_than(0)) {
    cat("\n")
    fitMixStr_start_Time <- Sys.time()
    j <- 0
    for (i in toRun[fitMixStrProg]) {
      j <- j + 1
      # Creating non-existent subdirectory
      thisDir = paste0(baseWd, "/FIT_MIX/", subDirArr[i])
      if (!dir.exists(thisDir)) {dir.create(thisDir, recursive = T)}
      # Writing files
      write(fitMixStrArr[i], fitMixStrPathArr[i]) %>% invisible
      fitMixStr_End_Time <- Sys.time()
      if (j%%100 == 0 || j == length(fitMixStrProg)) {
        cat(paste0("\r\t", strrep(" ", nchar(displayStr) - 2)))
        displayStr <- paste0("\r\tGenerating fit mixture control streams... ", j, "/", length(fitMixStrProg),
                             if (j == length(fitMixStrProg)) {"  |  Done!"} else {""},
                             "  |  Elapsed Time: ",
                             diffTime_to_display(fitMixStr_End_Time - fitMixStr_start_Time))
        cat(displayStr)
      }
    }
  } else {
    cat(" Done previously!")
  }
} else {
  cat(" Process skipped!")
}

#### Generating fit without mixture control streams ####

cat("\nGenerating fit without mixture control streams...")

if (!skip_generating_fit_nomix_stream) {
  # Generating initial estimates
  fitNoMix_IE_TVCL <- dt_param[, (TVCL * MIXP + TVCL * RCL * (1 - MIXP)) +
                                 sample(c(-1, 1), Ndat, replace = T) *
                                 runif(Ndat, IE_fix_minDev, IE_fix_maxDev) *
                                 (TVCL * MIXP + TVCL * RCL * (1 - MIXP))] %>% round(digits = 3)
  fitNoMix_IE_CVCL <- dt_param[, CVCL +
                                 runif(Ndat, IE_rnd_minDev, IE_rnd_maxDev) * CVCL] %>% round(digits = 3)
  fitNoMix_IE_TVVD <- dt_param[, TVVD + sample(c(-1, 1), Ndat, replace = T) *
                                 runif(Ndat, IE_fix_minDev, IE_fix_maxDev) * TVVD] %>% round(digits = 3)
  fitNoMix_IE_CVVD <- dt_param[, CVVD +
                                 runif(Ndat, IE_rnd_minDev, IE_rnd_maxDev) * CVVD] %>% round(digits = 3)
  fitNoMix_IE_TVKA <- dt_param[, TVKA + sample(c(-1, 1), Ndat, replace = T) *
                                 runif(Ndat, IE_fix_minDev, IE_fix_maxDev) * TVKA] %>% round(digits = 3)
  fitNoMix_IE_CVKA <- dt_param[, CVKA +
                                 runif(Ndat, IE_rnd_minDev, IE_rnd_maxDev) * CVKA] %>% round(digits = 3)
  fitNoMix_IE_W_PROP <- dt_param[, W_PROP +
                                   runif(Ndat, IE_rnd_minDev, IE_rnd_maxDev) * W_PROP] %>% round(digits = 3)
  # Creating strings
  fitNoMixStrArr <- paste0(
    if (is.na(fitNoMixStrTempSplit[1])) {""} else {fitNoMixStrTempSplit[1]}, fitNoMixProbNameArr,
    fitNoMixStrTempSplit[2],
    if (Sys.info()["sysname"] == "Windows") {
      paste0("..\\..\\DATA\\", datPathArrBS)
    } else {
      paste0("../../DATA/", datPathArr)
    },
    fitNoMixStrTempSplit[3], fitNoMix_IE_TVCL, fitNoMixStrTempSplit[4], fitNoMix_IE_CVCL,
    fitNoMixStrTempSplit[5], fitNoMix_IE_TVVD, fitNoMixStrTempSplit[6], fitNoMix_IE_CVVD,
    fitNoMixStrTempSplit[7], fitNoMix_IE_TVKA, fitNoMixStrTempSplit[8], fitNoMix_IE_CVKA,
    fitNoMixStrTempSplit[9], fitNoMix_IE_W_PROP,
    fitNoMixStrTempSplit[10], fitNoMixTxtNameArr,
    if (is.na(fitNoMixStrTempSplit[11])) {""} else {fitNoMixStrTempSplit[11]}
  )
  # Generating files
  dir.create(paste0(baseWd, "/FIT_NOMIX"), showWarnings = F, recursive = T)
  setwd(paste0(baseWd, "/FIT_NOMIX"))
  cat("\n\tChecking progress and deciding batch to run...")
  fitNoMixStrProg <- which(!file.exists(fitNoMixStrPathArr[toRun]))
  
  if (fitNoMixStrProg %>% length %>% is_greater_than(0)) {
    cat("\n")
    fitNoMixStr_start_Time <- Sys.time()
    j <- 0
    for (i in toRun[fitNoMixStrProg]) {
      j <- j + 1
      # Creating non-existent subdirectory
      thisDir = paste0(baseWd, "/FIT_NOMIX/", subDirArr[i])
      if (!dir.exists(thisDir)) {dir.create(thisDir, recursive = T)}
      # Writing files
      write(fitNoMixStrArr[i], fitNoMixStrPathArr[i]) %>% invisible
      fitNoMixStr_End_Time <- Sys.time()
      if (j%%100 == 0 || j == length(fitNoMixStrProg)) {
        cat(paste0("\r\t", strrep(" ", nchar(displayStr) - 2)))
        displayStr <- paste0("\r\tGenerating fit without mixture control streams... ", j, "/", length(fitNoMixStrProg),
                             if (j == length(fitNoMixStrProg)) {"  |  Done!"} else {""},
                             "  |  Elapsed Time: ",
                             diffTime_to_display(fitNoMixStr_End_Time - fitNoMixStr_start_Time))
        cat(displayStr)
      }
    }
  } else {
    cat(" Done previously!")
  }
} else {
  cat(" Process skipped!")
}

#### Generating fit mixture (true value) control streams ####

cat("\nGenerating fit mixture (true value) control streams...")

if (!skip_generating_fit_mixTV_stream) {
  # Creating strings
  fitMixTVStrArr <- paste0(
    if (is.na(fitMixTVStrTempSplit[1])) {""} else {fitMixTVStrTempSplit[1]}, fitMixTVProbNameArr,
    fitMixTVStrTempSplit[2],
    if (Sys.info()["sysname"] == "Windows") {
      paste0("..\\..\\DATA\\", datPathArrBS)
    } else {
      paste0("../../DATA/", datPathArr)
    },
    fitMixTVStrTempSplit[3], dt_param[, TVCL], fitMixTVStrTempSplit[4], dt_param[, RCL],
    fitMixTVStrTempSplit[5], dt_param[, CVCL],
    fitMixTVStrTempSplit[6], dt_param[, TVVD], fitMixTVStrTempSplit[7], dt_param[, CVVD],
    fitMixTVStrTempSplit[8], dt_param[, TVKA], fitMixTVStrTempSplit[9], dt_param[, CVKA],
    fitMixTVStrTempSplit[10], dt_param[, MIXP], fitMixTVStrTempSplit[11], dt_param[, W_PROP],
    fitMixTVStrTempSplit[12], fitMixTVTxtNameArr,
    if (is.na(fitMixTVStrTempSplit[13])) {""} else {fitMixTVStrTempSplit[13]}
  )
  # Generating files
  dir.create(paste0(baseWd, "/FIT_MIXTV"), showWarnings = F, recursive = T)
  setwd(paste0(baseWd, "/FIT_MIXTV"))
  cat("\n\tChecking progress and deciding batch to run...")
  fitMixTVStrProg <- which(!file.exists(fitMixTVStrPathArr[toRun]))
  
  if (fitMixTVStrProg %>% length %>% is_greater_than(0)) {
    cat("\n")
    fitMixTVStr_start_Time <- Sys.time()
    j <- 0
    for (i in toRun[fitMixTVStrProg]) {
      j <- j + 1
      # Creating non-existent subdirectory
      thisDir = paste0(baseWd, "/FIT_MIXTV/", subDirArr[i])
      if (!dir.exists(thisDir)) {dir.create(thisDir, recursive = T)}
      # Writing files
      write(fitMixTVStrArr[i], fitMixTVStrPathArr[i]) %>% invisible
      fitMixTVStr_End_Time <- Sys.time()
      if (j%%100 == 0 || j == length(fitMixTVStrProg)) {
        cat(paste0("\r\t", strrep(" ", nchar(displayStr) - 2)))
        displayStr <- paste0("\r\tGenerating fit mixture (true value) control streams... ", j, "/", length(fitMixTVStrProg),
                             if (j == length(fitMixTVStrProg)) {"  |  Done!"} else {""},
                             "  |  Elapsed Time: ",
                             diffTime_to_display(fitMixTVStr_End_Time - fitMixTVStr_start_Time))
        cat(displayStr)
      }
    }
  } else {
    cat(" Done previously!")
  }
} else {
  cat(" Process skipped!")
}

cat("\n\n")
cat("********************************\n")
cat("******** Running NONMEM ********\n")
cat("********************************\n")

#+++++++++++++++++++++++++++++++++++++++++#
#### +++++Simulation with NONMEM+++++ #####
#+++++++++++++++++++++++++++++++++++++++++#

cat("\nSimulating PK datasets with NONMEM...")

if (!skip_sim) {
  setwd(paste0(baseWd, "/DATA"))
  nchar_sim_prefix <- nchar(sim_psn_prefix)
  nchar_sim_path_Arr <- nchar(simStrPathArr) + 1
  
  cat("\n\tChecking progress and deciding batch to run...")
  simStrProg <- which(!file.exists(simLstPathArr[toRun]))
  
  if (simStrProg %>% length %>% is_greater_than(0)) {
    simList <- toRun[simStrProg]
    cat("\n")
    first <- 1
    
    Sys.unsetenv("GFORTRAN_STDOUT_UNIT")
    PK_sim_start_Time <- Sys.time()
    while (first <= length(simList)) {
      last <- cumsum(nchar_sim_path_Arr[simList[first:length(simList)]]) %>% add(nchar_sim_prefix) %>% is_less_than(8000) %>%
        which %>% max %>% add(first) %>% subtract(1)
      batStr <- paste0(sim_psn_prefix, " ", do.call(paste, c(list(sep = " "), as.list(simStrPathArr[simList[first:last]]))))
      PK_sim_end_Time <- Sys.time()
      cat(paste0("\r\t", strrep(" ", nchar(displayStr) - 2)))
      displayStr <- paste0("\r\tSimulating PK datasets with NONMEM... batch ", first, "-", last, "/", length(simList), "  |  Elapsed time: ",
                           diffTime_to_display(PK_sim_end_Time - PK_sim_start_Time))
      cat(displayStr)
      if (Sys.info()["sysname"] == "Windows") {
        shell(batStr, intern = T)
      } else if (Sys.info()["sysname"] == "Linux") {
        system(batStr, intern = T)
      } else {
        stop("Operating system is not Windows nor Linux!")
      }
      first <- last + 1
    }
    PK_sim_end_Time <- Sys.time()
    cat(paste0("\r\tSimulating PK datasets with NONMEM... batch ", last, "/", length(simList),
               "  |  Done!  |  Elapsed time: ", diffTime_to_display(PK_sim_end_Time - PK_sim_start_Time)))
  } else {
    cat(" Done previously!")
  }
} else {
  cat(" Process skipped!")
}


#++++++++++++++++++++++++++++++++++++++++++#
#### +++++Fit mixture with NONMEM+++++ #####
#++++++++++++++++++++++++++++++++++++++++++#

cat("\nFitting mixture with NONMEM...")

if (!skip_fit_mix) {
  setwd(paste0(baseWd, "/FIT_MIX"))
  nchar_fitMix_prefix <- nchar(fitMix_psn_prefix)
  nchar_fitMix_path_Arr <- nchar(fitMixStrPathArr) + 1
  
  cat("\n\tChecking progress and deciding batch to run...")
  fitMixStrProg <- which(!file.exists(fitMixLstPathArr[toRun]))
  
  if (fitMixStrProg %>% length %>% is_greater_than(0)) {
    fitMixList <- toRun[fitMixStrProg]
    cat("\n")
    first <- 1
    
    Sys.unsetenv("GFORTRAN_STDOUT_UNIT")
    fitMix_start_Time <- Sys.time()
    while (first <= length(fitMixList)) {
      last <- cumsum(nchar_fitMix_path_Arr[fitMixList[first:length(fitMixList)]]) %>%
        add(nchar_fitMix_prefix) %>% is_less_than(8000) %>%
        which %>% max %>% add(first) %>% subtract(1)
      batStr <- paste0(fitMix_psn_prefix, " ",
                       do.call(paste, c(list(sep = " "), as.list(fitMixStrPathArr[fitMixList[first:last]]))))
      fitMix_end_Time <- Sys.time()
      cat(paste0("\r\t", strrep(" ", nchar(displayStr) - 2)))
      displayStr <- paste0("\r\tFitting mixture with NONMEM... batch ", first, "-", last, "/", length(fitMixList),
                           "  |  Elapsed time: ", diffTime_to_display(fitMix_end_Time - fitMix_start_Time))
      cat(displayStr)
      if (Sys.info()["sysname"] == "Windows") {
        shell(batStr, intern = T)
      } else if (Sys.info()["sysname"] == "Linux") {
        system(batStr, intern = T)
      } else {
        stop("Operating system is not Windows nor Linux!")
      }
      first <- last + 1
    }
    fitMix_end_Time <- Sys.time()
    cat(paste0("\r\tFitting mixture with NONMEM... batch ", last, "/", length(fitMixList),
               "  |  Done!  |  Elapsed time: ", diffTime_to_display(fitMix_end_Time - fitMix_start_Time)))
  } else {
    cat(" Done previously!")
  }
} else {
  cat(" Process skipped!")
}

#++++++++++++++++++++++++++++++++++++++++++++++++++#
#### +++++Fit without mixture with NONMEM+++++ #####
#++++++++++++++++++++++++++++++++++++++++++++++++++#

cat("\nFitting without mixture with NONMEM...")

if (!skip_fit_nomix) {
  setwd(paste0(baseWd, "/FIT_NOMIX"))
  nchar_fitNoMix_prefix <- nchar(fitNoMix_psn_prefix)
  nchar_fitNoMix_path_Arr <- nchar(fitNoMixStrPathArr) + 1
  
  cat("\n\tChecking progress and deciding batch to run...")
  fitNoMixStrProg <- which(!file.exists(fitNoMixLstPathArr[toRun]))
  
  if (fitNoMixStrProg %>% length %>% is_greater_than(0)) {
    fitNoMixList <- toRun[fitNoMixStrProg]
    cat("\n")
    first <- 1
    
    Sys.unsetenv("GFORTRAN_STDOUT_UNIT")
    fitNoMix_start_Time <- Sys.time()
    while (first <= length(fitNoMixList)) {
      last <- cumsum(nchar_fitNoMix_path_Arr[fitNoMixList[first:length(fitNoMixList)]]) %>%
        add(nchar_fitNoMix_prefix) %>% is_less_than(8000) %>%
        which %>% max %>% add(first) %>% subtract(1)
      batStr <- paste0(fitNoMix_psn_prefix, " ",
                       do.call(paste, c(list(sep = " "), as.list(fitNoMixStrPathArr[fitNoMixList[first:last]]))))
      fitNoMix_end_Time <- Sys.time()
      cat(paste0("\r\t", strrep(" ", nchar(displayStr) - 2)))
      displayStr <- paste0("\r\tFitting without mixture with NONMEM... batch ", first, "-", last, "/", length(fitNoMixList),
                           "  |  Elapsed time: ", diffTime_to_display(fitNoMix_end_Time - fitNoMix_start_Time))
      cat(displayStr)
      if (Sys.info()["sysname"] == "Windows") {
        shell(batStr, intern = T)
      } else if (Sys.info()["sysname"] == "Linux") {
        system(batStr, intern = T)
      } else {
        stop("Operating system is not Windows nor Linux!")
      }
      first <- last + 1
    }
    fitNoMix_end_Time <- Sys.time()
    cat(paste0("\r\tFitting without mixture with NONMEM... batch ", last, "/", length(fitNoMixList),
               "  |  Done!  |  Elapsed time: ", diffTime_to_display(fitNoMix_end_Time - fitNoMix_start_Time)))
  } else {
    cat(" Done previously!")
  }
} else {
  cat(" Process skipped!")
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#### +++++Fit mixture (true value) with NONMEM+++++ #####
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++#

cat("\nFitting mixture (true value) with NONMEM...")

if (!skip_fit_mixTV) {
  setwd(paste0(baseWd, "/FIT_MIXTV"))
  nchar_fitMixTV_prefix <- nchar(fitMixTV_psn_prefix)
  nchar_fitMixTV_path_Arr <- nchar(fitMixTVStrPathArr) + 1
  
  cat("\n\tChecking progress and deciding batch to run...")
  fitMixTVStrProg <- which(!file.exists(fitMixTVLstPathArr[toRun]))
  
  if (fitMixTVStrProg %>% length %>% is_greater_than(0)) {
    fitMixTVList <- toRun[fitMixTVStrProg]
    cat("\n")
    first <- 1
    
    Sys.unsetenv("GFORTRAN_STDOUT_UNIT")
    fitMixTV_start_Time <- Sys.time()
    while (first <= length(fitMixTVList)) {
      last <- cumsum(nchar_fitMixTV_path_Arr[fitMixTVList[first:length(fitMixTVList)]]) %>%
        add(nchar_fitMixTV_prefix) %>% is_less_than(8000) %>%
        which %>% max %>% add(first) %>% subtract(1)
      batStr <- paste0(fitMixTV_psn_prefix, " ",
                       do.call(paste, c(list(sep = " "), as.list(fitMixTVStrPathArr[fitMixTVList[first:last]]))))
      fitMixTV_end_Time <- Sys.time()
      cat(paste0("\r\t", strrep(" ", nchar(displayStr) - 2)))
      displayStr <- paste0("\r\tFitting mixture (true value) with NONMEM... batch ", first, "-", last, "/", length(fitMixTVList),
                           "  |  Elapsed time: ", diffTime_to_display(fitMixTV_end_Time - fitMixTV_start_Time))
      cat(displayStr)
      if (Sys.info()["sysname"] == "Windows") {
        shell(batStr, intern = T)
      } else if (Sys.info()["sysname"] == "Linux") {
        system(batStr, intern = T)
      } else {
        stop("Operating system is not Windows nor Linux!")
      }
      first <- last + 1
    }
    fitMixTV_end_Time <- Sys.time()
    cat(paste0("\r\tFitting mixture (true value) with NONMEM... batch ", last, "/", length(fitMixTVList),
               "  |  Done!  |  Elapsed time: ", diffTime_to_display(fitMixTV_end_Time - fitMixTV_start_Time)))
  } else {
    cat(" Done previously!")
  }
} else {
  cat(" Process skipped!")
}

cat("\n\n")
cat("************************************\n")
cat("******** Extracting results ********\n")
cat("************************************\n")

#++++++++++++++++++++++++++++++++++++#
#### +++++Extracting results+++++ ####
#++++++++++++++++++++++++++++++++++++#

if (!is.null(testRun)) {
  cat("\n!!!!---- Since non-full run has been specified, only the test-run content will be extracted! ----!!!!")
  cat("\n!!!!---- Extracted results will be written to a non-final file! ----!!!!\n")
}

#### Checking completeness of run ####

cat("\nChecking completeness of run...")
cat("\n\tChecking simulation progress... ")
  setwd(paste0(baseWd, "/DATA"))
  simStrCompChk <- which(!file.exists(simLstPathArr[toRun]))
  if (!length(simStrCompChk)) {cat(" Done!")} else {cat(" !!!--- Results incomplete ---!!!")}
cat("\n\tChecking fit mixture progress... ")
  setwd(paste0(baseWd, "/FIT_MIX"))
  fitMixStrCompChk <- which(!file.exists(fitMixLstPathArr[toRun]))
  if (!length(fitMixStrCompChk)) {cat(" Done!")} else {cat(" !!!--- Results incomplete ---!!!")}
cat("\n\tChecking fit without mixture progress... ")
  setwd(paste0(baseWd, "/FIT_NOMIX"))
  fitNoMixStrCompChk <- which(!file.exists(fitNoMixLstPathArr[toRun]))
  if (!length(fitNoMixStrCompChk)) {cat(" Done!")} else {cat(" !!!--- Results incomplete ---!!!")}
cat("\n\tChecking fit mixture (true value) progress... ")
  setwd(paste0(baseWd, "/FIT_MIXTV"))
  fitMixTVStrCompChk <- which(!file.exists(fitMixTVLstPathArr[toRun]))
  if (!length(fitMixTVStrCompChk)) {cat(" Done!")} else {cat(" !!!--- Results incomplete ---!!!")}
setwd(baseWd)
if (is.null(testRun) && file.exists(paste0(projName, "_DATASET_EXTRACT.csv")) &&
    file.exists(paste0(projName, "_INDIVIDUAL_EXTRACT.csv")) && file.exists("MISSING_TEMP.csv")) {
  exist_tab_miss <- T
  tab_miss <- fread("MISSING_TEMP.csv", header = T)
  toExtract <- tab_miss[missingSimLst == 1 | missingFitMixLst == 1 | missingFitNoMixLst == 1 | missingFitMixTVLst == 1, datasetID]
} else {
  exist_tab_miss <- F
  toExtract <- toRun
}
if (any(length(simStrCompChk), length(fitMixStrCompChk), length(fitNoMixStrCompChk), length(fitMixTVStrCompChk))) {
  cat("\n\n\t !!!! Missing lst files identified! Missing files are written to \"~/MISSING_TEMP.csv\"!!!!")
  cat("\n\t !!!! This does not necessarily mean missing run. Unrectifiable runs can also have missing lst files. !!!!")
  cat("\n\t !!!! Extraction will continue... !!!!\n")
}
setwd(baseWd)
write_csv(data.table(
  datasetID = toRun,
  missingSimLst = toRun %in% toRun[simStrCompChk] %>% as.numeric(),
  missingFitMixLst = toRun %in% toRun[fitMixStrCompChk] %>% as.numeric(),
  missingFitNoMixLst = toRun %in% toRun[fitNoMixStrCompChk] %>% as.numeric(),
  missingFitMixTVLst = toRun %in% toRun[fitMixTVStrCompChk] %>% as.numeric()
), if (is.null(testRun)) {"MISSING_TEMP.csv"} else {"MISSING_TEMP_TEST.csv"})

#### Dataset information Extraction ####

cat("\nExtracting dataset level information...")

# > Extracting from simulated dataset ####
datMissingArr <- numeric()
setwd(paste0(baseWd, "/DATA"))
cat("\n")
j <- 0
extractDat_start_Time <- Sys.time()
arr_DAT_MIXP <- sapply(toExtract, function(i) {
  j <<- j + 1
  tryCatch({
    datMixp = fread(datPathArr[i], header = T, skip = 1)[, .SD[1], ID][, 2 - sum(SIMMIX)/.N]
  }, error = function(e) {
    assign("datMixp", NA, envir = parent.frame(4))
    datMissingArr <- c(datMissingArr, toExtract)
  })
  extractDat_End_Time <- Sys.time()
  if (j%%100 == 0 || j == length(toExtract)) {
    strDisplay = paste0("\r\tExtracting from datasets... ", j, "/", length(toExtract),
                        if (j == length(toExtract)) {"... Finishing up..."} else {""},
                        "  |  Elapsed Time: ",
                        diffTime_to_display(extractDat_End_Time - extractDat_start_Time))
    if (j == length(toExtract)) {maxCatLen <<- nchar(strDisplay) - 2}
    cat(strDisplay)
  }
  return(datMixp)
}) %>% data.table(DAT_MIXP = .)
extractDat_End_Time <- Sys.time()
cat(paste0("\r\t", strrep(" ", maxCatLen)))
strDisplay <- paste0("\r\tExtracting from datasets... ", j, "/", length(toExtract),
                     if (j == length(toExtract)) {"... Finished!"} else {""},
                     "  |  Elapsed Time: ",
                     diffTime_to_display(extractDat_End_Time - extractDat_start_Time))
cat(strDisplay)

# > Extracting from Mix.ext ####
mixExtMissingArr <- numeric()
setwd(paste0(baseWd, "/FIT_MIX"))
cat("\n")
j <- 0
extractMixExt_start_Time <- Sys.time()
tab_Mix_ext <- sapply(toExtract, function(i) {
  j <<- j + 1
  tryCatch({
    tryCatch({tabMixExt = fread(fitMixExtPathArr[i], header = T, skip = 1)},
             error = function(e) {stop("tabErr")})
    tryCatch({FE_Line = tabMixExt[ITERATION == -1000000000]; if (nrow(FE_Line) == 0) stop()},
             error = function(e) {stop("FE_Err")})
    tryCatch({RSE_Line = tabMixExt[ITERATION == -1000000001]; if (nrow(RSE_Line) == 0) stop()},
             error = function(e) {stop("RSE_Err")})
    FE = FE_Line[, c(THETA1, THETA2, THETA3, THETA4, THETA5, THETA6, THETA7, THETA8, THETA9, OBJ)]
    RSE = RSE_Line[, c(THETA1, THETA2, THETA3, THETA4, THETA5, THETA6, THETA7, THETA8, THETA9)] %>%
      divide_by(FE[1:9])
  }, error = function(e) {
    if (e$message == "tabErr" || e$message == "FE_Err") {
      assign("FE", rep(NA, 10), envir = parent.frame(4))
      assign("RSE", rep(NA, 9), envir = parent.frame(4))
      assign("mixExtMissingArr", c(mixExtMissingArr, i), envir = .GlobalEnv)
    } else if (e$message == "RSE_Err") {
      assign("FE", FE_Line[, c(THETA1, THETA2, THETA3, THETA4, THETA5, THETA6, THETA7, THETA8, THETA9, OBJ)],
             envir = parent.frame(4))
      assign("RSE", rep(NA, 9), envir = parent.frame(4))
    }
  })
  extractMixExt_End_Time <- Sys.time()
  if (j%%100 == 0 || j == length(toExtract)) {
    strDisplay = paste0("\r\tExtracting from Mix.ext... ", j, "/", length(toExtract),
                        if (j == length(toExtract)) {"... Finishing up..."} else {""},
                        "  |  Elapsed Time: ",
                        diffTime_to_display(extractMixExt_End_Time - extractMixExt_start_Time))
    if (j == length(toExtract)) {maxCatLen <<- nchar(strDisplay) - 2}
    cat(strDisplay)
  }
  return(c(FE, RSE))
}) %>% t %>% as.data.table %>% set_colnames(., c(
  "MIX_FE_TVCL", "MIX_FE_RCL", "MIX_FE_CVCL", "MIX_FE_TVVD", "MIX_FE_CVVD",
  "MIX_FE_TVKA", "MIX_FE_CVKA", "MIX_FE_LOGIT_MIXP", "MIX_FE_W_PROP", "MIX_OFV",
  "MIX_RSE_TVCL", "MIX_RSE_RCL", "MIX_RSE_CVCL", "MIX_RSE_TVVD", "MIX_RSE_CVVD",
  "MIX_RSE_TVKA", "MIX_RSE_CVKA", "MIX_RSE_LOGIT_MIXP", "MIX_RSE_W_PROP"
))
extractMixExt_End_Time <- Sys.time()
cat(paste0("\r\t", strrep(" ", maxCatLen)))
strDisplay <- paste0("\r\tExtracting from Mix.ext... ", j, "/", length(toExtract),
                     if (j == length(toExtract)) {"... Finished!"} else {""},
                     "  |  Elapsed Time: ",
                     diffTime_to_display(extractMixExt_End_Time - extractMixExt_start_Time))
cat(strDisplay)

# > Extracting from Mix.lst ####
mixLstMissingArr <- numeric()
setwd(paste0(baseWd, "/FIT_MIX"))
cat("\n")
j <- 0
extractMixLst_start_Time <- Sys.time()
str_Mix_lst <- sapply(toExtract, function(i) {
  j <<- j + 1
  tryCatch({
    tryCatch({strMixLst = read_file(fitMixLstPathArr[i])},
             error = function(e) {stop("strErr")})
    STATUS = c(
      regexpr("MINIMIZATION SUCCESSFUL", strMixLst) %>% is_greater_than(0) %>% any,
      regexpr("COVARIANCE MATRIX OF ESTIMATE", strMixLst) %>% is_greater_than(0) %>% any,
      regexpr("NEAR ITS BOUNDARY", strMixLst) %>% is_greater_than(0) %>% any,
      regexpr("COVARIANCE STEP ABORTED", strMixLst) %>% is_greater_than(0) %>% any,
      regexpr("DUE TO ROUNDING ERROR", strMixLst) %>% is_greater_than(0) %>% any,
      regexpr("PROBLEMS OCCURRED WITH THE MINIMIZATION", strMixLst) %>% is_greater_than(0) %>% any,
      regexpr("TERMINAT", strMixLst) %>% is_greater_than(0) %>% any
    )
  }, error = function(e) {
    assign("STATUS", rep(NA, 7), envir = parent.frame(4))
    assign("mixLstMissingArr", c(mixLstMissingArr, i), envir = .GlobalEnv)
  })
  extractMixLst_End_Time <- Sys.time()
  if (j%%100 == 0 || j == length(toExtract)) {
    strDisplay = paste0("\r\tExtracting from Mix.lst... ", j, "/", length(toExtract),
                        if (j == length(toExtract)) {"... Finishing up..."} else {""},
                        "  |  Elapsed Time: ",
                        diffTime_to_display(extractMixLst_End_Time - extractMixLst_start_Time))
    if (j == length(toExtract)) {maxCatLen <<- nchar(strDisplay) - 2}
    cat(strDisplay)
  }
  return(STATUS)
}) %>% t %>% as.data.table %>% set_colnames(., c(
  "MIX_MIN_SUCCESS", "MIX_COV_SUCCESS", "MIX_BOUNDARY", "MIX_COV_WARN",
  "MIX_ROUNDING_ERR", "MIX_MIN_PROB", "MIX_TERMINATION"
))
extractMixLst_End_Time <- Sys.time()
cat(paste0("\r\t", strrep(" ", maxCatLen)))
strDisplay <- paste0("\r\tExtracting from Mix.lst... ", j, "/", length(toExtract),
                     if (j == length(toExtract)) {"... Finished!"} else {""},
                     "  |  Elapsed Time: ",
                     diffTime_to_display(extractMixLst_End_Time - extractMixLst_start_Time))
cat(strDisplay)

# > Extracting from Mix.shk ####
mixShkMissingArr <- numeric()
setwd(paste0(baseWd, "/FIT_MIX"))
cat("\n")
j <- 0
extractMixShk_start_Time <- Sys.time()
tab_Mix_shk <- sapply(toExtract, function(i) {
  j <<- j + 1
  tryCatch({
    tryCatch({tabMixShk = fread(fitMixShkPathArr[i], header = T, skip = 1)},
             error = function(e) {stop("tabErr")})
    tryCatch({SHR1_ETA = tabMixShk[TYPE == 4 & SUBPOP == 1]; if (nrow(SHR1_ETA) == 0) stop()},
             error = function(e) {assign("SHR1_ETA", NULL, envir = parent.frame(4))})
    tryCatch({SHR2_ETA = tabMixShk[TYPE == 4 & SUBPOP == 2]; if (nrow(SHR2_ETA) == 0) stop()},
             error = function(e) {assign("SHR2_ETA", NULL, envir = parent.frame(4))})
    tryCatch({SHR1_EPS = tabMixShk[TYPE == 5 & SUBPOP == 1]; if (nrow(SHR1_EPS) == 0) stop()},
             error = function(e) {assign("SHR1_EPS", NULL, envir = parent.frame(4))})
    tryCatch({SHR2_EPS = tabMixShk[TYPE == 5 & SUBPOP == 2]; if (nrow(SHR2_EPS) == 0) stop()},
             error = function(e) {assign("SHR2_EPS", NULL, envir = parent.frame(4))})
    if (!is.null(SHR1_ETA)) {SHR1_ETA = SHR1_ETA[1, c(`ETA(1)`, `ETA(2)`, `ETA(3)`)/100]} else {SHR1_ETA = rep(NA, 3)}
    if (!is.null(SHR2_ETA)) {SHR2_ETA = SHR2_ETA[1, c(`ETA(1)`, `ETA(2)`, `ETA(3)`)/100]} else {SHR2_ETA = rep(NA, 3)}
    if (!is.null(SHR1_EPS)) {SHR1_EPS = SHR1_EPS[1, `ETA(1)`/100]} else {SHR1_EPS = NA}
    if (!is.null(SHR2_EPS)) {SHR2_EPS = SHR2_EPS[1, `ETA(1)`/100]} else {SHR2_EPS = NA}
  }, error = function(e) {
    assign("SHR1_ETA", rep(NA, 3), envir = parent.frame(4))
    assign("SHR2_ETA", rep(NA, 3), envir = parent.frame(4))
    assign("SHR1_EPS", NA, envir = parent.frame(4))
    assign("SHR2_EPS", NA, envir = parent.frame(4))
    if (e$message == "tabErr") {assign("mixShkMissingArr", c(mixShkMissingArr, i), envir = .GlobalEnv)}
  })
  extractMixShk_End_Time <- Sys.time()
  if (j%%100 == 0 || j == length(toExtract)) {
    strDisplay = paste0("\r\tExtracting from Mix.shk... ", j, "/", length(toExtract),
                        if (j == length(toExtract)) {"... Finishing up..."} else {""},
                        "  |  Elapsed Time: ",
                        diffTime_to_display(extractMixShk_End_Time - extractMixShk_start_Time))
    if (j == length(toExtract)) {maxCatLen <<- nchar(strDisplay) - 2}
    cat(strDisplay)
  }
  return(c(SHR1_ETA, SHR1_EPS, SHR2_ETA, SHR2_EPS))
}) %>% t %>% as.data.table %>% set_colnames(., c(
  "MIX_SHR1_CL", "MIX_SHR1_VD", "MIX_SHR1_KA", "MIX_SHR1_EPS_PROP",
  "MIX_SHR2_CL", "MIX_SHR2_VD", "MIX_SHR2_KA", "MIX_SHR2_EPS_PROP"
))
extractMixShk_End_Time <- Sys.time()
cat(paste0("\r\t", strrep(" ", maxCatLen)))
strDisplay <- paste0("\r\tExtracting from Mix.shk... ", j, "/", length(toExtract),
                     if (j == length(toExtract)) {"... Finished!"} else {""},
                     "  |  Elapsed Time: ",
                     diffTime_to_display(extractMixShk_End_Time - extractMixShk_start_Time))
cat(strDisplay)

# > Extracting from MixTV.ext ####
mixTVExtMissingArr <- numeric()
setwd(paste0(baseWd, "/FIT_MIXTV"))
cat("\n")
j <- 0
extractMixTVExt_start_Time <- Sys.time()
tab_MixTV_ext <- sapply(toExtract, function(i) {
  j <<- j + 1
  tryCatch({
    tryCatch({tabMixTVExt = fread(fitMixTVExtPathArr[i], header = T, skip = 1)},
             error = function(e) {stop("tabErr")})
    tryCatch({FE_Line = tabMixTVExt[ITERATION == -1000000000]; if (nrow(FE_Line) == 0) stop()},
             error = function(e) {stop("FE_Err")})
    MixTV_OBJ = FE_Line[, OBJ]
  }, error = function(e) {
    if (e$message == "tabErr") {assign("mixTVExtMissingArr", c(mixTVExtMissingArr, i), envir = .GlobalEnv)}
    assign("MixTV_OBJ", NA, envir = parent.frame(4))
  })
  extractMixTVExt_End_Time <- Sys.time()
  if (j%%100 == 0 || j == length(toExtract)) {
    strDisplay = paste0("\r\tExtracting from MixTV.ext... ", j, "/", length(toExtract),
                        if (j == length(toExtract)) {"... Finishing up..."} else {""},
                        "  |  Elapsed Time: ",
                        diffTime_to_display(extractMixTVExt_End_Time - extractMixTVExt_start_Time))
    if (j == length(toExtract)) {maxCatLen <<- nchar(strDisplay) - 2}
    cat(strDisplay)
  }
  return(MixTV_OBJ)
}) %>% data.table(MIXTV_OFV = .)
extractMixTVExt_End_Time <- Sys.time()
cat(paste0("\r\t", strrep(" ", maxCatLen)))
strDisplay <- paste0("\r\tExtracting from MixTV.ext... ", j, "/", length(toExtract),
                     if (j == length(toExtract)) {"... Finished!"} else {""},
                     "  |  Elapsed Time: ",
                     diffTime_to_display(extractMixTVExt_End_Time - extractMixTVExt_start_Time))
cat(strDisplay)

# > Extracting from NoMix.ext ####
noMixExtMissingArr <- numeric()
setwd(paste0(baseWd, "/FIT_NOMIX"))
cat("\n")
j <- 0
extractNoMixExt_start_Time <- Sys.time()
tab_NoMix_ext <- sapply(toExtract, function(i) {
  j <<- j + 1
  tryCatch({
    tryCatch({tabNoMixExt = fread(fitNoMixExtPathArr[i], header = T, skip = 1)},
             error = function(e) {stop("tabErr")})
    tryCatch({FE_Line = tabNoMixExt[ITERATION == -1000000000]; if (nrow(FE_Line) == 0) stop()},
             error = function(e) {stop("FE_Err")})
    tryCatch({RSE_Line = tabNoMixExt[ITERATION == -1000000001]; if (nrow(RSE_Line) == 0) stop()},
             error = function(e) {stop("RSE_Err")})
    FE = FE_Line[, c(THETA1, THETA3, THETA4, THETA5, THETA6, THETA7, THETA9, OBJ)]
    RSE = RSE_Line[, c(THETA1, THETA3, THETA4, THETA5, THETA6, THETA7, THETA9)] %>% divide_by(FE[1:7])
  }, error = function(e) {
    if (e$message == "tabErr" || e$message == "FE_Err") {
      assign("FE", rep(NA, 8), envir = parent.frame(4))
      assign("RSE", rep(NA, 7), envir = parent.frame(4))
      assign("noMixExtMissingArr", c(noMixExtMissingArr, i), envir = .GlobalEnv)
    } else if (e$message == "RSE_Err") {
      assign("FE", FE_Line[, c(THETA1, THETA3, THETA4, THETA5, THETA6, THETA7, THETA9, OBJ)],
             envir = parent.frame(4))
      assign("RSE", rep(NA, 7), envir = parent.frame(4))
    }
  })
  extractNoMixExt_End_Time <- Sys.time()
  if (j%%100 == 0 || j == length(toExtract)) {
    strDisplay = paste0("\r\tExtracting from NoMix.ext... ", j, "/", length(toExtract),
                        if (j == length(toExtract)) {"... Finishing up..."} else {""},
                        "  |  Elapsed Time: ",
                        diffTime_to_display(extractNoMixExt_End_Time - extractNoMixExt_start_Time))
    if (j == length(toExtract)) {maxCatLen <<- nchar(strDisplay) - 2}
    cat(strDisplay)
  }
  return(c(FE, RSE))
}) %>% t %>% as.data.table %>% set_colnames(., c(
  "NOMIX_FE_TVCL", "NOMIX_FE_CVCL", "NOMIX_FE_TVVD", "NOMIX_FE_CVVD",
  "NOMIX_FE_TVKA", "NOMIX_FE_CVKA", "NOMIX_FE_W_PROP", "NOMIX_OFV",
  "NOMIX_RSE_TVCL", "NOMIX_RSE_CVCL", "NOMIX_RSE_TVVD", "NOMIX_RSE_CVVD",
  "NOMIX_RSE_TVKA", "NOMIX_RSE_CVKA", "NOMIX_RSE_W_PROP"
))
extractNoMixExt_End_Time <- Sys.time()
cat(paste0("\r\t", strrep(" ", maxCatLen)))
strDisplay <- paste0("\r\tExtracting from NoMix.ext... ", j, "/", length(toExtract),
                     if (j == length(toExtract)) {"... Finished!"} else {""},
                     "  |  Elapsed Time: ",
                     diffTime_to_display(extractNoMixExt_End_Time - extractNoMixExt_start_Time))
cat(strDisplay)

# > Extracting from NoMix.lst ####
noMixLstMissingArr <- numeric()
setwd(paste0(baseWd, "/FIT_NOMIX"))
cat("\n")
j <- 0
extractNoMixLst_start_Time <- Sys.time()
str_NoMix_lst <- sapply(toExtract, function(i) {
  j <<- j + 1
  tryCatch({
    tryCatch({strNoMixLst = read_file(fitNoMixLstPathArr[i])},
             error = function(e) {stop("strErr")})
    STATUS = c(
      regexpr("MINIMIZATION SUCCESSFUL", strNoMixLst) %>% is_greater_than(0) %>% any,
      regexpr("COVARIANCE MATRIX OF ESTIMATE", strNoMixLst) %>% is_greater_than(0) %>% any,
      regexpr("NEAR ITS BOUNDARY", strNoMixLst) %>% is_greater_than(0) %>% any,
      regexpr("COVARIANCE STEP ABORTED", strNoMixLst) %>% is_greater_than(0) %>% any,
      regexpr("DUE TO ROUNDING ERROR", strNoMixLst) %>% is_greater_than(0) %>% any,
      regexpr("PROBLEMS OCCURRED WITH THE MINIMIZATION", strNoMixLst) %>% is_greater_than(0) %>% any,
      regexpr("TERMINAT", strNoMixLst) %>% is_greater_than(0) %>% any
    )
  }, error = function(e) {
    assign("STATUS", rep(NA, 7), envir = parent.frame(4))
    assign("noMixLstMissingArr", c(noMixLstMissingArr, i), envir = .GlobalEnv)
  })
  extractNoMixLst_End_Time <- Sys.time()
  if (j%%100 == 0 || j == length(toExtract)) {
    strDisplay = paste0("\r\tExtracting from NoMix.lst... ", j, "/", length(toExtract),
                        if (j == length(toExtract)) {"... Finishing up..."} else {""},
                        "  |  Elapsed Time: ",
                        diffTime_to_display(extractNoMixLst_End_Time - extractNoMixLst_start_Time))
    if (j == length(toExtract)) {maxCatLen <<- nchar(strDisplay) - 2}
    cat(strDisplay)
  }
  return(STATUS)
}) %>% t %>% as.data.table %>% set_colnames(., c(
  "NOMIX_MIN_SUCCESS", "NOMIX_COV_SUCCESS", "NOMIX_BOUNDARY", "NOMIX_COV_WARN",
  "NOMIX_ROUNDING_ERR", "NOMIX_MIN_PROB", "NOMIX_TERMINATION"
))
extractNoMixLst_End_Time <- Sys.time()
cat(paste0("\r\t", strrep(" ", maxCatLen)))
strDisplay <- paste0("\r\tExtracting from NoMix.lst... ", j, "/", length(toExtract),
                     if (j == length(toExtract)) {"... Finished!"} else {""},
                     "  |  Elapsed Time: ",
                     diffTime_to_display(extractNoMixLst_End_Time - extractNoMixLst_start_Time))
cat(strDisplay)

# > Extracting from NoMix.shk ####
noMixShkMissingArr <- numeric()
setwd(paste0(baseWd, "/FIT_NOMIX"))
cat("\n")
j <- 0
extractNoMixShk_start_Time <- Sys.time()
tab_NoMix_shk <- sapply(toExtract, function(i) {
  j <<- j + 1
  tryCatch({
    tryCatch({tabNoMixShk = fread(fitNoMixShkPathArr[i], header = T, skip = 1)},
             error = function(e) {stop("tabErr")})
    tryCatch({SHR_ETA = tabNoMixShk[TYPE == 4]; if (nrow(SHR_ETA) == 0) stop()},
             error = function(e) {assign("SHR_ETA", NULL, envir = parent.frame(4))})
    tryCatch({SHR_EPS = tabNoMixShk[TYPE == 5]; if (nrow(SHR_EPS) == 0) stop()},
             error = function(e) {assign("SHR_EPS", NULL, envir = parent.frame(4))})
    if (!is.null(SHR_ETA)) {SHR_ETA = SHR_ETA[1, c(`ETA(1)`, `ETA(2)`, `ETA(3)`)/100]} else {SHR_ETA = rep(NA, 3)}
    if (!is.null(SHR_EPS)) {SHR_EPS = SHR_EPS[1, `ETA(1)`/100]} else {SHR_EPS = NA}
  }, error = function(e) {
    assign("SHR_ETA", rep(NA, 3), envir = parent.frame(4))
    assign("SHR_EPS", NA, envir = parent.frame(4))
    if (e$message == "tabErr") {assign("noMixShkMissingArr", c(noMixShkMissingArr, i), envir = .GlobalEnv)}
  })
  extractNoMixShk_End_Time <- Sys.time()
  if (j%%100 == 0 || j == length(toExtract)) {
    strDisplay = paste0("\r\tExtracting from NoMix.shk... ", j, "/", length(toExtract),
                        if (j == length(toExtract)) {"... Finishing up..."} else {""},
                        "  |  Elapsed Time: ",
                        diffTime_to_display(extractNoMixShk_End_Time - extractNoMixShk_start_Time))
    if (j == length(toExtract)) {maxCatLen <<- nchar(strDisplay) - 2}
    cat(strDisplay)
  }
  return(c(SHR_ETA, SHR_EPS))
}) %>% t %>% as.data.table %>% set_colnames(., c(
  "NOMIX_SHR1_CL", "NOMIX_SHR1_VD", "NOMIX_SHR1_KA", "NOMIX_SHR1_EPS_PROP"
))
extractNoMixShk_End_Time <- Sys.time()
cat(paste0("\r\t", strrep(" ", maxCatLen)))
strDisplay <- paste0("\r\tExtracting from NoMix.shk... ", j, "/", length(toExtract),
                     if (j == length(toExtract)) {"... Finished!"} else {""},
                     "  |  Elapsed Time: ",
                     diffTime_to_display(extractNoMixShk_End_Time - extractNoMixShk_start_Time))
cat(strDisplay)

# >>> Writing dataset extract file ####
cat("\n\t> Writing dataset extract file ...")
writing_tab_dat_start_time <- Sys.time()
setwd(baseWd)
tab_dat <- cbind(
  data.table(datasetID = toExtract, DOSE = rep(dose, length(toExtract))), 
  dt_param[toExtract, .(SAMPLESIZE = SampleSize, SIM_TVCL = TVCL, SIM_RCL = RCL, SIM_CVCL = CVCL,
                        SIM_TVVD = TVVD, SIM_CVVD = CVVD, SIM_TVKA = TVKA, SIM_CVKA = CVKA,
                        SIM_MIXP = MIXP, SIM_W_PROP = W_PROP)],
  arr_DAT_MIXP,
  str_Mix_lst,
  tab_Mix_ext,
  tab_Mix_shk,
  str_NoMix_lst,
  tab_NoMix_ext,
  tab_NoMix_shk,
  tab_MixTV_ext
)
if (is.null(testRun)) {
  dat_extract_file_name <- paste0(projName, "_DATASET_EXTRACT.csv")
  if (file.exists(dat_extract_file_name) && file.exists(paste0(projName, "_INDIVIDUAL_EXTRACT.csv")) && file.exists("MISSING_TEMP.csv")) {
    cur_tab_dat <- fread(dat_extract_file_name, header = T)
    cur_tab_dat[match(tab_dat[, datasetID], cur_tab_dat[, datasetID])] <- tab_dat
    tab_dat <- cur_tab_dat
  }
} else {
  ExistingDatTest <- sapply(dir(), function(dr) {
    num = if (regexpr(paste0(projName, "_DATASET_EXTRACT_TEST"), dr) %>% as.numeric %>% equals(-1) %>% not) {
      substr(dr, nchar(dr) - 5, nchar(dr) - 4)
    } else {NULL}
    return(num)
  }, USE.NAMES = F) %>% unlist %>% as.numeric
  if (length(ExistingDatTest)) {
    maxExistingDatTest <- max(ExistingDatTest)
    if (maxExistingDatTest < 99) {
      DatTestSerial <- maxExistingDatTest + 1
    } else {
      DatTestSerial <- 99
    }
  } else {
    DatTestSerial <- 1
  }
  dat_extract_file_name <- paste0(projName, "_DATASET_EXTRACT_TEST", sprintf("%02d", DatTestSerial), ".csv")
}
write_csv(tab_dat, dat_extract_file_name)
writing_tab_dat_End_time <- Sys.time()
cat(paste0(" Done !  |  Elapsed Time: ", diffTime_to_display(writing_tab_dat_End_time - writing_tab_dat_start_time)))

#### Individual information Extraction ####

cat("\nExtracting individual level information...")

# > Extracting from Mix.phi ####
mixPhiMissingArr <- numeric()
setwd(paste0(baseWd, "/FIT_MIX"))
cat("\n")
j <- 0
extractMixPhi_start_Time <- Sys.time()
tab_Mix_phi <- lapply(toExtract, function(i) {
  j <<- j + 1
  tryCatch({
    thisN = dt_param[i, SampleSize]
    tryCatch({
      tabMixPhi = fread(fitMixPhiPathArr[i], header = T, skip = 1)
      if (nrow(tabMixPhi) != thisN) {stop()}
    }, error = function(e) {stop("tabErr")})
    OBJ_arr = tabMixPhi[, OBJ]
  }, error = function(e) {
    if (e$message == "tabErr") {
      assign("OBJ_arr", rep(NA, thisN), envir = parent.frame(4))
      assign("mixPhiMissingArr", c(mixPhiMissingArr, i), envir = .GlobalEnv)
    }
  })
  extractMixPhi_End_Time <- Sys.time()
  if (j%%100 == 0 || j == length(toExtract)) {
    strDisplay = paste0("\r\tExtracting from Mix.phi... ", j, "/", length(toExtract),
                        if (j == length(toExtract)) {"... Finishing up..."} else {""},
                        "  |  Elapsed Time: ",
                        diffTime_to_display(extractMixPhi_End_Time - extractMixPhi_start_Time))
    if (j == length(toExtract)) {maxCatLen <<- nchar(strDisplay) - 2}
    cat(strDisplay)
  }
  return(OBJ_arr)
}) %>% do.call(c, .) %>% data.table(MIX_PHI_OFV = .)
extractMixPhi_End_Time <- Sys.time()
cat(paste0("\r\t", strrep(" ", maxCatLen)))
strDisplay <- paste0("\r\tExtracting from Mix.phi... ", j, "/", length(toExtract),
                     if (j == length(toExtract)) {"... Finished!"} else {""},
                     "  |  Elapsed Time: ",
                     diffTime_to_display(extractMixPhi_End_Time - extractMixPhi_start_Time))
cat(strDisplay)

# > Extracting from Mix.phm ####
mixPhmMissingArr <- numeric()
setwd(paste0(baseWd, "/FIT_MIX"))
cat("\n")
j <- 0
extractMixPhm_start_Time <- Sys.time()
tab_Mix_phm <- lapply(toExtract, function(i) {
  j <<- j + 1
  tryCatch({
    thisN = dt_param[i, SampleSize]
    tryCatch({
      tabMixPhm = fread(fitMixPhmPathArr[i], header = T, skip = 1)
      if (nrow(tabMixPhm) != thisN * 2) {stop()}
    }, error = function(e) {stop("tabErr")})
    tab1 = tabMixPhm[SUBPOP == 1, .(MIX_PMIX1 = PMIX, MIX_PHM_OFV1 = OBJ)]
    tab2 = tabMixPhm[SUBPOP == 2, .(MIX_PMIX2 = PMIX, MIX_PHM_OFV2 = OBJ)]
    tab = cbind(tab1, tab2)
  }, error = function(e) {
    if (e$message == "tabErr") {
      tab = data.table(MIX_PMIX1 = rep(NA, thisN), MIX_PHM_OFV1 = rep(NA, thisN),
                       MIX_PMIX2 = rep(NA, thisN), MIX_PHM_OFV2 = rep(NA, thisN))
      assign("tab", tab, envir = parent.frame(4))
      assign("mixPhmMissingArr", c(mixPhmMissingArr, i), envir = .GlobalEnv)
    }
  })
  extractMixPhm_End_Time <- Sys.time()
  if (j%%100 == 0 || j == length(toExtract)) {
    strDisplay = paste0("\r\tExtracting from Mix.phm... ", j, "/", length(toExtract),
                        if (j == length(toExtract)) {"... Finishing up..."} else {""},
                        "  |  Elapsed Time: ",
                        diffTime_to_display(extractMixPhm_End_Time - extractMixPhm_start_Time))
    if (j == length(toExtract)) {maxCatLen <<- nchar(strDisplay) - 2}
    cat(strDisplay)
  }
  return(tab)
}) %>% rbindlist
extractMixPhm_End_Time <- Sys.time()
cat(paste0("\r\t", strrep(" ", maxCatLen)))
strDisplay <- paste0("\r\tExtracting from Mix.phm... ", j, "/", length(toExtract),
                     if (j == length(toExtract)) {"... Finished!"} else {""},
                     "  |  Elapsed Time: ",
                     diffTime_to_display(extractMixPhm_End_Time - extractMixPhm_start_Time))
cat(strDisplay)

# > Extracting from Mix.txt ####
mixTxtMissingArr <- numeric()
setwd(paste0(baseWd, "/FIT_MIX"))
cat("\n")
j <- 0
extractMixTxt_start_Time <- Sys.time()
tab_Mix_txt <- lapply(toExtract, function(i) {
  j <<- j + 1
  tryCatch({
    thisN = dt_param[i, SampleSize]
    tryCatch({
      tabMixTxt = fread(fitMixTxtPathArr[i], header = T, skip = 1)
      tabID = tabMixTxt[, .SD[1], ID]
      if (tabID[, ID %>% equals(seq(thisN)) %>% all %>% not]) {stop()}
    }, error = function(e) {stop("tabErr")})
    tab = tabID[, .(
      ID_in_set = ID, SIMMIX, SIM_CL, SIM_VD, SIM_KA,
      MIX_CL = CL, MIX_VD = VD, MIX_KA = KA
    )]
  }, error = function(e) {
    if (e$message == "tabErr") {
      tab = data.table(
        ID_in_set = rep(NA, thisN), SIMMIX = rep(NA, thisN),
        SIM_CL = rep(NA, thisN), SIM_VD = rep(NA, thisN), SIM_KA = rep(NA, thisN),
        MIX_CL = rep(NA, thisN), MIX_VD = rep(NA, thisN), MIX_KA = rep(NA, thisN)
      )
      assign("tab", tab, envir = parent.frame(4))
      assign("mixTxtMissingArr", c(mixTxtMissingArr, i), envir = .GlobalEnv)
    }
  })
  extractMixTxt_End_Time <- Sys.time()
  if (j%%100 == 0 || j == length(toExtract)) {
    strDisplay = paste0("\r\tExtracting from Mix.txt... ", j, "/", length(toExtract),
                        if (j == length(toExtract)) {"... Finishing up..."} else {""},
                        "  |  Elapsed Time: ",
                        diffTime_to_display(extractMixTxt_End_Time - extractMixTxt_start_Time))
    if (j == length(toExtract)) {maxCatLen <<- nchar(strDisplay) - 2}
    cat(strDisplay)
  }
  return(tab)
}) %>% rbindlist
extractMixTxt_End_Time <- Sys.time()
cat(paste0("\r\t", strrep(" ", maxCatLen)))
strDisplay <- paste0("\r\tExtracting from Mix.txt... ", j, "/", length(toExtract),
                     if (j == length(toExtract)) {"... Finished!"} else {""},
                     "  |  Elapsed Time: ",
                     diffTime_to_display(extractMixTxt_End_Time - extractMixTxt_start_Time))
cat(strDisplay)

# > Extracting from MixTV.phi ####
mixTVPhiMissingArr <- numeric()
setwd(paste0(baseWd, "/FIT_MIXTV"))
cat("\n")
j <- 0
extractMixTVPhi_start_Time <- Sys.time()
tab_MixTV_phi <- lapply(toExtract, function(i) {
  j <<- j + 1
  tryCatch({
    thisN = dt_param[i, SampleSize]
    tryCatch({
      tabMixTVPhi = fread(fitMixTVPhiPathArr[i], header = T, skip = 1)
      if (nrow(tabMixTVPhi) != thisN) {stop()}
    }, error = function(e) {stop("tabErr")})
    OBJ_arr = tabMixTVPhi[, OBJ]
  }, error = function(e) {
    if (e$message == "tabErr") {
      assign("OBJ_arr", rep(NA, thisN), envir = parent.frame(4))
      assign("mixTVPhiMissingArr", c(mixTVPhiMissingArr, i), envir = .GlobalEnv)
    }
  })
  extractMixTVPhi_End_Time <- Sys.time()
  if (j%%100 == 0 || j == length(toExtract)) {
    strDisplay = paste0("\r\tExtracting from MixTV.phi... ", j, "/", length(toExtract),
                        if (j == length(toExtract)) {"... Finishing up..."} else {""},
                        "  |  Elapsed Time: ",
                        diffTime_to_display(extractMixTVPhi_End_Time - extractMixTVPhi_start_Time))
    if (j == length(toExtract)) {maxCatLen <<- nchar(strDisplay) - 2}
    cat(strDisplay)
  }
  return(OBJ_arr)
}) %>% do.call(c, .) %>% data.table(MIXTV_PHI_OFV = .)
extractMixTVPhi_End_Time <- Sys.time()
cat(paste0("\r\t", strrep(" ", maxCatLen)))
strDisplay <- paste0("\r\tExtracting from MixTV.phi... ", j, "/", length(toExtract),
                     if (j == length(toExtract)) {"... Finished!"} else {""},
                     "  |  Elapsed Time: ",
                     diffTime_to_display(extractMixTVPhi_End_Time - extractMixTVPhi_start_Time))
cat(strDisplay)

# > Extracting from MixTV.phm ####
mixTVPhmMissingArr <- numeric()
setwd(paste0(baseWd, "/FIT_MIXTV"))
cat("\n")
j <- 0
extractMixTVPhm_start_Time <- Sys.time()
tab_MixTV_phm <- lapply(toExtract, function(i) {
  j <<- j + 1
  tryCatch({
    thisN = dt_param[i, SampleSize]
    tryCatch({
      tabMixTVPhm = fread(fitMixTVPhmPathArr[i], header = T, skip = 1)
      if (nrow(tabMixTVPhm) != thisN * 2) {stop()}
    }, error = function(e) {stop("tabErr")})
    tab1 = tabMixTVPhm[SUBPOP == 1, .(MIXTV_PMIX1 = PMIX, MIXTV_PHM_OFV1 = OBJ)]
    tab2 = tabMixTVPhm[SUBPOP == 2, .(MIXTV_PMIX2 = PMIX, MIXTV_PHM_OFV2 = OBJ)]
    tab = cbind(tab1, tab2)
  }, error = function(e) {
    if (e$message == "tabErr") {
      tab = data.table(MIXTV_PMIX1 = rep(NA, thisN), MIXTV_PHM_OFV1 = rep(NA, thisN),
                       MIXTV_PMIX2 = rep(NA, thisN), MIXTV_PHM_OFV2 = rep(NA, thisN))
      assign("tab", tab, envir = parent.frame(4))
      assign("mixTVPhmMissingArr", c(mixTVPhmMissingArr, i), envir = .GlobalEnv)
    }
  })
  extractMixTVPhm_End_Time <- Sys.time()
  if (j%%100 == 0 || j == length(toExtract)) {
    strDisplay = paste0("\r\tExtracting from MixTV.phm... ", j, "/", length(toExtract),
                        if (j == length(toExtract)) {"... Finishing up..."} else {""},
                        "  |  Elapsed Time: ",
                        diffTime_to_display(extractMixTVPhm_End_Time - extractMixTVPhm_start_Time))
    if (j == length(toExtract)) {maxCatLen <<- nchar(strDisplay) - 2}
    cat(strDisplay)
  }
  return(tab)
}) %>% rbindlist
extractMixTVPhm_End_Time <- Sys.time()
cat(paste0("\r\t", strrep(" ", maxCatLen)))
strDisplay <- paste0("\r\tExtracting from MixTV.phm... ", j, "/", length(toExtract),
                     if (j == length(toExtract)) {"... Finished!"} else {""},
                     "  |  Elapsed Time: ",
                     diffTime_to_display(extractMixTVPhm_End_Time - extractMixTVPhm_start_Time))
cat(strDisplay)

# > Extracting from MixTV.txt ####
mixTVTxtMissingArr <- numeric()
setwd(paste0(baseWd, "/FIT_MIXTV"))
cat("\n")
j <- 0
extractMixTVTxt_start_Time <- Sys.time()
tab_MixTV_txt <- lapply(toExtract, function(i) {
  j <<- j + 1
  tryCatch({
    thisN = dt_param[i, SampleSize]
    tryCatch({
      tabMixTVTxt = fread(fitMixTVTxtPathArr[i], header = T, skip = 1)
      tabID = tabMixTVTxt[, .SD[1], ID]
      if (tabID[, ID %>% equals(seq(thisN)) %>% all %>% not]) {stop()}
    }, error = function(e) {stop("tabErr")})
    tab = tabID[, .(MIXTV_CL = CL, MIXTV_VD = VD, MIXTV_KA = KA)]
  }, error = function(e) {
    if (e$message == "tabErr") {
      tab = data.table(MIXTV_CL = rep(NA, thisN), MIXTV_VD = rep(NA, thisN), MIXTV_KA = rep(NA, thisN))
      assign("tab", tab, envir = parent.frame(4))
      assign("mixTVTxtMissingArr", c(mixTVTxtMissingArr, i), envir = .GlobalEnv)
    }
  })
  extractMixTVTxt_End_Time <- Sys.time()
  if (j%%100 == 0 || j == length(toExtract)) {
    strDisplay = paste0("\r\tExtracting from MixTV.txt... ", j, "/", length(toExtract),
                        if (j == length(toExtract)) {"... Finishing up..."} else {""},
                        "  |  Elapsed Time: ",
                        diffTime_to_display(extractMixTVTxt_End_Time - extractMixTVTxt_start_Time))
    if (j == length(toExtract)) {maxCatLen <<- nchar(strDisplay) - 2}
    cat(strDisplay)
  }
  return(tab)
}) %>% rbindlist
extractMixTVTxt_End_Time <- Sys.time()
cat(paste0("\r\t", strrep(" ", maxCatLen)))
strDisplay <- paste0("\r\tExtracting from MixTV.txt... ", j, "/", length(toExtract),
                     if (j == length(toExtract)) {"... Finished!"} else {""},
                     "  |  Elapsed Time: ",
                     diffTime_to_display(extractMixTVTxt_End_Time - extractMixTVTxt_start_Time))
cat(strDisplay)

# > Extracting from NoMix.phi ####
noMixPhiMissingArr <- numeric()
setwd(paste0(baseWd, "/FIT_NOMIX"))
cat("\n")
j <- 0
extractNoMixPhi_start_Time <- Sys.time()
tab_NoMix_phi <- lapply(toExtract, function(i) {
  j <<- j + 1
  tryCatch({
    thisN = dt_param[i, SampleSize]
    tryCatch({
      tabNoMixPhi = fread(fitNoMixPhiPathArr[i], header = T, skip = 1)
      if (nrow(tabNoMixPhi) != thisN) {stop()}
    }, error = function(e) {stop("tabErr")})
    OBJ_arr = tabNoMixPhi[, OBJ]
  }, error = function(e) {
    if (e$message == "tabErr") {
      assign("OBJ_arr", rep(NA, thisN), envir = parent.frame(4))
      assign("noMixPhiMissingArr", c(noMixPhiMissingArr, i), envir = .GlobalEnv)
    }
  })
  extractNoMixPhi_End_Time <- Sys.time()
  if (j%%100 == 0 || j == length(toExtract)) {
    strDisplay = paste0("\r\tExtracting from NoMix.phi... ", j, "/", length(toExtract),
                        if (j == length(toExtract)) {"... Finishing up..."} else {""},
                        "  |  Elapsed Time: ",
                        diffTime_to_display(extractNoMixPhi_End_Time - extractNoMixPhi_start_Time))
    if (j == length(toExtract)) {maxCatLen <<- nchar(strDisplay) - 2}
    cat(strDisplay)
  }
  return(OBJ_arr)
}) %>% do.call(c, .) %>% data.table(NOMIX_PHI_OFV = .)
extractNoMixPhi_End_Time <- Sys.time()
cat(paste0("\r\t", strrep(" ", maxCatLen)))
strDisplay <- paste0("\r\tExtracting from NoMix.phi... ", j, "/", length(toExtract),
                     if (j == length(toExtract)) {"... Finished!"} else {""},
                     "  |  Elapsed Time: ",
                     diffTime_to_display(extractNoMixPhi_End_Time - extractNoMixPhi_start_Time))
cat(strDisplay)

# > Extracting from NoMix.txt ####
noMixTxtMissingArr <- numeric()
setwd(paste0(baseWd, "/FIT_NOMIX"))
cat("\n")
j <- 0
extractNoMixTxt_start_Time <- Sys.time()
tab_NoMix_txt <- lapply(toExtract, function(i) {
  j <<- j + 1
  tryCatch({
    thisN = dt_param[i, SampleSize]
    tryCatch({
      tabNoMixTxt = fread(fitNoMixTxtPathArr[i], header = T, skip = 1)
      tabID = tabNoMixTxt[, .SD[1], ID]
      if (tabID[, ID %>% equals(seq(thisN)) %>% all %>% not]) {stop()}
    }, error = function(e) {stop("tabErr")})
    tab = tabID[, .(
      NOMIX_CL = CL, NOMIX_VD = VD, NOMIX_KA = KA
    )]
  }, error = function(e) {
    if (e$message == "tabErr") {
      tab = data.table(
        NOMIX_CL = rep(NA, thisN), NOMIX_VD = rep(NA, thisN), NOMIX_KA = rep(NA, thisN)
      )
      assign("tab", tab, envir = parent.frame(4))
      assign("noMixTxtMissingArr", c(noMixTxtMissingArr, i), envir = .GlobalEnv)
    }
  })
  extractNoMixTxt_End_Time <- Sys.time()
  if (j%%100 == 0 || j == length(toExtract)) {
    strDisplay = paste0("\r\tExtracting from NoMix.txt... ", j, "/", length(toExtract),
                        if (j == length(toExtract)) {"... Finishing up..."} else {""},
                        "  |  Elapsed Time: ",
                        diffTime_to_display(extractNoMixTxt_End_Time - extractNoMixTxt_start_Time))
    if (j == length(toExtract)) {maxCatLen <<- nchar(strDisplay) - 2}
    cat(strDisplay)
  }
  return(tab)
}) %>% rbindlist
extractNoMixTxt_End_Time <- Sys.time()
cat(paste0("\r\t", strrep(" ", maxCatLen)))
strDisplay <- paste0("\r\tExtracting from NoMix.txt... ", j, "/", length(toExtract),
                     if (j == length(toExtract)) {"... Finished!"} else {""},
                     "  |  Elapsed Time: ",
                     diffTime_to_display(extractNoMixTxt_End_Time - extractNoMixTxt_start_Time))
cat(strDisplay)

# >>> Writing individual extract file ####
cat("\n\t> Writing individual extract file ...")
writing_tab_ind_start_time <- Sys.time()
tab_ind <- cbind(
  data.table(datasetID = rep(toExtract, dt_param[toExtract, SampleSize])),
  tab_Mix_txt,
  tab_Mix_phi,
  tab_Mix_phm,
  tab_NoMix_txt,
  tab_NoMix_phi,
  tab_MixTV_txt,
  tab_MixTV_phi,
  tab_MixTV_phm
)
setwd(baseWd)
if (is.null(testRun)) {
  ind_extract_file_name <- paste0(projName, "_INDIVIDUAL_EXTRACT.csv")
  if (file.exists(dat_extract_file_name) && file.exists(ind_extract_file_name) && file.exists("MISSING_TEMP.csv")) {
    cur_tab_ind <- fread(ind_extract_file_name, header = T)
    temp_rowID <- lapply(unique(tab_ind[, datasetID]), function(ID) {which(cur_tab_ind[, datasetID] == ID)}) %>% do.call(c, .)
    cur_tab_ind[temp_rowID] <- tab_ind
    tab_ind <- cur_tab_ind
    rm(temp_rowID)
  }
} else {
  IndTestSerial <- DatTestSerial
  ind_extract_file_name <- paste0(projName, "_INDIVIDUAL_EXTRACT_TEST", sprintf("%02d", IndTestSerial), ".csv")
}
write_csv(tab_ind, ind_extract_file_name)
writing_tab_ind_End_time <- Sys.time()
cat(paste0(" Done !  |  Elapsed Time: ", diffTime_to_display(writing_tab_ind_End_time - writing_tab_ind_start_time)))

cat("\n\n")
cat("****************************************************************\n")
cat("********************* Final time reporting *********************\n")
cat("****************************************************************\n")

#+++++++++++++++++++++++++++++++++++++++#
#### +++++Final time reporting+++++ #####
#+++++++++++++++++++++++++++++++++++++++#

total_End_Time <- Sys.time()
total_run_Time_str <- diffTime_to_display(total_End_Time - total_start_Time)
cat(paste0("\nRun finished!  |  Total elapsed time: ", total_run_Time_str))