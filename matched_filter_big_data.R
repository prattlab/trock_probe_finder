# This script aims to develop a matched filter for detection of optimal
# drop-off and reference probes along TP53 based on data describing
# frequency of mutations for different patients.

start_time <- Sys.time()

# STEP 0: LOADING LIBRARIES
library(pracma)
library(dplyr)
source('C:\\Users\\Ted Rock\\OneDrive - University of Pittsburgh\\Documents\\BU\\Rotations\\Pratt Lab\\scoring.R')


# STEP 1: PARAMETERS

probe_min <- 19 # Shortest length each probe can be
probe_max <- 31 # Longest length each probe can be
space_max <- 30 # Longest space in between probes
num_probes <- 7 # Number of probe locations to search for
ref_probe_size <- 21 # Pre-allocated reference probe length
a <- 3 # Controls weighting of score


# STEP 2: DATA SETUP

# Importing from CSV and selecting relevant columns
df <- read.csv('C:\\Users\\Ted Rock\\OneDrive - University of Pittsburgh\\Documents\\BU\\Rotations\\Pratt Lab\\cBioPortal_data_big.csv')
df <- df %>% 
  transmute(ID = ï..Sample_ID,
            start=Start.Pos,
            end=End.Pos,
            protein_change = Protein.Change,
            dbSNP = dbSNP,
            dbSNP_exist = (dbSNP != "" || is.na(dbSNP)))
df <- df[-c(5819:5837),]

# Defining a function that converts sample ID from CSV to patient ID
ID_samp_to_patient <- function(patients){
  patient_ID <- vector(mode = "character", length = length(patients))
  for (i in 1:length(patients)){
    IDs <- strsplit(patients[i], "-")[[1]]
    if (length(IDs) == 4 | length(IDs) == 5){patient_ID[i] <- IDs[3]}
    else if (length(IDs) == 6){patient_ID[i] <- IDs[4]}
  }
  patient_ID
}

# Implementing ID conversion function
patients <- ID_samp_to_patient(df$ID)

# Setting up vectors for lollipop
mutations <- vector(mode="integer", length=(max(df$end,na.rm=TRUE)-min(df$start)))
patients_mutated <- vector(mode="list", length=(max(df$end,na.rm=TRUE)-min(df$start)))
mutation_ID <- vector(mode="list", length=(max(df$end,na.rm=TRUE)-min(df$start)))
mutation_dbSNP <- vector(mode="integer", length=(max(df$end,na.rm=TRUE)-min(df$start)))

# Cleaning NA values out of mutations vector
for (i in 1:nrow(df)){
  if (is.na(df$end[i])){
    df$end[i] <- df$start[i]
  }
}

# Converting start/end positions to number of mutations/list of patients per base
for (i in min(df$start):max(df$end, na.rm=TRUE)){
  
  mutation_present <- (i >= df$start & i <= df$end)
  mutations[i-min(df$start)+1] <- sum(mutation_present)
  mutation_dbSNP[i-min(df$start)+1] <- df$dbSNP_exist[i-min(df$start)+1]
  
  if (sum(mutation_present) != 0){
    patients_mutated[[i-min(df$start)+1]] <- unique(patients[mutation_present])
    mutation_ID[[i-min(df$start)+1]] <- unique(df$protein_change[mutation_present])
  }
}

# STEP 3: REWARD VECTOR

# Creating a vector with high value for the middle points (minimum probe length)
# and lower value for distance away from the middle (moving towards maximum
# probe length)
ends = 1:(probe_max-probe_min)
reward = c(ends, rep(max(ends+1), probe_min), rev(ends))


# STEP 4: CONVOLUTION

# Convolving reward vector with mutations vector to see highest/lowest
# association
association = convolve(mutations, reward,
                       conj = TRUE,
                       type="open");

# Creating vector conv_location which matches each point in association to
# the location of the center base pair of the probe
conv_location = (2-ceiling(length(reward)/2)):(length(mutations)+floor(length(reward)/2))


## STEP 5: OPTIMIZING PROBE CENTERS

# Finding peaks along association vector
pks = findpeaks(association,
                zero = "0",
                minpeakheight = 0.15*max(association),
                minpeakdistance = 70,
                npeaks = num_probes)
I_DO <- pks[,2]

# Finding window around drop-off probe to look for center of reference
# probe and checking whether this window goes out of the length of the gene
open_window <- probe_min - 1 + 2*(probe_max - probe_min) + space_max
lower_bound <- I_DO - open_window
upper_bound <- I_DO + open_window
I_R <- vector(mode = "integer", length=num_probes)
for (i in 1:length(I_DO)){
  if (lower_bound[i] <= 0){
    lower_bound[i] = 1
  }
  if (upper_bound[i] > length(association)){
    upper_bound[i] = length(association)
  }
  I_R[i] = which.min(association[lower_bound[i]:upper_bound[i]])
}

# Converting I to location of probe center
I_DO = conv_location[I_DO]
I_R = conv_location[I_R + lower_bound - 1]


# STEP 6: OPTIMIZING PROBE SIDE LENGTHS

# probe_extensions holds each probe's optimal extension on the left and
# right in a row
left_side <- vector(mode="integer", length = num_probes)
right_side <- vector(mode="integer", length = num_probes)
probe_score <- vector(mode="numeric", length=num_probes)
probe_extensions <- data.frame(left_side, right_side)

# For each index in I_DO, extend along either side and see how this changes
# the score to maximize score
for (i in 1:num_probes){
    
  # left_right is a matrix with number of extra bases on left as columns
  # [0; 1; 2; 3 ... n] and extra bases on right as rows [0 1 2 3 ... n].
  # Values are the score for that setup
  left_right = matrix(nrow=probe_max-probe_min+1,
                      ncol=probe_max-probe_min+1)

  # j is bases on left side, k is on right side
  for (j in 0:(probe_max-probe_min)){
    for (k in 0:(probe_max-probe_min)){
        
      # Checking still under probe_max
      if (j + k > (probe_max-probe_min)) {next}

      # Scoring
      left_right[j+1, k+1] = scoring(mutations,
                                     patients_mutated,
                                     I_DO[i]-floor(probe_min/2)-j,
                                     I_DO[i]+floor(probe_min/2)+k,
                                     I_R[i]-floor(probe_min/2),
                                     I_R[i]+floor(probe_min/2),
                                     a)
    }
  }

  # Saving optimal values to probe_extensions
  left_right[is.na(left_right)] <- 0
  opt = which(left_right == max(left_right), arr.ind=TRUE)
  probe_extensions[i, 'left_side'] = opt[1] - 1
  probe_extensions[i, 'right_side'] = opt[2] - 1
  probe_extensions[i, 'probe_score'] = max(left_right)
}


# STEP 7: PUTTING IT ALL TOGETHER

# final_probes describes each probe pair's start and end and their overall score
drop_off_start <- I_DO - floor(probe_min/2) - probe_extensions[,'left_side']
drop_off_end <- I_DO + floor(probe_min/2) + probe_extensions[,'right_side']
ref_start <- I_R - floor(ref_probe_size/2)
ref_end <- I_R + floor(ref_probe_size/2)
score <- probe_extensions[,'probe_score']
D_m <- vector(mode="integer", length=num_probes)
R_m <- vector(mode="integer", length=num_probes)
D_p <- vector(mode="integer", length=num_probes)
R_p <- vector(mode="integer", length=num_probes)
D_L <- vector(mode="integer", length=num_probes)
R_L <- vector(mode="integer", length=num_probes)
dbSNP_ratio <- vector(mode="numeric", length=num_probes)
contains_R175H <- vector(mode="logical", length=num_probes)
contains_R248 <- vector(mode="logical", length=num_probes)
contains_R273 <- vector(mode="logical", length=num_probes)
for (i in 1:num_probes){
  D_m[i] = sum(mutations[drop_off_start[i]:drop_off_end[i]])
  R_m[i] = sum(mutations[ref_start[i]:ref_end[i]])
  D_p[i] = length(unique(unlist(patients_mutated[drop_off_start[i]:drop_off_end[i]])))
  R_p[i] = length(unique(unlist(patients_mutated[ref_start[i]:ref_end[i]])))
  D_L[i] = drop_off_end[i] - drop_off_start[i] + 1
  R_L[i] = ref_end[i] - ref_start[i] + 1
  dbSNP_ratio[i] = sum(mutation_dbSNP[drop_off_start[i]:drop_off_end[i]])/D_L[i]
  contains_R175H[i] = sum(grepl("R175H", mutation_ID[drop_off_start[i]:drop_off_end[i]], fixed=TRUE)) > 0
  contains_R248[i] = sum(grepl("R248", mutation_ID[drop_off_start[i]:drop_off_end[i]], fixed=TRUE)) > 0
  contains_R273[i] = sum(grepl("R273", mutation_ID[drop_off_start[i]:drop_off_end[i]], fixed=TRUE)) > 0
}
final_probes <- data.frame(drop_off_start, drop_off_end, ref_start, ref_end, score, D_m, R_m, D_p, R_p, D_L, R_L, dbSNP_ratio, contains_R175H, contains_R248, contains_R273)

# STEP 8: FINDING MIDPOINT TO PASS TO GBLOCK FINDER
midpoint_dbSNP <- vector(mode="list", length=num_probes)
for (i in 1:num_probes){
  
  # Finding approximate middle of drop off probe and converting to indices in df
  d <- 0
  midpoint_found = FALSE
  while (midpoint_found == FALSE){
    midpoint=round(mean(c(final_probes$drop_off_start[i], final_probes$drop_off_end[i])))+min(df$start)-1 + (2*(d%%2)-1)*ceiling(d/2)
    midpoint_df=which(df$start==midpoint)
    if (sum(df$dbSNP_exist[midpoint_df]) > 0) {midpoint_found=TRUE}
    d <- d + 1
  }
  
  # Finding dbSNP in df
  for (j in midpoint_df){
    if (df$dbSNP_exist[j]){
      midpoint_dbSNP[i] <- df$dbSNP[j]
      break
    }
  }
}

end_time <- Sys.time()

print(end_time - start_time)

