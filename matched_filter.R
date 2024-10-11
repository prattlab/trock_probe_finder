# This script aims to develop a matched filter for detection of optimal
# drop-off and reference probes along TP53 based on data describing
# frequency of mutations for different patients.

start_time <- Sys.time()

# STEP 0: LOADING LIBRARIES
library(pracma)
source('C:\\Users\\Ted Rock\\OneDrive - University of Pittsburgh\\Documents\\BU\\Rotations\\Pratt Lab\\scoring.R')


# STEP 1: PARAMETERS
probe_min <- 19 # Shortest length each probe can be
probe_max <- 31 # Longest length each probe can be
space_max <- 30 # Longest space in between probes
num_probes <- 5 # Number of probe locations to search for
ref_probe_size <- 21 # Pre-allocated reference probe length
a <- 3 # Controls weighting of score


# STEP 2: DATA SETUP

# Reading table from CSV and isolating vector of amino acid locations with
# mutations associated
df <- read.csv('C:\\Users\\Ted Rock\\OneDrive - University of Pittsburgh\\Documents\\BU\\Rotations\\Pratt Lab\\cbioPortal_data.csv')
df <- df[-c(549:557), ]
df <- df %>%
  transmute(Mutation=Protein..Change, Position=3*AA.Position)
mutation_locations = df[['Position']]

# Converting mutation_locations into a vector mutations with number of
# mutations at each location (find a nicer way to do this in the future)
mutations = vector(mode="integer", length=max(mutation_locations))
for (i in 1:max(mutation_locations)){
  mutation_present = mutation_locations == i
  mutations[i] = sum(mutation_present)
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
                minpeakheight = 0.4*max(association),
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
for (i in 1:num_probes){D_m[i] = sum(mutations[drop_off_start[i]:drop_off_end[i]])}
for (i in 1:num_probes){R_m[i] = sum(mutations[ref_start[i]:ref_end[i]])}
final_probes <- data.frame(drop_off_start, drop_off_end, ref_start, ref_end, D_m, R_m, score)

end_time <- Sys.time()

print(end_time - start_time)

