# This function scores pairs of drop-off and reference probes for use in
# optimizing an assay of TP53. Written for the Pratt Lab at BU.
# Edwin Rock 10/2/24

scoring <- function(mutations, patients, drop_off_start, drop_off_end, ref_start, ref_end, a){
  
  # The current score function is D_m/(a*R_m + D_L*R_L*abs(D_L+R_L+spacing - 89.9))
  # D_m is mutations covered by drop off probe
  # R_m is mutations covered by reference probe
  # D_L is length of drop off probe
  # R_L is length of reference probe
  # a is a weight to change priority of mutations covered vs probe length
  # spacing is the distance between the end of one probe and the start of the other
  
  D_m = sum(mutations[drop_off_start:drop_off_end])
  R_m = sum(mutations[ref_start:ref_end])
  D_L = drop_off_end - drop_off_start + 1
  R_L = ref_end - ref_start + 1
  spacing = abs(ref_start - drop_off_end) - 1
  if (ref_start > drop_off_end) spacing = ref_start - drop_off_end - 1
  else spacing = drop_off_start - ref_end - 1
  
  # D_m/(a*R_m + D_L*R_L*abs(D_L+R_L+spacing-89.9))

  # Trying a new scoring function that includes patients and changes spacing reward from linear abs to gaussian:
  # D_m*D_p*exp(-(P_L-90)^2/1000)/(a*R_m*R_p + D_L*R_L)
  # D_p is unique patients covered by drop off probe
  # R_p is unique patients covered by reference probe

  D_p = length(unique(unlist(patients[drop_off_start:drop_off_end])))
  R_p = length(unique(unlist(patients[ref_start:ref_end])))

  # D_m*D_p*exp(-(D_L+R_L+spacing-90)^2/1000)/(a*R_m*R_p + D_L*R_L)
  
  # Trying something new once again: instead of a Gaussian, the reward for pair length is a negative sigmoid. Lower pair
  # lengths should not really be punished, but higher ones still should be. Originally I was looking for a Poisson or
  # skewed Gaussian, but these were hard to tune well and there's really no reason to work to punish PL = 20 (especially
  # because that case wil likely never exist.)
  
  D_m*D_p/(a*R_m*R_p + D_L*R_L*(1+exp(0.2*(D_L+R_L+spacing-100))))
}