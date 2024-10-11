# This function scores pairs of drop-off and reference probes for use in
# optimizing an assay of TP53. Written for the Pratt Lab at BU.
# Edwin Rock 10/2/24

scoring <- function(mutations, drop_off_start, drop_off_end, ref_start, ref_end, a){
  
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
  
  D_m/(a*R_m + D_L*R_L*abs(D_L+R_L+spacing-89.9))
}