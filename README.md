Main script is matched_filter_big_data.R. It is designed for cBioPortal_data_big.csv, and uses the function scoring.R. If those two files are in the working directory, it should run without a problem. Parameters can be set at the top of the script.

The script uses matched filtering with a reward vector. The idea is best explained in the Google Drive under Ted Rock > Presentations > Rotation Presentation.pptx.

The script outputs information about each of the assays, as well as a dbSNP at the center of each drop-off probe, which can be used to find flank sequences in Entrez.

All plotting was done in MATLAB. I am attaching a rough MATLAB file I used for this, but be warned that this may require some more tweaking to get the exact plots you want. Data can be moved from R to MATLAB for plotting using the R package R.matlab and the function writeMat.

If you have any questions, please contact Ted Rock by email at tedrock@bu.edu, by phone at (443)-877-7446, or on Slack.
