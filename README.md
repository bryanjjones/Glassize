# Glassize
DSC glass transition analysis
Analyze text file containng data from a TA DSC Q1000 V9.9 Build 303
Differential scanning calorimetry data is anlayzed for glass transition (Tg) and onset of glass transition (Tg_onset).

Tg is identified as the highest rate of heat flow. Data is scanned for the most linear regions before and after the glass transition to find heat flow line for glassy and liquid states.
Tg_onset is identified as a deviation (10%, can be reassigned) from the glassy state toward the liquid state line.

Input files are standard output of TA DSC Q1000 V9.9 Build 303 instrument in text format. Files can be specified on with the -f flag or the script will anlayze everything in the ToProcess directory.
Output on command line is Tg and Tg_onset. As well as equations for glassy, transition, and liquid state lines. Plots are exported in png form showing heat flow vs temperature with Tg, g_onset, and glassy and liquid linear regions identified. Also a plot of the derivitave of heat flow is produced.
