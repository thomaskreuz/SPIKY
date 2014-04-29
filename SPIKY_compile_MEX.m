% This is the first file that should be run once the zip-package has been extracted.
% Once the MEX-files have been compiled you can run the main program SPIKY.

mex SPIKY_udists_MEX.c
mex SPIKY_ISI_MEX.c
mex SPIKY_SPIKE_MEX.c
mex SPIKY_realtimeSPIKE_MEX.c
mex SPIKY_futureSPIKE_MEX.c
mex SPIKY_SPIKEpico_MEX.c
mex SPIKY_realtimeSPIKEpico_MEX.c
mex SPIKY_futureSPIKEpico_MEX.c


%mex SPIKY_Victor_MEX.c
%[mex SPIKY_vanRossum_MEX.c]
%mex SPIKY_vanRossum_rearrange.c
%mex SPIKY_vanRossum_sort.c
