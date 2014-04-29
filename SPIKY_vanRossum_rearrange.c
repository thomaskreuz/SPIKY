/*	Rearanging spike trains in order to boost vR,
 note they (inputs) must not have any common element*/

#include "mex.h"
#include "math.h"
#include "limits.h"

void mexFunction(int nlhs, mxArray *plhs[],
int nrhs, const mxArray *prhs[])

{
    #define spikes_out plhs[0]
    
    #define spikesA_in prhs[0]
    #define spikesB_in prhs[1]
    #define length_in prhs[2]
    
    int num_spikesA, num_spikesB, *length, i, j;
    float *spikesA, *spikesB, *spikes1, *spikes2;
    mxArray *spikes1Pr, *spikes2Pr; 
    
    spikesA = (float *)mxGetPr(spikesA_in);
    spikesB = (float *)mxGetPr(spikesB_in);

    num_spikesA = mxGetN(spikesA_in);
    num_spikesB = mxGetN(spikesB_in); 
    
    length = (int *)mxGetPr(length_in);
    
    spikes_out = mxCreateCellMatrix(1, 2);
    
    spikes1Pr = mxCreateNumericArray(0, 0, mxSINGLE_CLASS, mxREAL);
    mxSetM(spikes1Pr, 1);
    mxSetN(spikes1Pr, num_spikesA);
    mxSetData(spikes1Pr, mxMalloc(sizeof(float) * num_spikesA));
    spikes1 = (float *)mxGetPr(spikes1Pr);      
    
    spikes2Pr = mxCreateNumericArray(0, 0, mxSINGLE_CLASS, mxREAL);
    mxSetM(spikes2Pr, 1);
    mxSetN(spikes2Pr, num_spikesB);
    mxSetData(spikes2Pr, mxMalloc(sizeof(float) * num_spikesB));
    spikes2 = (float *)mxGetPr(spikes2Pr);   
    
    i = 0;
    j = 0;
    
    while (i + j < *length) {
        
        if ((spikesA[i] < spikesB[j] && i < num_spikesA)||(j == num_spikesB)) {
            spikes1[i] = j+1;
            
            if (i < num_spikesA)
                i++;
        }
        else {
            spikes2[j] = i+1;
            
            if (j < num_spikesB)
                j++;
        }

    }
    
    mxSetCell(spikes_out, 0, spikes1Pr);
    mxSetCell(spikes_out, 1, spikes2Pr);
      
    return;
}