/*
 * SPIKE-distance (Kreuz) picewise constant version
 */

#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[],
int nrhs, const mxArray *prhs[])

{
    #define bi_spike_diffs_out plhs[0]
    
    #define num_pairs_in prhs[0]
    #define run_pico_lengths_ruc_in prhs[1]
    #define num_trains_in prhs[2]
    #define foll_spikes_in prhs[3]
    #define isi_pos_in prhs[4]
    #define prev_spikes_indy_in prhs[5]
    #define ints_in prhs[6]
    #define foll_spikes_indy_in prhs[7]
    #define prev_spikes_in prhs[8]
    #define udists_in prhs[9]
    
    int *num_pairs, *run_pico_lengths_ruc, *num_trains, *prev_spikes_indy, *foll_spikes_indy, pac = 0, sac, trac1, trac2, M;
    float *bi_spike_diffs, *udists, *udists2, *isi_pos, *foll_spikes, *prev_spikes, *ints;
    const mxArray *udistsPr, *udists2Pr;
    
    num_pairs = (int *)mxGetPr(num_pairs_in);
    run_pico_lengths_ruc = (int *)mxGetPr(run_pico_lengths_ruc_in);
    num_trains = (int *)mxGetPr(num_trains_in);
    foll_spikes = (float *)mxGetPr(foll_spikes_in);
    prev_spikes_indy = (int *)mxGetPr(prev_spikes_indy_in);
    isi_pos = (float *)mxGetPr(isi_pos_in);
    ints = (float *)mxGetPr(ints_in);
    foll_spikes_indy = (int *)mxGetPr(foll_spikes_indy_in);
    prev_spikes = (float *)mxGetPr(prev_spikes_in);
    
    bi_spike_diffs_out = mxCreateNumericArray(0, 0, mxSINGLE_CLASS, mxREAL);
    mxSetM(bi_spike_diffs_out, *num_pairs);
    mxSetN(bi_spike_diffs_out, *run_pico_lengths_ruc);
    mxSetData(bi_spike_diffs_out, mxMalloc(sizeof(float) * *num_pairs * *run_pico_lengths_ruc));
    bi_spike_diffs = (float *)mxGetPr(bi_spike_diffs_out);
    
    M  = mxGetM(udists_in);
    for(trac1 = 0; trac1 < *num_trains-1; ++trac1)
        for(trac2 = trac1 + 1; trac2 < *num_trains;  ++trac2) {
            pac++;
            
            udistsPr = mxGetCell(udists_in, trac2 * M + trac1);
            udists2Pr = mxGetCell(udists_in, trac1 * M + trac2);
            
            udists = (float *)mxGetPr(udistsPr);
            udists2 = (float *)mxGetPr(udists2Pr);
            
            for(sac = 0; sac < *run_pico_lengths_ruc; ++sac)
                bi_spike_diffs[(pac-1) + *num_pairs * sac] =
                ((udists[prev_spikes_indy[trac1 + *num_trains * sac] - 1]*(foll_spikes[trac1 + *num_trains * sac] - isi_pos[sac])
                + udists[foll_spikes_indy[trac1 + *num_trains * sac] - 1]*(isi_pos[sac] - prev_spikes[trac1 + *num_trains*sac]))
                /ints[trac1 + *num_trains*sac]*ints[trac2 + *num_trains*sac]
                + (udists2[prev_spikes_indy[trac2 + *num_trains * sac] - 1]*(foll_spikes[trac2 + *num_trains*sac] - isi_pos[sac])
                +  udists2[foll_spikes_indy[trac2 + *num_trains * sac] - 1]*(isi_pos[sac] - prev_spikes[trac2 + *num_trains*sac]))
                /ints[trac2 + *num_trains*sac]*ints[trac1 + *num_trains*sac])
                /((ints[trac1 + *num_trains*sac] + ints[trac2 + *num_trains*sac])*(ints[trac1 + *num_trains*sac] + ints[trac2 + *num_trains*sac])/2);
        }
    return;
}