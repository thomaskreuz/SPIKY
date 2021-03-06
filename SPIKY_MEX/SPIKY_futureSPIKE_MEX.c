/*
 * Future SPIKE-distance (only known spikes are ones that follow) (Kreuz)
 */

#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[],
int nrhs, const mxArray *prhs[])

{
    #define bi_spike_diffs_future_t_out plhs[0]
    
    #define num_pairs_in prhs[0]
    #define run_pili_lengths_ruc_in prhs[1]
    #define num_trains_in prhs[2]
    #define folli_in prhs[3]
    #define foll_indy_in prhs[4]
    #define udists_in prhs[5]
    
    int *num_pairs, *run_pili_lengths_ruc, *num_trains, *foll_indy, pac = 0, sac, trac1, trac2, M;
    float  *bi_spike_diffs_future_t, *udists, *udists2, *folli, fnormy;
    const mxArray *udistsPr, *udists2Pr;
    
    num_pairs = (int *)mxGetPr(num_pairs_in);
    run_pili_lengths_ruc = (int *)mxGetPr(run_pili_lengths_ruc_in);
    num_trains = (int *)mxGetPr(num_trains_in);
    folli = (float *)mxGetPr(folli_in);
    foll_indy = (int *)mxGetPr(foll_indy_in);
    
    bi_spike_diffs_future_t_out = mxCreateNumericArray(0, 0, mxSINGLE_CLASS, mxREAL);
    mxSetM(bi_spike_diffs_future_t_out, *num_pairs);
    mxSetN(bi_spike_diffs_future_t_out, *run_pili_lengths_ruc);
    mxSetData(bi_spike_diffs_future_t_out, mxMalloc(sizeof(float) * *num_pairs * *run_pili_lengths_ruc));
    bi_spike_diffs_future_t = (float *) mxGetPr(bi_spike_diffs_future_t_out);
    
    M  = mxGetM(udists_in);
    
    for(trac1 = 0; trac1 < *num_trains-1; ++trac1)
        for(trac2 = trac1 + 1; trac2 < *num_trains;  ++trac2) {
            pac++;
            
            udistsPr = mxGetCell(udists_in, trac2 * M + trac1);
            udists2Pr = mxGetCell(udists_in, trac1 * M + trac2);
            
            udists = (float *)mxGetPr(udistsPr);
            udists2 = (float *)mxGetPr(udists2Pr);
            
            for(sac = 0; sac < *run_pili_lengths_ruc; ++sac) {
                fnormy = folli[trac1 + *num_trains * sac] + folli[trac2 + *num_trains * sac];
                if (folli[trac1 + *num_trains * sac ] < folli[trac2 + *num_trains * sac])
                    bi_spike_diffs_future_t[(pac-1) + *num_pairs * sac] = (folli[trac2 + *num_trains * sac] - folli[trac1 + *num_trains * sac]
                    + udists2[foll_indy[trac2 + *num_trains * sac] - 1])/(2*fnormy + !fnormy);
                else
                    bi_spike_diffs_future_t[(pac-1) + *num_pairs * sac] = (folli[trac1 + *num_trains * sac] - folli[trac2 + *num_trains * sac]
                    + udists[foll_indy[trac1 + *num_trains * sac] - 1])/(2*fnormy + !fnormy);
            }
        }
    return;
}