function locvals = get_locvals(Mean_PP_trunc,Mean_PP_trunc_av,PPEparams,amp)

    [val_mins,loc_mins] = flmm(Mean_PP_trunc,Mean_PP_trunc_av,"min",PPEparams,amp);
    [val_maxs,loc_maxs] = flmm(Mean_PP_trunc,Mean_PP_trunc_av,"max",PPEparams,amp);
    
    N_loc_max_av        = length(loc_maxs);
    N_loc_min_av        = length(loc_mins);
    
    assert(N_loc_min_av == N_loc_max_av,...
    "The number of local minima and maxima should be equal")
    
    locvals.loc_mins    = loc_mins;
    locvals.loc_maxs    = loc_maxs;
    locvals.val_mins    = val_mins;
    locvals.val_maxs    = val_maxs;

end