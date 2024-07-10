function pp_flat = get_ppflat_pd(PPEparams,Axis,las,tx,ft,amp)

    tx_flat             = tx;
    tx_flat.pd          = set_pd(90,Axis.symbrate);
    PPEparams_flat      = set_PPEparams(PPEparams,ft,amp,tx_flat);
    pp_flat             = get_pp(PPEparams_flat,Axis,las,tx_flat,ft,amp);