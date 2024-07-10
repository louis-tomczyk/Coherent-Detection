function noise_power = get_noise_from_pp(ai,varargin)

    if nargin == 1
        noise_power = mean(ai);
    else
        zloss       = varargin{1};
        PPEparams   = varargin{2};
        amp         = varargin{3};

        Nsteps_span = PPEparams.link.nsteps_fibre/amp.Nspan;
        loc_span    = floor(zloss/Nsteps_span);
    
        ai_trunc    = [ai(1:loc_span*Nsteps_span-1),ai((loc_span+1)*Nsteps_span+1:end)];
        noise_power = mean(ai_trunc);
    end