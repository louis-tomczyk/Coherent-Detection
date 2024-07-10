function LinkParams = set_link(amps,varargin)

    namps   = length(fieldnames(amps));
    Ltot    = 0;

    for k = 1:namps
        Ltot = Ltot+amps.().length*amp
    end