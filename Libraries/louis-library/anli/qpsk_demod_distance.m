function out = qpsk_demod_distance(in,Axis)

    n_samp  = size(in.field,1);
    n_polar = size(in.field,2);
    out     = in;

    %%% downsampling if necessary
    if n_samp > Axis.Nsymb
        % getting the number of samples per symbol
        Nsps    = n_samp/Axis.Nsymb;
        assert(is_integer(Nsps) == 1,"Number of samples per symbol should be an integer")

        % downsample to get 1 sample per symbol
        in          = ds(in,1,Nsps,1);
        n_samp      = Axis.Nsymb;
        out         = rmfield(out,"field");
        out.field   = zeros(n_samp,n_polar);
    end
    
    m   = [-3,-1,+1,+3];
    phi = exp(1i*m*pi/4);

    %%% demodulation
    for j = 1:n_polar
        for k = 1:n_samp
            amp             = abs(in.field(k,j));
            th              = fastexp(angle(in.field(k,j)));
            distance        = get_distance(phi,th,'L1');
            [~,imin]        = min(distance);
            m_hat           = m(imin);
            out.field(k,j)  = amp.*th*fastexp(-m_hat*pi/4);
        end
    end