function apo = set_gausswin(N,center,res)

    x     = 1:N;
    sigma = res/(2*log(2))*sqrt(2);
    apo   = exp(-(x-center).^2/2/sigma.^2);
    apo   = apo.';