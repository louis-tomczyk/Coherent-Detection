function Axis = get_axis(fs,Nsamp)

    Axis.df     = fs/Nsamp;
    Axis.freq   = ((-Nsamp/2:1:(Nsamp/2-1)*1)+1/2)*Axis.df;
    
    Axis.dt     = 1/fs;
    Axis.time   = (0:(Nsamp-1))*Axis.dt;

    Axis.fs     = fs;
    Axis.Nsamp  = Nsamp;

    Axis.tmax   = Axis.time(end);

    Axis.type   = "SI";

end