function BER = ber(patbin,patbinhat)

    err = biterr(patbin,patbinhat,'row-wise');
    BER = sum(err~=0)/length(err);

end