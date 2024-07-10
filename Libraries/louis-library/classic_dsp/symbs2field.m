function E = symbs2field(symbols)

    [Nsymb,Nt] = size(symbols);
    E = zeros(1,Nsymb.*Nt);

    for k = 1:Nsymb
        E((k-1)*Nt+1:k*Nt) = symbols(k,:);
    end

end