function R = resizer(Rm,acc,sc)

[nrlo, nclo] = size(Rm);
Nlo = nrlo*nclo;

% hi-res dimensions:
[nrhi, nchi] = size(acc);
Nhi = nrhi * nchi;

% size of blocks
br = nrhi / nrlo;
bc = nchi / nclo;

% Scaling factor based on ratio of cell areas
scaling = Nhi / Nlo;

% Create hi-res runoff map
R = zeros(size(acc),'single');

% fill it with blocks of lo-res runoff
for r=1:nrlo
    r1 = (r-1)*br+1;
    r2 = r*br;
    for c=1:nclo
        c1 = (c-1)*bc+1;
        c2 = c*bc;
        if sc==1
            R(r1:r2, c1:c2) = Rm(r,c) / scaling;
        else
            R(r1:r2, c1:c2) = Rm(r,c);
        end
    end;
end;

end