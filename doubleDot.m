function out = doubleDot(A,B)

out = 0;
for ii = 1 : 3
    for jj = 1 : 3
        out = out + A(ii,jj) * B(jj,ii);
    end
end
