N=160;
fid= fopen('fort.99','r');

for i=1:N
    for j=1:i
        A=fgetl(fid);
        L=sscanf(A,'%d %d %*c%g%*c%g%*c');
        M(L(1),L(2))=complex(L(3),L(4));
        M(L(2),L(1))=conj( M(L(1),L(2)) );
    end
end

fclose(fid);