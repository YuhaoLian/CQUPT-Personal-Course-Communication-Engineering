load constellation

[a,b]=size(constmatrix);

NoSats=b/6;

c=[];
listSats=[];
newSatMatrix=[];
for ii=1:NoSats
    c=find(constmatrix(:,1+6*(ii-1))> 5);
    if size(c)>0, 
        listSats=[listSats; ii]; 
        newSatMatrix=[newSatMatrix, constmatrix(:,1+6*(ii-1):3+6*(ii-1))];
    end
    c=[];
end

save newSatMatrix newSatMatrix