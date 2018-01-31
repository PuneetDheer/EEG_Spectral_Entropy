function [ EIGEN, Cup  ] = CORRELATION( SPECTRUM )

C= corr(SPECTRUM,'type','Pearson');

C(isnan(C)) = 0;    
C(isinf(C)) = 0; 

EIGEN = sort(eig(C));

tmp = ones(size(C));
tmp = triu(tmp,1);  % upper matrix of tmp
Cup = C(tmp==1);

end

