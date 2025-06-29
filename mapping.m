function W_int = mapping(A,Ns)
% statistic rounding
tol=1e-4;
[r,c]=size(A);
W_int=zeros(r,c);
for i=1:r
    for j=1:c     
        x=A(i,j);
        alphabet=[floor(x),floor(x)+1];
        prob=[1-(x-floor(x)), x-floor(x)];
        W_int(i,j)=randsrc(1,1,[alphabet;prob]); % stochastic rounding

    end
end
account=0;
while rank(W_int,tol) < Ns && account <=1e+3
    account=account+1;
    for i=1:r
        for j=1:c 
            x=A(i,j);
            alphabet=[floor(x),floor(x)+1];
            prob=[1-(x-floor(x)), x-floor(x)];
            W_int(i,j)=randsrc(1,1,[alphabet;prob]); % stochastic rounding
        end
    end
    ccc=1;
end%->while
ccc=1;
end