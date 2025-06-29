function A = rounding(A,delta)
% statistic rounding
[r,c]=size(A);
for i=1:r
    for j=1:c     
        x=A(i,j);
        if abs(x)<=delta
            A(i,j)=0;
        elseif abs(1-x)<=delta
            A(i,j)=1;
        end
    end
end

end