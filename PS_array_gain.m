function [Score_true,Score_apx,BSR]=PS_array_gain(ofdm_paras,N,AtDs)
B=ofdm_paras.B;
Nc=ofdm_paras.Nc;
fc=ofdm_paras.fc;
Score=0;
M=10000;
for k=1:Nc
    fk=fc+(k-(Nc+1)/2)*B/Nc;
    xi=1-fk/fc;
%     tic
    temp=0;
    for m=1:M
        temp=temp+abs(sin(N*pi*AtDs*xi*m/M)/sin(pi*AtDs*xi*m/M))/N;
    end
    temp=temp/M;

    Score=Score+temp;
%     t1=toc

%     tic
%     syms x
%     expr=abs(sin(N*pi*AtDs*xi*x)/sin(pi*AtDs*xi*x))/N;
%     Score=vpaintegral(expr,[0,1])
%     t2=toc

end

Score_true=Score/Nc;

BSR=N*B/fc*AtDs/8;
syms x
fobj=abs(sinc(4*x));
Score_apx=2/(3*BSR)*vpaintegral(fobj,[0,BSR],'RelTol', 1e-4, 'AbsTol', 0)+1/3;
Score_apx=round(Score_apx,6);

end