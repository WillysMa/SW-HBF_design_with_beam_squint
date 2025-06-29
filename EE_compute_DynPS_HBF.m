function EE_PS_HBF=EE_compute_DynPS_HBF(SE_PS_HBF,Nr,Nt,N_rf,bit_ps)
factor2=1e-3;
P_LNA=136*factor2;
P_PA=268*factor2;
P_SP=19.5*factor2;
P_C=19.5*factor2;
P_SW=5*factor2;

P_M=22*factor2;
P_LO=4*factor2;
P_LPF=14*factor2;
P_BBamp=5*factor2;

Pcommon=Nt*P_PA+Nr*P_LNA;
P_ADC=560*factor2;
switch bit_ps
    case 4
        P_PS=42*factor2;
    case 2
        P_PS=21*factor2; 
    case 1
        P_PS=10*factor2;     
end
P_RF=P_M+P_LO+P_LPF+P_BBamp;

part1=(Nt+Nr)*(P_SW+P_PS)+2*N_rf*(P_RF+2*P_ADC)+Pcommon;
part2=N_rf*(P_SP + P_C);

Consp_PS=part1+part2;

EE_PS_HBF=SE_PS_HBF/Consp_PS;
end