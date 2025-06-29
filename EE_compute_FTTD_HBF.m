
function EE_TTD_HBF=EE_compute_FTTD_HBF(SE_FTTD_HBF,Nr,Nt,N_rf,Nfdt,Nfdr)
factor2=1e-3;
P_LNA=136*factor2;
P_PA=268*factor2;
P_SP=19.5*factor2;
P_C=19.5*factor2;

P_M=22*factor2;
P_LO=4*factor2;
P_LPF=14*factor2;
P_BBamp=5*factor2;

Pcommon=Nt*P_PA+Nr*P_LNA;
P_ADC=560*factor2;
P_SW=5*factor2;
P_RF=P_M+P_LO+P_LPF+P_BBamp;

P_FTTD=63*factor2;% ttd power consumption

part1=Pcommon+(Nt+Nr)*P_SW+2*N_rf*(P_RF+2*P_ADC)+N_rf*(Nfdt+Nfdr)*P_FTTD;
part2=N_rf*(Nfdt*P_SP + Nfdr*P_C);

Consp_FTTD=part1+part2;

EE_TTD_HBF=SE_FTTD_HBF/Consp_FTTD;
end