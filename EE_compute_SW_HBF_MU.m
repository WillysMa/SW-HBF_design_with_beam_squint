function EE_SW_HBF=EE_compute_SW_HBF_MU(SE_SW_HBF,Nt,N_rf)
factor2=1e-3;
P_LNA=136*factor2;
P_PA=268*factor2;
P_SP=19.5*factor2;
P_C=19.5*factor2;

P_M=22*factor2;
P_LO=4*factor2;
P_LPF=14*factor2;
P_BBamp=5*factor2;

Pcommon=Nt*P_PA;
P_ADC=560*factor2;
P_SW=5*factor2;
P_RF=P_M+P_LO+P_LPF+P_BBamp;

part1=Nt*N_rf*P_SW+2*N_rf*(P_RF+2*P_ADC)+Pcommon;
part2=N_rf*P_SP + (Nt+N_rf)*P_C;

Consp_SW=part1+part2;

EE_SW_HBF=SE_SW_HBF./Consp_SW;
end