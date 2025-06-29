function EE_DC=EE_compute_DC(SE_DC,Nr,Nt)
factor2=1e-3;
P_LNA=136*factor2;
P_PA=268*factor2;

P_M=22*factor2;
P_LO=4*factor2;
P_LPF=14*factor2;
P_BBamp=5*factor2;

Pcommon=Nt*P_PA+Nr*P_LNA;
P_ADC=560*factor2;
P_RF=P_M+P_LO+P_LPF+P_BBamp;

Consp_DC=Pcommon+(Nr+Nt)*(P_RF+2*P_ADC);

EE_DC=SE_DC/Consp_DC;
end