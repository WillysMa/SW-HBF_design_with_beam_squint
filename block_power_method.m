function block_eigen_vector=block_power_method(G,N_rf,Np)
[N,~]=size(G);
V=ones(N,N_rf);
for i=1:Np
   B=G*V;
   [Q,~]=qr(B);
   V=Q(:,1:N_rf);
end
block_eigen_vector=V;
end
% complexity:Np(2N^2 N_rf-N N_rf+ 2N N_rf^2)