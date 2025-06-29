function SE=compute_rate_v3(Y,BeaMtx,~)
[U,~]=qr(BeaMtx,0);
SE=norm(Y*U,'fro')^2;

end