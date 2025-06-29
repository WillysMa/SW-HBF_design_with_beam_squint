function SE=compute_rate_v2(X,BeaMtx,~)

[U,~]=qr(BeaMtx,0);

SE=-norm(X*U,'fro')^2;

end