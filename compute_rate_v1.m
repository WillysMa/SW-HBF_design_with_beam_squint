function SE=compute_rate_v1(X,BeaMtx,coef)
[U,~]=qr(BeaMtx,0);
Values=eig(U'*X*U);
SE=real(sum(log2(1+coef*Values)));
end