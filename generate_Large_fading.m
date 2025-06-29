function large_fading=generate_Large_fading(Nu,radius,r_h)
LSE = 10.^(sqrt(64)*randn(1,Nu)/10); 
US_CDT = rand(Nu, 2)*2 - ones(Nu, 2);
for j=1:Nu
    x_pos = US_CDT(j, 1);
    y_pos = US_CDT(j, 2)*(2/sqrt(3) - (1/sqrt(3))*abs(x_pos));
    US_CDT(j,:) = [x_pos, y_pos]*radius;
end
large_fading=zeros(1,Nu);
for ii=1:Nu
    dis_term=(r_h/norm(US_CDT(ii,:)))^3.7;%一个小区基站到一个小区用户的距离
    large_fading(ii)=LSE(ii)*dis_term;
end

end