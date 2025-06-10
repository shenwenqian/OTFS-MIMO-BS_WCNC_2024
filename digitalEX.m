function Q = digitalEX(H,SNR,Ntx,Nrx,M,N,Ptot)

S=M*N;
I = eye(Ntx*S,Ntx*S);
Q2 = 0;

sigma_n = Ptot/(S*SNR);


for s=1:S
    K = [];
    K = H(:,(s-1)*Ntx+1:s*Ntx);

    [~,Seff1,V1] = svd(K);
    M1 = min(size(Seff1,1),size(Seff1,2));
    for i = 1:M1
        W1(1,(s-1)*M1+i) = sigma_n/(Seff1(i,i))^2;
    end
end
p1 = waterfilling(W1,Ptot);
P1 = zeros(S,M1);
for i = 1:S
    for j = 1:M1
        P1(i,j) = p1(1,(i-1)*M1+j);
    end
end

for s = 1:S
    p1 = diag(P1(s,:));
    if rank(p1) == 0
        continue;
    end
    K = [];
    K = H(:,(s-1)*Ntx+1:s*Ntx);

    [~,Seff1,V1] = svd(K);
    F_d1 = V1*p1^(1/2);
    Q2 = Q2+(K*F_d1*F_d1'*(K)')/sigma_n;

end
Q = log2(det(I+Q2*Nrx*Ntx));
end