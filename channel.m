function [H1,H2] = channel(M,N,L,Nrx,Ntx,c,fc,fs,d)
%% Channel matrix generation 

theta = zeros(L,1);
fai = zeros(1,1);
tao = zeros(1,1);
doppler = zeros(L,1);
T = zeros(M*N,M*N);
for i=1:M*N
    T(i,i)= exp(1j*2*pi*i/(M*N));
end
J = circshift(eye(M*N,M*N),1);
FN =dftmtx(N)/sqrt(N);
FM =dftmtx(M)/sqrt(M);
FNH = conj(FN);
FMH = conj(FM);
Hbs = [];
Hqt = [];
v = round(rand(1,L)*0.8*M);
tao = round(rand(1,L)*0.8*N);
hp = rand(1,L);%Channel taps
for i=1:Ntx
    Hbs_temp = [];
    Hqt_temp = [];
    fai = pi*rand(1,1) - pi/2;
    for j=1:Nrx
        theta = pi*rand(1,1) - pi/2; 
        Gtx=zeros(M,M);
        Grx=zeros(M,M);
        Gtxn=zeros(M,M);
        Grxn=zeros(M,M);
        Ho = zeros(M*N,M*N);
        Hef = zeros(M*N,M*N);

        detataotx = (i-1)*sin(fai)*d/c;
        detataorx = (j-1)*sin(theta)*d/c;

        for ii=1:M
             Gtx(ii,ii) =exp(-1j*2*pi*detataotx*fc*(2+(ii-1)*fs/fc));
             Grx(ii,ii) =exp(-1j*2*pi*detataorx*fc*(2+(ii-1)*fs/fc));

        end

        for l=1:L
            Ht = hp(l)*exp(-1j*2*pi*tao(l)*fc)*(J^(tao(l)))*(T^(v(l)));
            Ho = kron(FN,eye(M,M))*Ht*kron(FNH,eye(M,M))+Ho;
            Hef = kron(FN,FMH*Grx*FM)*Ht*kron(FNH,FMH*Gtx*FM)+Hef;

        end
        Hbs_temp = [Hbs_temp; Hef;];
        Hqt_temp = [Hqt_temp; Ho;];

    end
    Hbs = [Hbs Hbs_temp ];
    Hqt = [Hqt Hqt_temp ];
end
%% Channel matrix rearrangement
Hbs_re = [];
Hqt_re = [];
for i=1:M*N
    Hbs_temp = [];
    Hqt_temp = [];
    for j=1:M*N
        Ho = zeros(Nrx,Ntx);
        Hef = zeros(Nrx,Ntx);
        for ii=1:Nrx
            for jj=1:Ntx
                Ho(ii,jj)=Hqt((ii-1)*M*N+i,(jj-1)*M*N+j);
                Hef(ii,jj)=Hbs((ii-1)*M*N+i,(jj-1)*M*N+j);
            end
        end

        Hbs_temp = [Hbs_temp; Hef];
        Hqt_temp = [Hqt_temp; Ho];

    end
    Hbs_re = [Hbs_re Hbs_temp ];
    Hqt_re = [Hqt_re Hqt_temp ];   
end



H1 = Hbs_re;
H2 = Hqt_re;
%% Channel matrix normalization
G = norm(H1, 'fro');
G2 = norm(H2, 'fro');
H1 = H1/G;%with and without bs
H2 = H2/G2;%without bs



end