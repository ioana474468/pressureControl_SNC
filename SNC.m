clc,clear,close all

%% PREGATIRE EXPERIMENT IDENTIFICARE
load raspIndicial.mat
t=out.comanda.Time(134:end)-13.3;
u=out.comanda.Data(134:end);
y=out.simout.Data(134:end)-3;
info=stepinfo(y,t,'SettlingTimeThreshold',0.05)

plot(t,u)
hold on
plot(t,y)
title("Comanda si raspunsul indicial")
xlabel("Timp (s)")
legend("Comanda","Raspuns")

%% REALIZARE SI ANALIZA EXPERIMENT IDENTIFICARE
N=9; p=2;
L_u0=2*N*p;
L_spab=p*(2^N-1);
Ts=0.3; Fs=1/Ts;
u_spab=SPAB_generator(N,p,Ts,60,10);
load dateExper.mat
y_spab=out.simout.Data(1:end-1);

plotFreq(u_spab(L_u0+1:L_u0+L_spab),Fs,'s',[],[-Fs/2 Fs/2 0 2000])
plotFreq(y_spab(L_u0+1:L_u0+L_spab),Fs,'s',[],[-Fs/2 Fs/2 0 2000])

%% IDENTIFICARE SI VALIDARE MODEL MATLAB
data=iddata(y_spab(L_u0+1:end),u_spab(L_u0+1:end),Ts)
data_d=detrend(data,0);
[b,a]=butter(1,0.2);
data_f=iddata(filter(b,a,data_d.y),data_d.u,Ts);

plotFreq(data.y(1:L_spab),Fs,'s',[],[-Fs/2 Fs/2 0 2000])
sgtitle("Iesirea nefiltrata")
plotFreq(data_f.y(1:L_spab),Fs,'s',[],[-Fs/2 Fs/2 0 2000])
sgtitle("Iesirea filtrata")

eData=data_f(1:L_spab);
vData=data_f(L_spab+1:end);
advice(eData)
nk=min(delayest(eData),1)

% ARX
M=struc(1:10,1:10,nk);
V=arxstruc(eData,vData,M);
ordin=selstruc(V,0);
%selstruc(V,'plot')
arx_best=arx(eData,ordin);

figure
compare(vData,arx_best)
figure
plot(resid(vData,arx_best))
figure
resid(vData,arx_best,'corr')

% ARMAX
modele_armax=cell(1,20);
MSE_armax=zeros(1,20);
for i=1:20
    na=randi(10);
    nb=randi(10);
    nc=randi(10);
    modele_armax{i}=armax(eData,[na nb nc nk]);
    MSE_armax(i)=mean(resid(vData,modele_armax{i}).y.^2);
end
[~,idx_armax]=min(MSE_armax);
armax_best=modele_armax{idx_armax};

figure
compare(vData,armax_best)
figure
plot(resid(vData,armax_best))
figure
resid(vData,armax_best,'corr')

% BJ
modele_bj=cell(1,20);
MSE_bj=zeros(1,20);
for i=1:20
    nb=randi(10);
    nc=randi(10);
    nd=randi(10);
    nf=randi(10);
    modele_bj{i}=bj(eData,[nb nc nd nf nk]);
    MSE_bj(i)=mean(resid(vData,modele_bj{i}).y.^2);
end
[~,idx_bj]=min(MSE_bj);
bj_best=modele_bj{idx_bj};

figure
compare(vData,bj_best)
figure
plot(resid(vData,bj_best))
figure
resid(vData,bj_best,'corr')

% OE
modele_oe=cell(1,20);
MSE_oe=zeros(1,20);
for i=1:20
    nb=randi(10);
    nf=randi(10);
    modele_oe{i}=oe(eData,[nb nf nk]);
    MSE_oe(i)=mean(resid(vData,modele_oe{i}).y.^2);
end
[~,idx_oe]=min(MSE_oe);
oe_best=modele_oe{idx_oe};

figure
compare(vData,oe_best)
figure
plot(resid(vData,oe_best))
figure
resid(vData,oe_best,'corr')

%% Corectie modele
load modele.mat
c=153;
K=info.SettlingMax/c;

% ARX
y_arx=step(arx_best);
figure
subplot(1,2,1)
step(c*arx_best);
title("Model ARX initial")
subplot(1,2,2)
step(c*arx_best*(K/y_arx(end)))
title("Model ARX corectat")

% ARMAX
y_armax=step(armax_best);
figure
subplot(1,2,1)
step(c*armax_best);
title("Model ARMAX initial")
subplot(1,2,2)
step(c*armax_best*(K/y_armax(end)))
title("Model ARMAX corectat")

% BJ
y_bj=step(bj_best);
figure
subplot(1,2,1)
step(c*bj_best);
title("Model BJ initial")
subplot(1,2,2)
step(c*bj_best*(K/y_bj(end)))
title("Model BJ corectat")

% OE
y_oe=step(oe_best);
figure
subplot(1,2,1)
step(c*oe_best);
title("Model OE initial")
subplot(1,2,2)
step(c*oe_best*(K/y_oe(end)))
title("Model OE corectat")

%% Stabilitate modele
figure
pzmap(arx_best,armax_best,bj_best,oe_best)
legend('Location','northwest')

figure
nyquist(arx_best,armax_best,bj_best,oe_best)
legend('Location','northwest')

figure
bode(arx_best,armax_best,bj_best,oe_best)
legend('Location','northwest')

%% MODEL ALES CORECTAT
nB=numel(bj_best.B)-1;
nA=numel(bj_best.F)-1;
d=0;

y_model=step(bj_best);
delay=tf([0 1],1,Ts,'Variable','z^-1');
model=tf(bj_best.B*(K/y_model(end)),bj_best.F,Ts,'Variable','z^-1')%*delay

%% Performante
% Reglare
zeta=0.9; %0.6;
tt=12;
wn=4/zeta/tt;
Hcl=c2d(tf(wn^2,[1 2*wn*zeta wn^2]),Ts,'zoh')
P=Hcl.Denominator{1}

% Urmarire
zeta=0.95; %0.9
tt=6; %8
wn=4/zeta/tt;
GT=c2d(tf(wn^2,[1 2*wn*zeta wn^2]),Ts,'zoh')

%% REGLARE + URMARIRE RST fara integrator
B=model.Numerator{1}(1:nB+d+1)
nB=numel(B)-d-1;
A=model.Denominator{1}
nA=numel(A)-1;

nS=nB+d-1;
nR=nA-1;
M=[convmtx(A',nB+d) convmtx(B',nA)]
P=[P zeros(1,size(M, 1)-length(P))];
x=M\P'

S=x(1:nS+1)'
R=x(nS+2:end)'
T=P/sum(B)

%% REGLARE + URMARIRE RST cu integrator
int=tf(1,[1 -1],Ts,'Variable','z^-1');
model_int=model*int
B=model_int.Numerator{1}(1:nB+d+1)
nB=numel(B)-d-1;
A_int=model_int.Denominator{1}
nA_int=numel(A_int)-1;

nSint=nB+d-1;
nR=nA_int-1;
M=[convmtx(A_int',nB+d) convmtx(B',nA_int)]
P=[P zeros(1,size(M, 1)-length(P))];
x=M\P'

Sint=x(1:nSint+1)';
R=x(nSint+2:end)'
S=tf(Sint,1,Ts)/int;
S=S.Numerator{1}
T=P/sum(B)

%save 'RST.mat' GT R S T model

%% Validare proces real
load RST\dateRST_r1.mat
figure
plot(out.comanda)
hold on
plot(out.simout)
xlabel("Timp (s)")
ylabel("")
legend("Comanda","Iesire")

load RST\dateRST_r2.mat
figure
plot(out.comanda)
hold on
plot(out.simout)
xlabel("Timp (s)")
ylabel("")
legend("Comanda","Iesire")

load RST\dateRST_test.mat
figure
plot(out_test.comanda)
hold on
plot(out_test.simout)
xlabel("Timp (s)")
ylabel("")
legend("Comanda","Iesire")

%% Robustete
A_tf=tf(A,1,Ts,'Variable','z^-1');
S_tf=tf(S,1,Ts,'Variable','z^-1');
P_tf=tf(P,1,Ts,'Variable','z^-1');
Hpy=A_tf*S_tf/P_tf
bode(Hpy)
deltaM=1/norm(Hpy,'inf')