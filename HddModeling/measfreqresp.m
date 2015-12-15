clc; clear; close all;

bodeOpt = bodeoptions;
bodeOpt.grid = 'on';
bodeOpt.FreqUnits = 'Hz';
bodeOpt.PhaseWrapping = 'off';
bodeOpt.xlim = [1e0 1e5];

%% Let's start from simplest one
m = 1; b = 0.2; k = 100;
s = tf('s');
Den = m*s^2 + b*s + k;
G_X1 = 1/Den;

X = [m b k].';
nf = 1000;
Wrad = logspace(-1,6,nf);
Wrad = Wrad.';
measfrd = frd(G_X1,Wrad);

N = 2; RD = 2;
model = fitfrd(measfrd,N,RD);
disp('Real system')
zpk(G_X1)
disp('fitfrd id sys')
zpk(model)
modelFrd = frd(model,Wrad);

figurename('2nd order system: meas, model'); bode(measfrd,'r',modelFrd,'k:');
%% 
measmag = abs(reshape(measfrd.ResponseData,[],1));
measagl = angle(reshape(measfrd.ResponseData,[],1));

% response = freqRespPara(X,Wrad);

lambda = 1e-1; nrmlzFlag = 1;
if nrmlzFlag
    fun = @(x) abs(freqRespPara(x.*X,Wrad)) - measmag + lambda*(angle(freqRespPara(x.*X,Wrad))-measagl);
    X0= ones(size(X)) + rand(size(X))*0.8;
    lb= zeros(size(X)) + eps;
    ub = ones(size(X)) * 3;
else
    fun = @(x) abs(freqRespPara(x,Wrad)) - measmag + lambda*(angle(freqRespPara(x,Wrad))-measagl);
    X0= X + rand(size(X))*0.8;
    lb = [0.5 0 10].';
    ub = [1.5 1 1000].';
end
options = optimoptions('lsqnonlin','Display','iter');
numiter = 10;
% parpool(5);
X0t = rand(size(X,1),numiter).*repmat(ub-lb,1,numiter) + repmat(lb,1,numiter);
for i = 1:numiter
    disp(i);
    [x(:,i),resnorm(:,i)] = lsqnonlin(fun,X0t(:,i),lb, ub,options);
%     [x(:,i),resnorm(:,i)] = fmincon(fun,X0t(:,i),lb, ub,options);
end
[~,best]=min(resnorm);

[x(:,best) X0t(:,best)]
[norm(fun(x(:,best))) norm(fun(X0t(:,best)))]


%%
figurename('2nd order system: meas, optimized');
loglog(Wrad,measmag,'-b');hold all;
loglog(Wrad,fun(x(:,best))+ measmag,'--r');

if nrmlzFlag
    [~,idfrd] = freqRespPara(x(:,best).*X,Wrad);
else
    [~,idfrd] = freqRespPara(x(:,best),Wrad);
end
figurename('2nd order system: meas, id');
bode(measfrd,'-b',idfrd,'--r');

%% Mechanical system created by hand
s = tf('s');
m1 = 1.5; m2 = 0.5; m3 = 0.0003; m4 = 0.002; 
b1 = 0.1; b2 = 0.2;  b3 = 0.1; b4 = 0.1;
k1 = 1e10; k2 = 1e4; k3 = 1e6; k4 = 1e3;
Den = m1*m2*s^4 + (b1*m2 + b2*m1 + b2*m2)*s^3 + (b1*b2 + k1*m2 + k2*m1 + k2*m2)*s^2 + (b1*k2 + b2*k1)*s + k1*k2;
F_X1 = (m2*s^2 + b2*s + k2)/Den;
F_X2 = (b2*s + k2)/Den;

% Mechnical system created in Simulink
% Single-Stage
% Triple-Stage
% triplestageSimu,singlestageSimu,dualstageSimu
ns = 2;
switch ns
    case 1 
        stageSimu = 'singlestageSimu';
    case 2 
        stageSimu = 'dualstageSimu';
    case 3 
        stageSimu = 'triplestageSimu';
end
[A, B, C, D] = linmod(stageSimu,[],[],[1e-5 0 1]);

sys = ss(A,B,C,D);
figure;
switch ns
    case 1
        bode(sys,'-b',F_X2,'--r',bodeOpt);
    case 2
        bode(sys(1,2),sys(2,2),bodeOpt);
        legend('Fv-P2','Fm-P2');
    case 3
        bode(sys(1,3),sys(2,3),sys(3,3),bodeOpt);
        legend('Fv-P3','Fm-P3','Ft-P3');
        % bode(sys(1,1),sys(1,2),sys(2,1),sys(2,2),bodeOpt);
        % legend('Fv-P1','Fv-P2','Fm-P1','Fm-P2');
end


%% Dual-stage, how to model PZT
syms m1 m2 m3 b1 b2 b3 k1 k2 k3 F1 F2 F3 x1 x2 x3
syms p
syms y1 y2
eq1 = ('y1 = x1 - x2');
eq2 = ('y2  = x2 -  x3');
eq3 = ('m1*p^2*x1 = F1 - k1*x1 - b1*p*x1 - k2*y1 - b2*p*y1');
eq4 = ('m2*p^2*x2 = F2 + k2*y1 + b2*p*y1 - k3*y2 - b2*p*y2');
eq5 = ('m3*p^2*x3 = F3 + k3*y2 + b3*p*y2');

[x1,x2,x3,y1,y2] = solve(eq1,eq2,eq3,eq4,eq5,x1,x2,x3,y1,y2);
x1 = collect(x1,{p,F1,F2,F3});
x2 = collect(x2,{p,F1,F2,F3});
x3 = collect(x3,{p,F1,F2,F3});
y1 = collect(y1,{p,F1,F2,F3});
y2 = collect(y2,{p,F1,F2,F3});

%% Find frequency response use transfer function with parameters
% single-stage with E-Block

syms m1 b1 k1 x1 m2 b2 k2 x2 u y p
eq1 = ('y = x1 - x2');
eq2 = ('m1*p^2*x1 = u - b1*p*x1 - k1*x1 - b2*p*y - k2*y');
eq3 = ('m2*p^2*x2 = + b2*p*y + k2*y');
[x1,x2,y] = solve(eq1,eq2,eq3,x1,x2,y);
x1 = collect(x1,{p,u});
x2 = collect(x2,{p,u});
y = collect(y,{p,u});

F_X1S= collect(x1/u,u);
F_X2S= collect(x2/u,u);

X = [0.001 0.01 3e3 1 0.1 1e6].';

nf = 1000;
Wrad = logspace(-1,6,nf);
Wrad = Wrad.';
response = double(subs(subs(F_X1S,[m1 b1 k1 m2 b2 k2].',X), p, 1i*Wrad));
%
magdB = 20*log10(abs(response));
figurename('magdB'); semilogx(Wrad,magdB);

measfrd = frd(F_X1,Wrad);
measmag = abs(reshape(measfrd.ResponseData,[],1));

N = 4; RD = 2;
model = fitfrd(measfrd,N,RD);
disp('real system')
zpk(F_X1)
disp('id system')
zpk(model)
modelFrd = frd(model,Wrad);

figurename('meas, model'); bode(measfrd,'r',modelFrd,'k:');

%% construct fun
clear x resnorm;
fun = @(x) abs(double(subs(subs(F_X1S,[m1 b1 k1 m2 b2 k2].',x.*X), p, 1i*Wrad))) - measmag;
%%
% singlestageSimu;
% H = frsinglestage(X,W);

% fun = @(p) abs(frsinglestage(p.*X,W))-abs(measFr);

% options = optimoptions('lsqnonlin','Display','iter');
options = optimoptions('lsqnonlin');
X0=ones(size(X))+rand(size(X))*0.8;
lb= zeros(size(X))+eps;
ub = ones(size(X))*3;

numiter = 20;
parpool(5);
X0t = rand(size(X,1),numiter).*repmat(ub-lb,1,numiter) + repmat(lb,1,numiter);
parfor i = 1:numiter
    disp(i);
    [x(:,i),resnorm(:,i)] = lsqnonlin(fun,X0t(:,i),lb, ub,options);
end
[~,best]=min(resnorm);

[x(:,best) X0t(:,best)]
[norm(fun(x(:,best))) norm(fun(X0t(:,best)))]


figure;
% loglog(abs(measfrd));hold all;
% loglog(abs( frsinglestage(optimresults1.x .* X, W)))

loglog(Wrad,abs(measmag),'-b');hold all;
loglog(Wrad,abs(fun(x(:,best))+ measmag),'--r')


