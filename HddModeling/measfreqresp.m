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

figurename('2nd order system: meas, model'); bode(measfrd,'r',modelFrd,'k:'); grid;
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

%
figurename('2nd order system: meas, optimized');
loglog(Wrad,measmag,'-b');hold all;
loglog(Wrad,fun(x(:,best))+ measmag,'--r');

if nrmlzFlag
    [~,idfrd] = freqRespPara(x(:,best).*X,Wrad);
else
    [~,idfrd] = freqRespPara(x(:,best),Wrad);
end
figurename('2nd order system: meas, id');
bode(measfrd,'-b',idfrd,'--r'); grid;

%% Mechanical system created by hand
format long 
s = tf('s');
m1 = 1.5; m2 = 0.5; m3 = 0.0003; m4 = 0.002; 
b1 = 0.2; b2 = 0.3;  b3 = 0.2; b4 = 0.1;
k1 = 1e9; k2 = 1e9; k3 = 1e6; k4 = 1e3;
Den = m1*m2*s^4 + (b1*m2 + b2*m1 + b2*m2)*s^3 + (b1*b2 + k1*m2 + k2*m1 + k2*m2)*s^2 + (b1*k2 + b2*k1)*s + k1*k2;
Den1 = m1*m2*s^4 + (m2*b1 + m1*b2)*s^3 + (m1*k2 + m2*k1 + b1*b2)*s^2 + (b1*k2 + b2*k1)*s + (k1*k2);

F_X1 = (m2*s^2 + b2*s + k2)/Den;
zpk(F_X1)
F_X2 = (b2*s + k2)/Den;
zpk(F_X2)
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

F_X1S = collect(x1/u,u);
F_X2S = collect(x2/u,u);

X = [0.001 0.01 3e3 1 0.1 1e6].';

nf = 1000;
Wrad = logspace(-1,6,nf);
Wrad = Wrad.';
response = double(subs(subs(F_X1S,[m1 b1 k1 m2 b2 k2].',X), p, 1i*Wrad));
%
magdB = 20*log10(abs(response));
figurename('magdB'); semilogx(Wrad,magdB); grid ;

measfrd = frd(F_X1,Wrad);
measmag = abs(reshape(measfrd.ResponseData,[],1));
weight = 1./measmag;
weight1 = measfrd;
weight1.ResponseData(1,1,:) = weight;
N = 4; RD = 2;
model = fitfrd(measfrd,N,RD,weight1);
disp('real system')
zpk(F_X1)
disp('id system')
zpk(model)
modelFrd = frd(model,Wrad);

figurename('meas, model'); bode(measfrd,'r',modelFrd,'k:'); grid;
legend('meas','model');

%% construct fun
clear x resnorm;
fun = @(x) abs(double(subs(subs(F_X1S,[m1 b1 k1 m2 b2 k2].',x.*X), p, 1i*Wrad))) - measmag;
%%
% singlestageSimu;
% H = frsinglestage(X,W);

% fun = @(p) abs(frsinglestage(p.*X,W))-abs(measFr);

options = optimoptions('lsqnonlin','Display','iter');
% options = optimoptions('lsqnonlin');
X0=ones(size(X))+rand(size(X))*0.8;
% lb= zeros(size(X))+eps;
lb= ones(size(X));
ub = ones(size(X))*3;

numiter = 10;
% parpool(5);
X0t = rand(size(X,1),numiter).*repmat(ub-lb,1,numiter) + repmat(lb,1,numiter);
parfor i = 1:numiter
    disp(i);
    [x(:,i),resnorm(:,i)] = lsqnonlin(fun,X0t(:,i),lb, ub,options);
end
[~,best]=min(resnorm);

[x(:,best) X0t(:,best)]
[norm(fun(x(:,best))) norm(fun(X0t(:,best)))]
norm(fun(ones(size(X))))

figure;
% loglog(abs(measfrd));hold all;
% loglog(abs( frsinglestage(optimresults1.x .* X, W)))

loglog(Wrad,abs(measmag),'-b');hold all;
loglog(Wrad,abs(fun(x(:,best))+ measmag),'--r');
grid;
