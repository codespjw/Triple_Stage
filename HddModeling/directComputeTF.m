%% clear workspace
clc; clear; close all;

% bodeOptions
bodeOpt = bodeoptions;
bodeOpt.grid = 'on';
bodeOpt.FreqUnits = 'Hz';
bodeOpt.PhaseWrapping = 'off';
bodeOpt.xlim = [1e0 1e5];

%% Public Parameters
% system para
m1n = 1.0; m2n = 0.2; m3n = 0.02; m4n = 0.02;
b1n = 0.5; b2n = 0.2; b3n = 0.1; b4n = 0.1;
k1n = 1e10; k2n = 1e5; k3n = 1e7; k4n = 1e8;
% for freqresp
nf = 1000;
Wrad = logspace(-1,6,nf);
%% Single-Stage
syms m1 m2 b1 b2 k1 k2 x1 x2 F
syms p
syms y
eq1 = ('y = x1 - x2');
eq2 = ('m1*p^2*x1 = F - b1*p*x1 - k1*x1 - b2*p*y - k2*y');
eq3 = ('m2*p^2*x2 = + b2*p*y + k2*y');
[x1,x2,y] = solve(eq1,eq2,eq3,x1,x2,y);
x1 = collect(x1,{p,F});
x2 = collect(x2,{p,F});
y = collect(y,{p,F});

X = [m1n m2n b1n b2n k1n k2n].';

F_X1S = collect(x1/F,F);
F_X2S = collect(x2/F,F);

response = double(subs(subs(F_X2S,[m1 m2 b1 b2 k1 k2].',X), p, 1i*Wrad));
S.VcmFrdModel = frd(response,Wrad);
figurename('Single-Stage: F -> PES'); 
bode(S.VcmFrdModel,bodeOpt);

%% Dual-stage, how to model PZT
syms m1 m2 m3 b1 b2 b3 k1 k2 k3 F1 F2 x1 x2 x3
syms p
syms y1 y2
eq1 = ('y1 = x1 - x2');
eq2 = ('y2  = x2 -  x3');
eq3 = ('m1*p^2*x1 = F1 - k1*x1 - b1*p*x1 - k2*y1 - b2*p*y1');
eq4 = ('m2*p^2*x2 = F2 + k2*y1 + b2*p*y1 - k3*y2 - b3*p*y2');
eq5 = ('m3*p^2*x3 = - F2 + k3*y2 + b3*p*y2');

[x1,x2,x3,y1,y2] = solve(eq1,eq2,eq3,eq4,eq5,x1,x2,x3,y1,y2);
x1 = collect(x1,{p,F1,F2});
x2 = collect(x2,{p,F1,F2});
x3 = collect(x3,{p,F1,F2});
y1 = collect(y1,{p,F1,F2});
y2 = collect(y2,{p,F1,F2});

X = [m1n m2n m3n b1n b2n b3n k1n k2n k3n].';

% VCM input to PES: F1 --> x3
F1_X3S = collect(x3/F1,p);

response = double(subs(subs(F1_X3S,[m1 m2 m3 b1 b2 b3 k1 k2 k3 F2].',[X;0]), p, 1i*Wrad));
D.VcmFrdModel = frd(response,Wrad);
figurename('Dual-Stage: F1 -> PES'); 
bode(D.VcmFrdModel,bodeOpt);

% MA input to PES: F2 --> x3
F2_X3S = collect(x3/F2,p);

response = double(subs(subs(F2_X3S,[m1 m2 m3 b1 b2 b3 k1 k2 k3 F1].',[X;0]), p, 1i*Wrad));
D.MaFrdModel = frd(response,Wrad);
figurename('Dual-Stage: F2 -> PES'); 
bode(D.MaFrdModel,bodeOpt);

%% Triple-stage, how to model PZT
syms m1 m2 m3 m4 b1 b2 b3 b4 k1 k2 k3 k4 F1 F2 F3 x1 x2 x3 x4
syms p
syms y1 y2 y3
eq1 = ('y1 = x1 - x2');
eq2 = ('y2  = x2 -  x3');
eq3 = ('y3  = x3 -  x4');
eq4 = ('m1*p^2*x1 = F1 - k1*x1 - b1*p*x1 - k2*y1 - b2*p*y1');
eq5 = ('m2*p^2*x2 = F2 + k2*y1 + b2*p*y1 - k3*y2 - b3*p*y2');
eq6 = ('m3*p^2*x3 = - F2 + F3 + k3*y2 + b3*p*y2 - k4*y3 - b4*p*y3');
eq7 = ('m4*p^2*x4 = - F3 + k4*y3 + b4*p*y3');

[x1,x2,x3,x4,y1,y2,y3] = solve(eq1,eq2,eq3,eq4,eq5,eq6,eq7,x1,x2,x3,x4,y1,y2,y3);
x1 = collect(x1,{p,F1,F2,F3});
x2 = collect(x2,{p,F1,F2,F3});
x3 = collect(x3,{p,F1,F2,F3});
x4 = collect(x4,{p,F1,F2,F3});
y1 = collect(y1,{p,F1,F2,F3});
y2 = collect(y2,{p,F1,F2,F3});
y3 = collect(y3,{p,F1,F2,F3});

X = [m1n m2n m3n m4n b1n b2n b3n b4n k1n k2n k3n k4n].';

% VCM input to PES: F1 --> x4
F1_X4S = collect(x4/F1,p);

response = double(subs(subs(F1_X4S,[m1 m2 m3 m4 b1 b2 b3 b4 k1 k2 k3 k4 F2 F3].',[X;0;0]), p, 1i*Wrad));
T.VcmFrdModel = frd(response,Wrad);
figurename('Triple-Stage: F1 -> PES'); 
bode(T.VcmFrdModel,bodeOpt);

% MA1 input to PES: F2 --> x4
F2_X4S = collect(x4/F2,p);
response = double(subs(subs(F2_X4S,[m1 m2 m3 m4 b1 b2 b3 b4 k1 k2 k3 k4 F1 F3].',[X;0;0]), p, 1i*Wrad));
T.Ma1FrdModel = frd(response,Wrad);
figurename('Triple-Stage: F2 -> PES'); 
bode(T.Ma1FrdModel,bodeOpt);

% MA2 input to PES: F3 --> x4
F3_X4S = collect(x4/F3,p);
response = double(subs(subs(F3_X4S,[m1 m2 m3 m4 b1 b2 b3 b4 k1 k2 k3 k4 F1 F2].',[X;0;0]), p, 1i*Wrad));
T.Ma2FrdModel = frd(response,Wrad);
figurename('Triple-Stage: F3 -> PES'); 
bode(T.Ma2FrdModel,bodeOpt);

%% Plot together
figurename('Comparison: VCM');
bode(S.VcmFrdModel,D.VcmFrdModel,T.VcmFrdModel,bodeOpt);
legend('S','D','T')

figurename('Comparison: Ma1');
bode(D.MaFrdModel,T.Ma1FrdModel,bodeOpt);
legend('D','T')

