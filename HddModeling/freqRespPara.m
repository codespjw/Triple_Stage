function [resp,respFrd] = freqRespPara(x,w)
m = x(1);
b = x(2);
k = x(3);

resp = 1./(m*(1i*w).^2 + b*(1i*w) +k);
sys = tf([1],[m b k]);
respFrd = frd(sys,w);