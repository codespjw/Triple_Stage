clc; clear; close all;

%%
optb = bodeoptions;
optb.PhaseWrapping = 'off';
optb.grid = 'on';
optb.FreqUnits = 'Hz';
optb.xlim = [1e1 1e5];
%%
load vcm_tf

Fs = 348*120;
freq = vcm_tf.freq; % 48Hz ~ 41760Hz
response = vcm_tf.xferFunction(:,4);
idx = ~isnan(response);
response = response(idx);
freq = freq(idx);

mag = abs(response);
respFrd = frd(response,freq,'FrequencyUnit','Hz');
weight = 1./mag;

minFreq = min(freq);
maxFreq = 6500;

idx = freq > maxFreq | freq < minFreq;
weight(idx) = 0;
opt = n4sidOptions('Focus',weight);
n = 15;
model = n4sid(respFrd,n,opt);
weight1 = respFrd;
weight1.ResponseData(1,1,:) = weight;
model1 = fitfrd(respFrd,4,2,weight1);
zpk(model1)
figure;
bode(respFrd,model1,optb);
legend('meas.','model')
%%

figurename('vcm_tf');
subplot(211),
semilogx(vcm_tf.freq,20*log10(abs(vcm_tf.xferFunction))); grid on;
subplot(212),
semilogx(vcm_tf.freq,180/pi*angle(vcm_tf.xferFunction)); grid on;

%%
load ma_tf
figurename('ma_tf');
subplot(211),
semilogx(ma_tf.freq,20*log10(abs(ma_tf.xferFunction(:,1:1:10)))); grid on;
subplot(212),
semilogx(ma_tf.freq,180/pi*angle(ma_tf.xferFunction(:,1:1:10))); grid on;

%%
Fs = 348*120;
freq = ma_tf.freq; % 96Hz ~ 41760Hz
response = ma_tf.xferFunction(:,7);
idx = ~isnan(response);
response = response(idx);
freq = freq(idx);

mag = abs(response);
respFrd = frd(response,freq,'FrequencyUnit','Hz');
weight = 1./mag;

% minFreq = min(freq);
% maxFreq = max(freq);

minFreq = 300;
maxFreq = 40000;

idx = freq > maxFreq | freq < minFreq;
weight(idx) = 0;
opt = n4sidOptions('Focus',weight);
n = 15;
model = n4sid(respFrd,n,opt);
weight1 = respFrd;
weight1.ResponseData(1,1,:) = weight;
model1 = fitfrd(respFrd,4,2,weight1);
zpk(model1)

figure;
bode(respFrd,model1,optb);
legend('meas.','model')

%%
load MPK_cont_notch
figurename('sysCv');
bode(sysCv,optb)

figurename('sysCm');
bode(sysCm,optb)

figurename('sysNm');
bode(sysNm,optb)

figurename('sysNv');
bode(sysNv,optb)

figurename('sysF');
bode(sysF,optb)




