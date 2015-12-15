%% Estimating frequency response for a Simulink model:

% Create input signal for simulation:
input = frest.Sinestream('Frequency',logspace(-3,2,30));

% Open the Simulink model:
watertank

% Specify portion of model to estimate:
io(1)=linio('watertank/PID Controller',1,'input');
io(2)=linio('watertank/Water-Tank System',1,'openoutput');

% Specify the steady state operating point for the estimation.
watertank_spec = operspec('watertank');
op = findop('watertank',watertank_spec);

% Estimate frequency response of specified blocks:
sysest = frestimate('watertank',op,io,input);
bode(sysest)

%% Validate exact linearization results using estimated frequency response of a Simulink model:

% Open the Simulink model:
watertank

% Specify portion of model to estimate:
io(1)=linio('watertank/PID Controller',1,'input');
io(2)=linio('watertank/Water-Tank System',1,'output');

% Specify operating point for linearization and estimation:
watertank_spec = operspec('watertank');
op = findop('watertank',watertank_spec);

% Linearize the model:
sys = linearize('watertank',op,io);

% Estimate the frequency response of the watertank model
input = frest.Sinestream('Frequency',logspace(-1,2,10));
[sysest,simout] = frestimate('watertank',op,io,input);

% Compare linearization and estimation results in frequency domain:
frest.simView(simout,input,sysest,sys)