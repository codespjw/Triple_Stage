function H = frsinglestage(X,W)

% X, a vector that define the system params
% W, frequency in Hz

% assign variables

for i = 1:numel(X)/3
    eval(['m', num2str(i), '=X((i-1)*3+1);']);
    eval(['b', num2str(i), '=X((i-1)*3+2);']);
    eval(['k', num2str(i), '=X((i-1)*3+3);']);
end

set_param('singlestageSimu/m1','mass',num2str(m1));
set_param('singlestageSimu/b1','D',num2str(b1));
set_param('singlestageSimu/k1','spr_rate',num2str(k1));
set_param('singlestageSimu/m2','mass',num2str(m2));
set_param('singlestageSimu/b2','D',num2str(b2));
set_param('singlestageSimu/k2','spr_rate',num2str(k2));

[A, B, C, D] = linmod('singlestageSimu',[],[],[1e-5 0 1]);

sys = ss(A,B,C,D);

H = freqresp(sys,W,'Hz');
H = reshape(H,[],1);

