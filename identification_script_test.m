numtaps = 50; mu = 0.01;
X = zeros([numtaps 1]); 
P = zeros([numtaps 1]); %Filter Koeffizienten
uk = zeros([4 1]); %Eingangssignal
yk = 0;
numiter = 500;
for k = 1:numiter,
uk = [randn(1); uk(1:end-1)]; % simulate actual plant - shift in random values
yk = 0.8187*yk + 0.1042*uk(3) + 0.1042*0.7402*uk(4); %Ausgang der Strecke

X = [uk(1); X(1:end-1)]; % simulate plant model - shift in same value as in uk
yhat = P'*X; %Faltung mit dem Filter

mu = min(mu,1/(X'*X)); % adaptive learning rate
P = P + 2*mu*(yk-yhat)*X; % adapt P
end
stem(P,'.');