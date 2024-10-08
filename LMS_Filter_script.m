% clear all

%Simulation parameters
t_total = 50; %s
freq = 1000; %Hz zu hohe Frequenz führt zu schlechter Konvergenz
samples = t_total*freq;

%Galvo Strecke und Regler
S_galvo = tf(5.263e09, [1 3.247e04 2.261e07 0]);
R_pid = pid(0.4, 9.2e-4, 9.2e-5, 0);

strecke = feedback(S_galvo*R_pid, 1); %Geschlossener Regelkreis von Galvo und PID
time_vector = linspace(0,t_total,samples);


% input_signal = sin(time_vector);
input_signal = ((rand(samples,1))-0.5)*100;
% input signal filter
input_filter = tf(1, [0.5 1]);
input_signal = lsim(input_filter,input_signal, time_vector);
%Ausgangssignal simulieren
output_signal = lsim(strecke,input_signal,time_vector);

%Filterparameter
num_taps = 75; % Anzahl der Filterkoeffizienten
w = zeros(num_taps, 1); % Initialisierung mit Nullen
mu = 0.08/num_taps;

%Adaptiven Filter an jedem Zeitpunkt berechnen und sein Ausgangssignal an
%diesem Zeitpunkt
adaptive_filter_out = zeros(samples,1);
for t = 1:length(time_vector)
    % w = adapt_filter_at(w, input_signal, output_signal, mu, t);
    u = output_signal(t);
    filter_out = fir_filter(input_signal,w,t);
    e = u -filter_out;
    w = lms(w,e,input_signal,t,mu);
    adaptive_filter_out(t) = filter_out;
    
end
%Fehler zwischen Soll Signal und Filter Ausgang
error_out = output_signal-adaptive_filter_out;
num_taps
end_error = mean(error_out(end-samples*0.1:end).^2)

%plotten
figure 
subplot(2,1,1);
hold on
plot(time_vector,adaptive_filter_out, 'x');
plot(time_vector,output_signal,'LineWidth',2);
title('Verlauf der Strecke und des adaptiven Filters')
ylabel('Ausgansamplitude')
xlabel('Zeit [s]')
hold off
legend('Filter Ausgang','Strecken Ausgang')
subplot(2,1,2);

semilogy(error_out.^2)
title(['Verlauf des Fehlers zwischen Strecke und Filter mit Lernrate mu = ',num2str(mu), ' und Filter taps = ', int2str(num_taps)])
ylabel('Fehler^2')
xlabel('Zeit [s]')

%% Adapt inverted System
mu_inv = mu /2;
n_inv = 20;
wc = w; %copy w
c = zeros(n_inv,1); %c should be double the length of w
c_out = zeros(samples,1);
w_out = zeros(samples,1);
%precalculate output of plant filter
for k= 1:samples
    w_out(k) = fir_filter(input_signal,wc,k);
end
for t = n_inv:length(time_vector)
    input_delayed = input_signal(t-(n_inv)+1); %delay input signal by length(c) steps for non-minimumphase system
    c_out_cur = fir_filter(w_out,c,t);%cal
    e = input_delayed - c_out_cur; 
    c = lms(c,e,input_signal,t,mu_inv);
    c_out(t) = c_out_cur;
    
end
figure
subplot(2,1,1);
hold on
plot(c_out(n_inv:end),'DisplayName','Inverted Filter Out shifted')
plot(input_signal(1:end-n_inv),'DisplayName','Input Signal')
legend;
hold off

error_out_inv = input_signal-c_out;
end_error_inv = mean(error_out_inv(end-samples*0.1:end).^2)
subplot(2,1,2);
semilogy(error_out_inv.^2)
title(['Verlauf des Fehlers zwischen Eingangssignal und Ausgang des invertierten Filters mit Lernrate mu = ',num2str(mu), ' und Filter taps = ', int2str(num_taps)])
ylabel('Fehler^2')
xlabel('Zeit [s]')

%% fir_filter function test
% c = [-0.02010411882885732
% -0.05842798004352509
% -0.061178403647821976
% -0.010939393385338943
% 0.05125096443534972
% 0.033220867678947885
% -0.05655276971833928
% -0.08565500737264514
% 0.0633795996605449
% 0.31085440365663597
% 0.4344309124179415
% 0.31085440365663597
% 0.0633795996605449
% -0.08565500737264514
% -0.05655276971833928
% 0.033220867678947885
% 0.05125096443534972
% -0.010939393385338943
% -0.061178403647821976
% -0.05842798004352509
% -0.02010411882885732
% ]
% 
% imp = zeros(length(c),1);
% imp(1) = 1;
% out = zeros(length(c),1);
% for i = 1:length(c)
%     out(i) = fir_filter(imp,c,i);
% end
% plot(out)

function output = fir_filter(input,coeff,position)
%calculate the current filter output for the input[position]
N = length(coeff);
output = 0;
for i = 1:N
    if position - i + 1 > 0 
        output = output + coeff(i) * input(position - i + 1);
    end
end
end

%Perform LMS Algorithm
%input_signal needs to include the last length(w) values
function w1 = lms(w0,error,input_signal,index,mu)
w1 = w0;
n = length(w0);
if index-n > 0
    xn = input_signal(index:-1:index-n+1);
    w1 = w0 + 2 * mu * error * xn;
end
end


%%Adapt filter over whole input and desired signal
function adaptedCoeff = adapt_filter_whole(startCoeff, inputSignal, desiredOutputSignal, mu)
    N = length(inputSignal);
    c =  startCoeff;
    for index = 1:N
        u = desiredOutputSignal(index);
        filter_out = fir_filter(inputSignal,c,index);
        e =  u - filter_out;
        c = lms(c,e,inputSignal,index,mu);
    end
    adaptedCoeff = c;
end

%%adapt filter until a position in the input and desired signal
function adaptedCoeff = adapt_filter_until(startCoeff, inputSignal, desiredOutputSignal, mu, position)

    c =  startCoeff;
    for index = 1:position
        u = desiredOutputSignal(index);
        filter_out = fir_filter(inputSignal,c,index);
        e =  u - filter_out;
        c = lms(c,e,inputSignal,index,mu);
    end
    adaptedCoeff = c;
end

function adaptedCoeff = adapt_filter_at(startCoeff, inputSignal, desiredOutputSignal, mu, position)

    c =  startCoeff;
    n = length(c);
    u = desiredOutputSignal(position);
    filter_out = fir_filter(inputSignal,c,position);
    e =  u - filter_out;
    if position - n > 0
    for index = position-n:position       
        c = lms(c,e,inputSignal,index,mu);
    end
    end
    adaptedCoeff = c;
end
