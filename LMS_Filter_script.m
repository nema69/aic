% clear all

%Simulation parameters
t_total = 30; %s
freq = 50000; %Wenn die Frequenz erhöht wird, müssen auch die Anzahl der 
% Filterkoeffizienten erhöht werden, damit die Strecke noch abgebildet 
% werden kann. 
samples = t_total * freq;
%Galvo Strecke und Regler
% S_galvo = tf(5.263e09, [1 3.247e04 2.261e07 0]);
S_galvo = tf(5.7e09, [1.1 2.9e04 2e07 0]);
R_pid = pid(5, 5e-1, 9.2e-4, 9.2e-5);

strecke = feedback(S_galvo*R_pid, 1); %Geschlossener Regelkreis von Galvo und PID
time_vector = linspace(0,t_total,samples);


% input_signal = sin(time_vector);
input_signal = ((rand(samples,1))-0.5);
% input signal filter
% input_filter = tf(1, [1/1000 1]);
% input_signal = lsim(input_filter,input_signal, time_vector);
%Ausgangssignal simulieren
output_signal = lsim(strecke,input_signal,time_vector);

%Filterparameter
num_taps = 750; % Anzahl der Filterkoeffizienten
w = zeros(num_taps, 1); % Initialisierung mit Nullen
% mu = 0.0008;
learnRate = 0.1;
learnRate = 6*min(learnRate,1/(input_signal'*input_signal));

%Adaptiven Filter an jedem Zeitpunkt berechnen und sein Ausgangssignal an
%diesem Zeitpunkt
adaptive_filter_out = zeros(samples,1);
for t = 1:length(time_vector)
    % w = adapt_filter_at(w, input_signal, output_signal, mu, t);
    u = output_signal(t);
    filter_out = fir_filter(input_signal,w,t);
    e = u -filter_out;
    w = lms(w,e,input_signal,t,learnRate);
    adaptive_filter_out(t) = filter_out;
    
end
%Fehler zwischen Soll Signal und Filter Ausgang
error_out = output_signal-adaptive_filter_out;
end_error = mean(error_out(end-samples*0.1:end).^2)

%plotten
figure 
subplot(3,1,1);
hold on
plot(time_vector,adaptive_filter_out, 'x');
plot(time_vector,output_signal,'LineWidth',2);
title('Verlauf der Strecke und des adaptiven Filters')
ylabel('Ausgansamplitude')
xlabel('Zeit [s]')
hold off
legend('Filter Ausgang','Strecken Ausgang')
subplot(3,1,2);

semilogy(error_out.^2)
title(['Verlauf des Fehlers zwischen Strecke und Filter mit Lernrate mu = ',num2str(learnRate), ' und Filter taps = ', int2str(num_taps)])
ylabel('Fehler^2')
xlabel('Zeit [s]')
subplot(3,1,3);
stem(w,'.');

t_end = 0.015;
t_teststep = linspace(0,t_end,t_end*freq);
teststep = ones(t_end*freq,1);
teststep(1:t_end*freq*0.1) = 0;

[y1, tout] = lsim(strecke,teststep, t_teststep);
y2 = y1;
for i = 1:t_end*freq
    y2(i) = fir_filter(teststep, w, i);
end
figure
hold on
plot(t_teststep,y2, 'DisplayName','Adapted Filter Step Response','Marker','x');
plot(t_teststep,y1, 'DisplayName','System Step Response');
plot(t_teststep, teststep, 'DisplayName','Step')
legend
hold off
%% Adapt inverted System
mu_inv = 0.01; %learning rate for plant inversion
delta = 0; %timeshift 
n_inv = 150+delta; %length of filter c

wc = w; %copy w
c = zeros(n_inv,1); 
c_out = zeros(samples,1);
w_out = zeros(samples,1);

%Reference model
s = tf('s');
H_ref = 1/((s*8e-5+1)^5);
input_signal_ref = lsim(H_ref, input_signal,time_vector);

input_signal_delayed = zeros(samples,1);
input_signal_delayed(delta+1:end)=input_signal_ref(1:end-delta);


%precalculate output of plant filter
for k= 1:samples
    w_out(k) = fir_filter(input_signal,wc,k);
end
mu_inv = 4*min(mu_inv,1/(w_out'*w_out));
for t = delta+1:length(time_vector)
    c_out_cur = fir_filter(w_out,c,t);
    e = input_signal_delayed(t) - c_out_cur; 
    c = lms(c,e,input_signal_delayed,t,mu_inv);
    c_out(t) = c_out_cur;
end   

figure
subplot(3,1,1);
hold on
plot(input_signal_delayed,'DisplayName','Input Signal')
plot(c_out,'DisplayName','Inverted Filter Out shifted','Marker','x')
legend;
hold off

error_out_inv = input_signal_delayed-c_out;
end_error_inv = mean(error_out_inv(end-samples*0.1:end).^2)
subplot(3,1,2);
semilogy(error_out_inv.^2)
title(['Verlauf des Fehlers zwischen Eingangssignal und Ausgang des invertierten Filters mit Lernrate mu = ',num2str(learnRate), ' und Filter taps = ', int2str(num_taps)])
ylabel('Fehler^2')
xlabel('Zeit [s]')
subplot(3,1,3);
stem(c,'.')

%% Feed Forward Control
triangle = teststep;

for i = 1:length(teststep)
    triangle(i) = i*0.001;
end
teststep_filt = teststep;
triangle_filt = triangle;
cinv = c(end:-1:1);
for i = 1:length(teststep)
    teststep_filt(i) = fir_filter(teststep, c, i);
    triangle_filt(i) = fir_filter(triangle, c, i);
end
figure
hold on
plot(teststep,'DisplayName','Referenzgröße ungefiltert')
plot(teststep_filt,'DisplayName','Referenzgröße gefiltert')
legend
hold off
%Diskreter PID Regler
% PID_d = c2d(R_pid,1/freq);
% Galvo_d = c2d(S_galvo,1/freq);
% 
% 
% 
% ctrl_ff = H_fir_c*feedback(Galvo_d*PID_d,1);
% [y1, tout ] = step(ctrl_ff);
% ctrl = feedback(Galvo_d*PID_d,1);
% [y2, tout] = step(ctrl,tout);
% figure
% hold on
% plot(tout,y1);
% plot(tout,y2);
figure
subplot(2,1,1)
hold on
ctrl_out = lsim(strecke, teststep, t_teststep);
ff_out = lsim(strecke, teststep_filt, t_teststep);
plot(t_teststep, ff_out,'DisplayName','Sprungantwort mit Vorsteuerung')
plot(t_teststep, ctrl_out,'DisplayName','Sprungantwort ohne Vorsteuerung')
legend
hold off
saveVectorsAsTypFile(t_teststep, ff_out, 'FF')

subplot(2,1,2)
hold on
ctrl_out = lsim(strecke, triangle, t_teststep);
ff_out = lsim(strecke, triangle_filt, t_teststep);
plot(t_teststep, triangle,'DisplayName','Linie Führungsgröße')
plot(t_teststep, ff_out,'DisplayName','Linie mit Vorsteuerung')
plot(t_teststep, ctrl_out,'DisplayName','Linie ohne Vorsteuerung')
legend
hold off


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
% imp = ones(length(c),1);
% imp(1) = 0;
% out = zeros(length(c),1);
% for i = 1:length(c)
%     out(i) = fir_filter(imp,c,i);
% end
% figure
% plot(out)
%%
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

function saveVectorsAsTypFile(vector1, vector2, filename)
    % Prüfen, ob die Vektoren die gleiche Länge haben
    if length(vector1) ~= length(vector2)
        error('Die Vektoren müssen die gleiche Länge haben.');
    end

    % Öffnen der Datei zum Schreiben
    fileID = fopen([filename, '.typ'], 'w');

    % Schreiben des Header-Teils in die Datei
    fprintf(fileID, '#let data = (\n');

    % Schreiben der Daten in das angegebene Format
    for i = 1:length(vector1)
        fprintf(fileID, '  (%d,%d)', vector1(i), vector2(i));
        % Ein Komma und ein Zeilenumbruch hinzufügen, außer für das letzte Element
        if i < length(vector1)
            fprintf(fileID, ',\n');
        else
            fprintf(fileID, '\n');
        end
    end

    % Schließen der Klammer und Datei
    fprintf(fileID, ')\n');
    fclose(fileID);
end