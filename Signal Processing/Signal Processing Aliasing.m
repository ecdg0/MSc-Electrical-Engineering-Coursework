% Universidad de Guadalajara - Master of Science in Electrical Engineering
% Erick Christopher Davalos Gonzalez

% Digital Signal Processing
% Project 1. Aliasing 

%Consider the following signal:
% xa(t) = sin(2*pi*Fo*t)

% The above signal can be described in it's sampling form:
% x(n) = xa(nT) = sin(2*pi*(Fo/Fs)*n) where Fs = 1/T

% a)
% Sampling the signal in the interval 0<= n <=99, Fs = 5kHz Fo = [0.5,
% 2,3,4.5] kHz

Fs = 5;  % Sampling Frequency in KHz
F0_values = [0.5, 2, 3, 4.5];  % Fundamental Frequencies
n = 0:99;  % Samples
t = n / Fs;  % Time for each sample
t0 = linspace(0,(max(t))/4,1000);


% Creating a figure to contain subpolts
figure;
for i = 1:length(F0_values)
    F0 = F0_values(i);
    y = sin(2*pi*(F0/Fs)*n);  % Sampled signal
    subplot(2, 2, i);  % 2x2 plot
    % Graph
    plot(t, y, 'k', 'LineWidth', 1.5);  
    hold on;  
    % Point and line for each sample
    stem(t, y, '--r', 'LineWidth', 1, 'Marker', '.','MarkerSize', 15,'Color','r');
    xlim([0 5]);

    title(['F_0 = ', num2str(F0), ' kHz']);
    xlabel('Time (ms)');
    ylabel('Amplitude');
    hold off;
    
end
sgtitle('Sampled Signal and Samples');  % Title

% Creating other figure to contain plot original signal vs sampled signal.
figure;
for i = 1:length(F0_values)

    F0 = F0_values(i);
    y0 = sin(2*pi*(F0)*t0);  % Original Signal
    y = sin(2*pi*(F0/Fs)*n);  % Sampled signal

    subplot(2, 2, i);  % 2x2 plot

    % Graph
    plot(t0, y0, 'k', 'LineWidth', 1.5);  
    hold on; 
    plot(t,y, 'r', 'LineWidth', 1.5);  
    stem(t, y,'LineStyle','none','Marker', '.','MarkerSize', 15,'Color','r');
    xlim([0 4])

    title(['F_0 = ', num2str(F0), ' kHz']);
    xlabel('Time (ms)');
    ylabel('Amplitude');
    hold off;

 end
 sgtitle('Original Signal vs Sampled Signal');  % Title


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% b)
Fs = 50;  % Sampling Frequency in KHz
F0 = 2;   % Fundamental Frequency in Khz
n = 0:99;  % Samples
t = n / Fs;  % Time


% Graph
x = sin(2*pi*(F0/Fs)*n);
figure;

subplot(2, 2, 1);  
y0 = sin(2*pi*(F0)*t0); %original signal
plot(t, x, '.-r', 'LineWidth', 3);
hold on;
plot(t0, y0, 'k', 'LineWidth', 2); %plotting original signal
xlim([0 1])
stem(t, x,'--r', 'LineWidth', 1, 'Marker', '.','MarkerSize', 15,'Color','r');
title(['x(n) Signal F_0 = ',num2str(F0),'kHz y F_s =',num2str(Fs),' kHz']);

xlabel('Time (ms)');
ylabel('Amplitude');
hold off;

% taking even samples
y = x(1:2:end);  % taking even samples
t_y = t(1:2:end);  % time of even samples

% Graph y(n)
subplot(2, 2, 2);  
plot(t_y, y, 'k', 'LineWidth', 1.5);
hold on;
stem(t_y, y,'--r', 'LineWidth', 1, 'Marker', '.','MarkerSize', 15,'Color','r');
xlim([0 1])
title('y(n) obtained from the even samples of x(n)');
xlabel('Time (ms)');
ylabel('Amplitude');
hold off;
sgtitle('Analysis of x(n) even samples');