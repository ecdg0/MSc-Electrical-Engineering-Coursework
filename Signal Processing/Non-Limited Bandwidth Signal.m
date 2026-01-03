clc;
clear;
close all;

% Parameters
 
FsT = 100e3;         
dtT = 1 / FsT;
Tobs = 2 * pi;       % Observation time 
NT = Tobs / dtT;
nT = 0:NT-1;
%tT =linspace(0 , Tobs , NT);
tT = nT * dtT;

% Generate the square wave (analog signal)
xT = zeros(1, length(tT)); % Initialize xT
for k = 1:length(tT)
    if mod(tT(k), 2 * pi) < pi
        xT(k) = 1;
    else
        xT(k) = 0;
    end
end

% Sampling parameters
Fs = 10 / pi;         % 10 samples
dt = 1 / Fs;
N = round(Tobs / dt); 
n = 0:N-1;
t = (n + 0.5)*dt;

% Generate the sampled square wave
x = zeros(1, N); % Initialize sampled signal
for k = 1:N
    if mod(t(k), 2 * pi) < pi
        x(k) = 1;
    else
        x(k) = 0;
    end
end

% Plot the analog and sampled signals
figure(1)
set(gcf, 'defaultAxesFontName', 'Times', 'DefaultAxesFontSize', 16);
% plot(tT, xT, 'k', 'LineWidth', 2); % Analog signal
% hold on;
xlim([0,Tobs]);
stem(t, x, 'r-.', 'LineWidth', 2); % Sampled signal
% hold off;
set(gcf, 'pos', [100 100 500 350]), grid on;
xlabel('Time [s]'), ylabel('Amplitude');


% Compute the DFT
X = zeros(1, N);
tic
for k = 1:N
    for m = 1:N
        X(k) = X(k) + x(m) * exp(-1j * 2 * pi * (k - 1) * (m - 1) / N);
    end
end
toc

% Plot the magnitude spectrum (discrete)
figure(2)
stem(n, abs(X), 'fill', '--', 'MarkerFaceColor', 'k', 'LineWidth', 2);
set(gcf, 'pos', [140 140 500 350]), grid on;
xlabel('Time [n]'), ylabel('Magnitude');

% Plot the magnitude spectrum vs. frequency
freq = (0:N-1) * Fs / N;

figure(3)
stem(freq, abs(X), 'fill', '--', 'MarkerFaceColor', 'k', 'LineWidth', 2);
set(gcf, 'pos', [240 240 500 350]), grid on;
xlabel('Frequency [Hz]'), ylabel('Magnitude');

figure(4)
set(gcf,'defaultAxesFontName','Times','DefaultAxesFontSize',16);
plot(freq,abs(X),'k','linewidth',2);
set(gcf,'pos',[160 160 500 350]), grid on;
xlabel('Frequency [Hz]'),ylabel('Magnitude');

% Reconstruct the signal from its DFT

xf = zeros(1, N);
for m = 1:N
    for k = 1:N
        xf(m) = xf(m) + X(k) * exp(1j * 2 * pi * (k - 1) * (m - 1) / N);
    end
    xf(m) = xf(m) / N;
end

% Plot the sampled and reconstructed signals
figure(5)
set(gcf,'defaultAxesFontName','Times','DefaultAxesFontSize',16);
plot(n, x, 'r-.', 'LineWidth', 2);
hold on;
plot(n, real(xf), '--b', 'LineWidth', 2);
hold off;
set(gcf, 'pos', [200 200 500 350]), grid on;
xlabel('Time [n]'), ylabel('Amplitude');
legend('Sampled', 'Reconstructed');


% Plot the original and interpolated signal
x_interp = zeros(size(tT)); 

for i = 1:length(x_interp)
    x_interp(i) = sum(x .* sinc(Fs * (tT(i) - t)));
end

figure(6)
set(gcf,'defaultAxesFontName','Times','DefaultAxesFontSize',16);
plot(tT, xT, 'k', 'LineWidth', 2); % Original signal
hold on;

plot(tT, x_interp, 'b', 'LineWidth', 2); % Interpolated signal
hold off;
xlim([0,Tobs]);
ylim([-0.3,1.3])
set(gcf, 'pos', [250 250 500 350]), grid on;
xlabel('Time [s]');
ylabel('Amplitude');
legend('Original Signal', 'Interpolated Signal');

%%%%

tic
xfft=fft(x);
toc

figure(7)
set(gcf,'defaultAxesFontName','Times','DefaultAxesFontSize',16);
plot(freq,abs(X),'k',freq,abs(xfft),'--r','LineWidth',2);
set(gcf,'pos',[60 60 500 350]), grid on;
xlabel('Frequency[Hz]'),ylabel('Magnitude');

%%%

% Fourier Series coef from the DFT
c_n = X / N; % Scale DFT coefficients to get Fourier Series coef
% Fourier Series
x_fourier = zeros(1, length(tT));
omega_0 = 2 * pi / Tobs; % Fundamental frequency
x_fourier = c_n(1); % DC component


x_fourier = c_n(1); % DC
for k = 2:N/2
    x_fourier = x_fourier + 2*real(c_n(k) * exp(1j * (k - 1) * omega_0 * (tT-0.5*dt)));
end


figure(10);
stem(0:N-1, abs(c_n), 'fill', '--', 'MarkerFaceColor', 'k', 'LineWidth', 2);
xlabel('Frequency Index (k)');
ylabel('|c_n|');
grid on;

% Plot 
figure(8)
set(gcf, 'defaultAxesFontName', 'Times', 'DefaultAxesFontSize', 16);
plot(tT, xT, 'k', 'LineWidth', 2); % Original
hold on;
plot(tT, x_fourier, 'r', 'LineWidth', 2);
hold off;
xlim([0,Tobs]);
ylim([-0.3,1.3])
set(gcf, 'pos', [300 300 500 350]), grid on;
xlabel('Time [s]');
ylabel('Amplitude');
legend('Original Signal', 'Fourier Series Reconstruction');

