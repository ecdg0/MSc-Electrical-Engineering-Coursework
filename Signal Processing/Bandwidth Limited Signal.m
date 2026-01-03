%Erick Davalos
%Universidad de Guadalajara
DC = 0.2;
FsT = 100e3;
dtT = 1/FsT;
Tobs = (1/40);
NT = Tobs/dtT;
nT = 0:NT-1;
tT = nT*dtT;
dfT = FsT/NT;

xT = DC + 0.4*cos(2*pi*120*tT)+ 0.3*cos(2*pi*240*tT) +0.2*cos(2*pi*600*tT)+0.3*cos(2*pi*1200*tT);
%xT = DC + 0.4*cos(2*pi*120*tT)+ 0.3*cos(2*pi*245*tT) +0.2*cos(2*pi*600*tT)+0.3*cos(2*pi*1205*tT);
%%
    
Fs = 3.2e3; %frecuencia de muestreo 
dt = 1/Fs;
N = Tobs/dt;
n = 0:N-1;
t = n*dt;
df = Fs/N;
freq = (0:N-1)*Fs/N;

x = DC + 0.4*cos(2*pi*120*t)+ 0.3*cos(2*pi*240*t) +0.2*cos(2*pi*600*t)+0.3*cos(2*pi*1200*t);
%x = DC + 0.4*cos(2*pi*120*t)+ 0.3*cos(2*pi*242*t) +0.2*cos(2*pi*600*t)+0.3*cos(2*pi*1205*t);


figure(8)
set(gcf,'defaultAxesFontName','Times','DefaultAxesFontSize',16);
plot(tT,xT,'k','LineWidth',2);
set(gcf,'pos',[80 80 500 350]), grid on;
xlabel('Time[s]'),ylabel('Amplitude');
legend('Analogica');

figure(1)
set(gcf,'defaultAxesFontName','Times','DefaultAxesFontSize',16);
plot(tT,xT,'k',t,x,'-.r','LineWidth',2);
set(gcf,'pos',[100 100 500 350]), grid on;
xlabel('Time[s]'),ylabel('Amplitude');
legend('Analogica','Sampled');

figure(2)
set(gcf,'defaultAxesFontName','Times','DefaultAxesFontSize',16);
plot(n,x,'-.r','LineWidth',2);
set(gcf,'pos',[120 120 500 350]), grid on;
xlabel('Time [n]'),ylabel('Amplitude');
legend('Sampled');

tic
X = zeros(1,N);
for k = 1:N
    for m = 1:N
        X(k) = X(k) + x(m)*exp(-1j*2*pi*(k-1)*(m-1)/N);
    end
end
toc

%espectro muextreado, eje en "n"
figure(3)
set(gcf,'defaultAxesFontName','Times','DefaultAxesFontSize',16);
stem(n,abs(X),'fill','--','MarkerFaceColor','k','linewidth',2);
set(gcf,'pos',[140 140 500 350]), grid on;
xlabel('X(k)'),ylabel('Magnitud');

figure(9)
set(gcf,'defaultAxesFontName','Times','DefaultAxesFontSize',16);
stem(freq,abs(X),'fill','--','MarkerFaceColor','k','linewidth',2);
set(gcf,'pos',[240 240 500 350]), grid on;
xlabel('Frequency Hz'),ylabel('Magnitud');

%espectro continuo, eje en "n"
figure(4)
set(gcf,'defaultAxesFontName','Times','DefaultAxesFontSize',16);
plot(n,abs(X),'k','linewidth',2);
set(gcf,'pos',[160 160 500 350]), grid on;
xlabel('X(k)'),ylabel('Magnitud');


% ESPECTRo continuo y eje en "Frecuencia"
figure(5)
set(gcf,'defaultAxesFontName','Times','DefaultAxesFontSize',16);
plot(freq,abs(X),'k','linewidth',2);
set(gcf,'pos',[180 180 500 350]), grid on;
xlabel('Frequency [Hz]'),ylabel('Magnitud');

xf = zeros(1,N);
for m = 1:N
    for k = 1:N
        xf(m) = xf(m) + X(k)*exp(1j*2*pi*(k-1)*(m-1)/N);
    end
    xf(m) = xf(m)/N;
end

figure(6)
set(gcf,'defaultAxesFontName','Times','DefaultAxesFontSize',16);
plot(n,x,'-.r',n,real(xf),'--b','LineWidth',2);
set(gcf,'pos',[200 200 500 350]), grid on;
xlabel('Time [n]'),ylabel('Amplitude');
legend('Sampled','Reconstructed');

figure(7)
set(gcf,'defaultAxesFontName','Times','DefaultAxesFontSize',16);
plot(tT,xT,'k',t,x,'-.r',t,real(xf),'--b','LineWidth',2);
set(gcf,'pos',[220 220 500 350]), grid on;
xlabel('Time [n]'),ylabel('Amplitude');
legend('Analog','Sampled','Reconstructed');


% Ideal Interpolator
t_interp = tT; % Time points for the interpolated signal
x_interp = zeros(size(t_interp)); % Initialize the interpolated signal

for i = 1:length(t_interp)
    % Ideal reconstruction using the sinc interpolation
    x_interp(i) = sum(real(xf) .* sinc(Fs * (t_interp(i) - t)));
end

% Plot the results
figure(10)
set(gcf,'defaultAxesFontName','Times','DefaultAxesFontSize',16);
plot(tT, xT, 'k', 'LineWidth', 2); % Original 
hold on;
% plot(t, real(xf), 'b', 'LineWidth', 2); % Reconstructed DFT
plot(t_interp, x_interp, '-.r', 'LineWidth', 2); % interpolated signal
hold off;
set(gcf,'pos',[250 250 500 350]), grid on;
xlabel('Time [s]');
ylabel('Amplitude');
%legend('Original Signal', 'Reconstructed from DFT', 'Interpolated Signal');

%% fft

tic
xfft=fft(x);
toc

figure(11)
set(gcf,'defaultAxesFontName','Times','DefaultAxesFontSize',16);
plot(freq,abs(X),'k',freq,abs(xfft),'--r','LineWidth',2);
set(gcf,'pos',[60 60 500 350]), grid on;
xlabel('Frequency[Hz]'),ylabel('Magnitude');
