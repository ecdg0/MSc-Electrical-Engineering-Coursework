%Parámetros de muestre
clear all;
fs = 56; 
N = 114;  
L=N/2;
t = (0:N-1)/fs;     
dt = 1/fs;
% Parámetros de la señal
A1 = 1;             
sigma1 = -0.1;      
phi1 = 0;           
A2 = 0.25;           
sigma2 = -0.125;     
phi2 = pi/8;
y = A1 * exp(sigma1 * t) .* cos(2*pi*1*t + phi1) + A2 * exp(sigma2 * t) .* cos(2*pi*7*t + phi2);

figure;
plot(t, y, 'b-', 'LineWidth', 1.5);
xlabel('Tiempo (s)');
ylabel('y(t)');
title('Señal muestreada');

nRows = N-L; 
nCols = L;

H = zeros(nRows-1, nCols-1);
for i = 1:nRows-1
    for j = 1:nCols-1
        H(i,j) = y(i+j-1);
    end
end

% disp(H);


H1 = zeros(nRows-1, nCols-1);
for i = 2:nRows
    for j = 2:nCols
        H1(i-1,j-1) = y(i+j-2);
    end
end

% disp(H1);

[U, S, V] = svd(H);
n = 4; % modos 


Uc = U(:,1:n);
Vc = V(:,1:n);
Sc = S(1:n,1:n);

A = (Sc)^(-1/2)*transpose(Uc)*H1*Vc*(Sc)^(-1/2);

[D EIG] = eig(A);
disp(EIG);


lamda= log(EIG)/dt;
w=imag(lamda)/(2*pi);
amortiguamiento =real(lamda);

val = diag(EIG);

Z=[];
for j = 0:N-1
    for i = 1:n
            Z(j+1,i) = (val(i))^j;
    end
end

B=pinv(Z)*y'
S=abs(B)
ang=angle(B)
ang_grado=ang/pi*180
Amplitude=2*abs(B)
