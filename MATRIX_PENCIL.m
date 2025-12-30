close all;
%Parámetros de muestre
clear all;
fs = 56; 
N = 114;  
L=N/2;
t = (0:N-1)/fs;     

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

H = zeros(nRows, nCols-1);
for i = 1:nRows
    for j = 1:nCols
        H(i,j) = y(i+j-1);
    end
end

% disp(H);

[U, S, V] = svd(H);
n = 4; % modos 

V1 = V(1:L-1,1:n);
V2 = V(2:L,1:n);


Y1 = transpose(V1)*V1;

Y2 = transpose(V2)*V1;

[MatrixL,eigenvalues] = eig(inv(Y1)*Y2);
deltat = 1/fs;
lambda = log(eigenvalues)/deltat;
disp(lambda);
eigdiag = diag(eigenvalues);

for i=0:N-1
    for j=1:n
        Z(i+1,j) = (eigdiag(j))^i;
    end
end



ycol = transpose (y);

B = pinv(Z)*ycol
amp = 2*abs(B)
angulo = angle(B)