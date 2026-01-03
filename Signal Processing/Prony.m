% Parámetros de muestreo
fs = 56;          
N = 114;           
t = (0:N-1)/fs;    

% Parámetros de la señal
A1 = 1;            
sigma1 = -0.1;     
phi1 = 0;          
A2 = 0.25;          
sigma2 = -0.125;    
phi2 = pi/8;
y = A1 * exp(sigma1 * t) .* cos(2*pi*1*t + phi1) + A2 * exp(sigma2 * t) .* cos(2*pi*7*t + phi2);
n=4;
figure;
plot(t, y, 'b-', 'LineWidth', 1.5);
xlabel('Tiempo (s)');
ylabel('y(t)');
title('Señal muestreada');

nRows = 109;
nCols = 4;

T = zeros(nRows, nCols);
for i = 1:nRows
    for j = 1:nCols
        T(i,j) = y((nCols-j)+i);
    end
end

disp(T);
for i=1:(length(T))
    b(i,1) = y(1,i+nCols);
end
disp(b);

a = inv(transpose(T)*T)*transpose(T)*b;
disp(a);
disp(size(a));

z = roots([1 -a']);
fprintf('Raices\n');
disp (z);
deltat = 1/fs;

lambda = log(z)/deltat;
fprintf('valores propios\n');
disp(lambda);

ncols = 114;
nfil = 114;


Z2= [];
for j = 0:N-1
    for i = 1:n
        Z2(j+1,i) = (z(i))^j;
    end
end

B = pinv(Z2)*y'

Babs = abs(B)
Bang = angle(B)

Amplitude = 2* Babs
