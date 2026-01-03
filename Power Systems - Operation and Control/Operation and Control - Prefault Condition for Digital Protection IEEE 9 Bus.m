clear; clc;
%% 1) READ Y-BUS FROM Ybus.txt
ybusFile = 'Ybus.txt';
nBus = 1000;
fid = fopen(ybusFile,'r');
if fid<0
    error('Cannot open %s', ybusFile);
end
fmt = 'Ybus(%d , %d) = %f+ j*( %f );';
dataY = textscan(fid, fmt, 'HeaderLines',2);  % if your file requires skipping lines
fclose(fid);
rowIdx   = dataY{1};
colIdx   = dataY{2};
realPart = dataY{3};
imagPart = dataY{4};
Yvals = complex(realPart, imagPart);
Ybus_sparse = sparse(rowIdx, colIdx, Yvals, nBus, nBus);
Ybus = full(Ybus_sparse);
G = real(Ybus);
B = imag(Ybus);

%% 2) READ BUS DATA FROM data.csv
csvFile = 'data.csv';
opts = detectImportOptions(csvFile,'NumHeaderLines',1);
T = readtable(csvFile,opts);
if height(T)~= nBus
    error('data.csv has %d rows, mismatch nBus=%d',height(T),nBus);
end

userBusType = [1;2;2;3;3;3;3;3;3];  % Slack=1, PV=2, PQ=3
varNames = T.Properties.VariableNames;

if ismember('Vpu', varNames),       Vpu     = T.Vpu;       else, Vpu     = zeros(nBus,1); end
if ismember('LoadMW', varNames),    loadMW  = T.LoadMW;    else, loadMW  = zeros(nBus,1); end
if ismember('LoadMvar', varNames),  loadMvar= T.LoadMvar;  else, loadMvar= zeros(nBus,1); end
if ismember('GenMW', varNames),     genMW   = T.GenMW;      else, genMW   = zeros(nBus,1); end
if ismember('GenMvar', varNames),   genMvar = T.GenMvar;    else, genMvar = zeros(nBus,1); end
% Replace any NaNs with 0
Vpu(isnan(Vpu))=0; loadMW(isnan(loadMW))=0; loadMvar(isnan(loadMvar))=0;
genMW(isnan(genMW))=0; genMvar(isnan(genMvar))=0;
% Build busData: [type, Pd(MW), Qd(Mvar), Pg(MW), Qg(Mvar), Vset]
busData = zeros(nBus,6);
for i=1:nBus
    bType = userBusType(i);
    switch bType
        case 1  % Slack
            busData(i,:) = [1, 0,0, 0,0,  Vpu(i)];
        case 2  % PV
            busData(i,:) = [2, 0,0, genMW(i), 0, Vpu(i)];
        case 3  % PQ
            busData(i,:) = [3, loadMW(i), loadMvar(i), 0,0, 0];
        otherwise
            error('Invalid busType %d at i=%d', bType, i);
    end
end
disp('=== busData: [type, Pd, Qd, Pg, Qg, Vset] ===');
disp(busData);
%% 3) SETUP
baseMVA = 100;
PdMW = busData(:,2);
QdMW = busData(:,3);
PgMW = busData(:,4);
QgMW = busData(:,5);
Vset = busData(:,6);
Pd = PdMW/baseMVA;
Qd = QdMW/baseMVA;
Pg = PgMW/baseMVA;
Qg = QgMW/baseMVA;
busType = busData(:,1);
pvIdx = find(busType==2);
pqIdx = find(busType==3);
%% 4) INITIAL GUESS
theta = zeros(nBus,1);  % in radians
V     = ones(nBus,1);
for i=1:nBus
    switch busType(i)
        case 1
            theta(i)=0;         % Slack angle=0
            V(i)=Vset(i);
        case 2
            theta(i)=0;         % guess angle
            V(i)=Vset(i);       % known voltage setpoint
        case 3
            theta(i)=0;
            V(i)=1.0;          % typical initial guess
    end
end
maxIter = 50;
tol = 1e-6;
%% 5) NEWTON–RAPHSON ITERATIONS
for iter=1:maxIter
   %compute Pcalc, Qcalc from current states
    Pcalc = zeros(nBus,1);
    Qcalc = zeros(nBus,1);
    for k=1:nBus
        for m=1:nBus
            dth = theta(k)-theta(m);
            Pcalc(k) = Pcalc(k) + V(k)*V(m)*( G(k,m)*cos(dth) + B(k,m)*sin(dth) );
            Qcalc(k) = Qcalc(k) + V(k)*V(m)*( G(k,m)*sin(dth) - B(k,m)*cos(dth) );
        end
    end
    %  Slack bus => no  eqn
    %  PV bus => 1 eqn 
    %  PQ bus => 2 eqns
    dP = [];  pIdx = [];
    dQ = [];  qIdx = [];
    for k=1:nBus
        if busType(k)==2
            % mismatch for P only
            dP(end+1,1) = Pcalc(k) - (Pg(k) - Pd(k));
            pIdx(end+1) = k;
        elseif busType(k)==3
            % mismatch for P, Q
            dP(end+1,1) = Pcalc(k) - (Pg(k) - Pd(k));
            dQ(end+1,1) = Qcalc(k) - (Qg(k) - Qd(k));
            pIdx(end+1) = k;
            qIdx(end+1) = k;
        end
    end
    dF = [dP; dQ];
    mis = max(abs(dF));
    fprintf('Iter %2d, mismatch=%g\n', iter, mis);
    if mis<tol
        fprintf('Converged in %d iterations, mismatch=%.4g\n',iter, mis);
        break;
    end
    %---- (C) build the sub-Jacobians
    nP = length(pIdx);
    nQ = length(qIdx);
    angleVars = [pvIdx; pqIdx];  % these have unknown angles
    voltVars  = pqIdx;           % only PQ have unknown volt
    Jtt = zeros(nP,nP);
    JtV = zeros(nP,nQ);
    JVt = zeros(nQ,nP);
    JVV = zeros(nQ,nQ);
    % partial derivatives references:
    for r=1:nP
        kBus = pIdx(r);
        for c=1:nP
            mBus = angleVars(c);
            if mBus==kBus
                Jtt(r,c) = - Qcalc(kBus) - B(kBus,kBus)*V(kBus)^2;
            else
                dth = theta(kBus)-theta(mBus);
                Jtt(r,c) =   V(kBus)*V(mBus)*( G(kBus,mBus)*sin(dth) - B(kBus,mBus)*cos(dth));
            end
        end
        for c=1:nQ
            mBus = voltVars(c);
            if mBus==kBus
                JtV(r,c) = (Pcalc(kBus)/V(kBus)) + G(kBus,kBus)*V(kBus);
            else
                dth = theta(kBus)-theta(mBus);
                JtV(r,c) =  V(kBus)*( G(kBus,mBus)*cos(dth) + B(kBus,mBus)*sin(dth) );
            end
        end
    end
    % >> Fill JVt, JVV
    for r=1:nQ
        kBus = qIdx(r);
        for c=1:nP
            mBus = angleVars(c);
            if mBus==kBus
                JVt(r,c) =   Pcalc(kBus) - G(kBus,kBus)*V(kBus)^2;
            else
                dth = theta(kBus)-theta(mBus);
                JVt(r,c) = - V(kBus)*V(mBus)*( G(kBus,mBus)*cos(dth) + B(kBus,mBus)*sin(dth) );
            end
        end
        for c=1:nQ
            mBus = voltVars(c);
            if mBus==kBus
                JVV(r,c) = (Qcalc(kBus)/V(kBus)) - B(kBus,kBus)*V(kBus);
            else
                dth = theta(kBus)-theta(mBus);
                JVV(r,c) =  V(kBus)*( G(kBus,mBus)*sin(dth) - B(kBus,mBus)*cos(dth) ); % CHANGED SIGN HERE BACK TO POSITIVE
            end
        end
    end
    % Combine sub-blocks
    J = [ Jtt  JtV
          JVt  JVV ];
    %---- (D) solve
    dX = J \ (-dF);   % "Newton step" => x_new = x_old + dX
                     % so mismatch = 0 => J * dX = - dF

    dTheta = dX(1:nP);
    dVolt  = dX(nP+1 : end);
    % update angles
    for iA=1:nP
        busA = angleVars(iA);
        theta(busA) = theta(busA) + dTheta(iA);
    end
    % update volt
    for iV=1:nQ
        busV = voltVars(iV);
        V(busV) = V(busV) + dVolt(iV);
    end
end
if iter==maxIter
    fprintf('Reached maxIter=%d, final mismatch=%.4g\n', maxIter, mis);
end
%% 6) PRINT FINAL
[Pfinal, Qfinal] = deal(zeros(nBus,1), zeros(nBus,1));
for k=1:nBus
    for m=1:nBus
        dth = theta(k)-theta(m);
        Pfinal(k) = Pfinal(k) + V(k)*V(m)*( G(k,m)*cos(dth) + B(k,m)*sin(dth) );
        Qfinal(k) = Qfinal(k) + V(k)*V(m)*( G(k,m)*sin(dth) - B(k,m)*cos(dth) );
    end
end
thetaDeg = theta*180/pi;
disp('===== Final Bus Results =====');
for k=1:nBus
    fprintf('Bus %d: V=%.4f pu, angle=%.4f deg, P=%.3f MW, Q=%.3f Mvar\n',...
        k, V(k), thetaDeg(k), Pfinal(k)*baseMVA, Qfinal(k)*baseMVA);
end


%% BNRANCH FLOWS 

disp('===== Branch Flows (i->j) =====');
baseMVA = 100; % Base MVA
flowPairs = [
    4	1
    2	7
    9	3
    5	4
    6	4
    7	5
    9	6
    7	8
    8	9
];

for idx = 1:size(flowPairs,1)
    iB = flowPairs(idx,1);
    jB = flowPairs(idx,2);

    % Skip if there is no actual branch admittance between iB and jB
    if abs(G(iB,jB)) < 1e-12 && abs(B(iB,jB)) < 1e-12
        continue;
    end

    dth = theta(iB) - theta(jB);

    % branch admittance
    Gij = -G(iB, jB); 
    Bij = -B(iB, jB);


    % Power Flow from bus iB to jB
    Pij = V(iB)^2 * Gij ...
          - V(iB)*V(jB) * ( Gij*cos(dth) + Bij*sin(dth) );

    % Reactive Power Flow from bus iB to jB 
    Qij = - V(iB)^2 * Bij ...
          - V(iB)*V(jB) * ( Gij*sin(dth) - Bij*cos(dth) );

    Bc = 0;  % Initialize charging susceptance for this branch

    if ( (iB == 7 && jB == 8) || (iB == 8 && jB == 7) )
        Bc = 0.149;
    elseif ( (iB == 8 && jB == 9) || (iB == 9 && jB == 8) )
        Bc = 0.209;
    elseif ( (iB == 7 && jB == 5) || (iB == 5 && jB == 7) )
        Bc = 0.306;
    elseif ( (iB == 5 && jB == 4) || (iB == 4 && jB == 5) )
        Bc = 0.176;
    elseif ( (iB == 4 && jB == 6) || (iB == 6 && jB == 4) )
        Bc = 0.158;
    elseif ( (iB == 9 && jB == 6) || (iB == 6 && jB == 9) )
        Bc = 0.358;
    end

    if Bc > 0
        Qij_before_shunt_comp = Qij; 
        Qij = Qij - (Bc/2)*V(iB)^2;

    end

    % Convert flows from per unit to MW/Mvar
    Pij_MW  = Pij  * baseMVA;
    Qij_Mvar = Qij * baseMVA;

    fprintf('%d->%d : P = %.3f MW, Q = %.3f Mvar\n', ...
            iB, jB, Pij_MW, Qij_Mvar);
end