% ===== Public-release KF script  =====
clc; clear; close all; format long;

%% ---------------------------- Configuration ----------------------------
c = 299792458;                          % Speed of light [m/s]
sampleSbg = 20*10;                      % Downsample factor for SBG log
ts = 0.005 * sampleSbg;                 % Sampling interval [1s]
addpath('common')
obsTxtPath  = '..\observations\car_observations_data.txt';
trajTxtPath = '..\trajectory\CAR.txt';

%% --------------------------- Base station positions ----------------------------
baseLlh = [ ...
    40+1/60+33.22365/3600, 117+8/60+56.60000/3600, 124.508;
    40+1/60+17.20441/3600, 117+8/60+48.55499/3600, 110.562;
    40+1/60+33.29400/3600, 117+9/60+27.22288/3600, 111.040;
    40+1/60+27.25856/3600, 117+9/60+05.65661/3600, 108.825;
    40+1/60+10.04148/3600, 117+9/60+36.81872/3600, 104.366 ];
basePosUtm = geo2utm_batch(baseLlh, 5);              % 5x3 [E,N,U] in meters

%% ------------------------------ Load RTK/SBG trajectory ------------------------------
optsTraj = detectImportOptions(trajTxtPath, 'FileType', 'text');
optsTraj.DataLines = [59360, Inf];
trajTbl = readtable(trajTxtPath, optsTraj);

requiredVars = {'UTCTime','Latitude','Longitude','HeightEllipsoid'};
if ~all(ismember(requiredVars, trajTbl.Properties.VariableNames))
    error('Missing columns: %s', strjoin(setdiff(requiredVars, trajTbl.Properties.VariableNames), ', '));
end

utcTime = trajTbl.UTCTime(1:sampleSbg:end);
latDeg  = trajTbl.Latitude(1:sampleSbg:end);
lonDeg  = trajTbl.Longitude(1:sampleSbg:end);
hEllM   = trajTbl.HeightEllipsoid(1:sampleSbg:end);

t0 = utcTime(1);
trajTimeSec = seconds(utcTime(:) - t0);             % Nx1, starts from 0
trajAbsSec  = seconds(utcTime(:));                  % Nx1, absolute seconds (for alignment)

posLlh = [latDeg, lonDeg, hEllM];
posUtm = geo2utm_batch(posLlh, 5)';                 % 3xN

rxPosE = posUtm(1,:);                               % 1xN
rxPosN = posUtm(2,:);
rxPosU = posUtm(3,:);

rxVelE = [diff(rxPosE)/ts, 0];
rxVelN = [diff(rxPosN)/ts, 0];
rxVelU = [diff(rxPosU)/ts, 0];

numTraj = numel(trajTimeSec);

%% ------------------------------ True geometric ranges ------------------------------
rhoTrue = zeros(5, numTraj);
for k = 1:numTraj
    r = [rxPosE(k), rxPosN(k), rxPosU(k)];
    for i = 1:5
        rhoTrue(i,k) = norm(r - basePosUtm(i,:));
    end
end

%% ------------------------------ Load observations from TXT ------------------------------
% TXT columns:
% 1) timeSec[s]
% 2-6) codePseudorange1..5 [m]
% 7-11) pseudorangeRate1..5 [m/s]
optsObs = detectImportOptions(obsTxtPath, 'FileType', 'text', 'Delimiter', '\t');
obsTbl = readtable(obsTxtPath, optsObs);

obsMat = table2array(obsTbl);
obsTimeSec = obsMat(:,1);               % Nx1 (already shifted at export stage)
obsCodeP   = obsMat(:,2:6);             % Nx5
obsRhoRate = obsMat(:,7:11);            % Nx5
obsCodeP5xN   = obsCodeP.';             % 5xN
obsRhoRate5xN = obsRhoRate.';           % 5xN

%% ------------------------------ Time alignment (trajectory -> observation) ------------------------------
maxObsTime = max(obsTimeSec);
validIdx = trajAbsSec <= maxObsTime;
numEpoch = sum(validIdx);

zMatched = zeros(11, numEpoch);         % 11 = 5 PR + 5 PRR + 1 height
matchIdx = zeros(1, numEpoch);

for k = 1:numEpoch
    [~, matchIdx(k)] = min(abs(obsTimeSec - trajAbsSec(k)));
    zMatched(1:5,  k) = obsCodeP5xN(:, matchIdx(k));
    zMatched(6:10, k) = obsRhoRate5xN(:, matchIdx(k));
end

measureTimeSec = trajAbsSec(1:numEpoch);            % 1..K 
measureCodeP   = zMatched(1:5,  :);                 % 5xK
measureRhoRate = zMatched(6:10, :);                 % 5xK

rhoTrueAligned = rhoTrue(:, 1:numEpoch);            % 5xK
rhoDiff = measureCodeP - rhoTrueAligned;            % 5xK

%% ------------------------------ Initial clock bias/drift ------------------------------
kInit = 1;

% True state at the aligned trajectory epoch kInit
posInit = [rxPosE(kInit), rxPosN(kInit), rxPosU(kInit)];        % [m]
velInit = [rxVelE(kInit), rxVelN(kInit), rxVelU(kInit)];        % [m/s]
distTrue0      = zeros(5,1);   % [m]
rangeRateTrue0 = zeros(5,1);   % [m/s]

for bs = 1:5
    delta = posInit - basePosUtm(bs,:);      % user - base
    d = norm(delta);
    distTrue0(bs) = d;

    if d > 0
        u = delta / d;
    else
        u = [0,0,0];
    end

    rangeRateTrue0(bs) = dot(velInit, u);
end

measRho0     = measureCodeP(:, kInit);       % 5x1 [m]
measRhoRate0 = measureRhoRate(:, kInit);     % 5x1 [m/s]

initClockBiasM    = measRho0     - distTrue0;        % 5x1 [m]
initClockDriftMps = measRhoRate0 - rangeRateTrue0;   % 5x1 [m/s]
%% ------------------------------ Kalman filter model ------------------------------
paraRub.h0   = 9.9e-21; paraRub.h_2 = 1.2e-22;
paraTCXO.h0  = 9.7e-21; paraTCXO.h_2 = 7.2e-20;
paraRub.h0   = 9.7/45*1e-21; paraRub.h_2 = 1/45*1e-22;
paraTCXO.h0  = 9.5/45*1e-21; paraTCXO.h_2 = 7/45*1e-20;

Fpv = [1 0 0 ts 0  0;
       0 1 0 0 ts 0;
       0 0 1 0  0 ts;
       0 0 0 1 0  0;
       0 0 0 0 1  0;
       0 0 0 0 0  1];
Fclk = [1 ts; 0 1];
F = blkdiag(Fpv, Fclk, Fclk, Fclk, Fclk, Fclk);

Sf_s = 2*pi^2*paraRub.h_2;  St_s = paraRub.h0/2;
Sf_r = 2*pi^2*paraTCXO.h_2; St_r = paraTCXO.h0/2;

Q_s = [St_s*ts + Sf_s*ts^3/3, Sf_s*ts^2/2;
       Sf_s*ts^2/2,           Sf_s*ts];
Q_r = [St_r*ts + Sf_r*ts^3/3, Sf_r*ts^2/2;
       Sf_r*ts^2/2,           Sf_r*ts];

Qclk = blkdiag(Q_s, Q_s, Q_s, Q_s, Q_s) + repmat(Q_r, 5, 5);
Qclk = Qclk * c^2;

qx = 2^2;  qy = 2^2;  qz = 1^2;

Qpv = [qx*ts^3/3 0 0 qx*ts^2/2 0 0;
       0 qy*ts^3/3 0 0 qy*ts^2/2 0;
       0 0 qz*ts^3/3 0 0 qz*ts^2/2;
       qx*ts^2/2 0 0 qx*ts 0 0;
       0 qy*ts^2/2 0 0 qy*ts 0;
       0 0 qz*ts^2/2 0 0 qz*ts];

Q = blkdiag(Qpv, Qclk);

% Measurement noise: 5 PR + 5 PRR + 1 height
R = diag([0.5^2*ones(1,5), 0.1^2*ones(1,5), 0.2^2]);

%% ------------------------------ State initialization ------------------------------
% State: [pos(3), vel(3), (bias,drift)x5] => 16 states
xEst = zeros(16, numEpoch);
xEst(:,1) = [rxPosE(1); rxPosN(1); rxPosU(1); rxVelE(1); rxVelN(1); rxVelU(1); ...
             reshape([initClockBiasM.'; initClockDriftMps.'], 10, 1)];

%% ------------------------------ ML-based covariance initialization ------------------------------
k0 = 1;
k1 = k0 + 2;

SigmaIni = blkdiag(diag(0.5^2*ones(1,3)), diag(0.5^2*ones(1,3)), ...
                   diag(0.3^2*ones(1,5)), diag(0.3^2*ones(1,5)));

Tini = trajAbsSec(k1) - trajAbsSec(k0);

r0 = [rxPosE(k0); rxPosN(k0); rxPosU(k0)];
r1 = [rxPosE(k1); rxPosN(k1); rxPosU(k1)];

z0 = measureCodeP(:, k0);
z1 = measureCodeP(:, k1);

numBs = 5;
dr0 = zeros(numBs,1); dr1 = zeros(numBs,1);
hr0 = zeros(numBs,3); hr1 = zeros(numBs,3);

for i = 1:numBs
    s = basePosUtm(i,:).';
    dr0(i) = norm(r0 - s);
    dr1(i) = norm(r1 - s);
    hr0(i,:) = ((r0 - s)/dr0(i)).';
    hr1(i,:) = ((r1 - s)/dr1(i)).';
end

cdtHat    = z1 - dr1;
cdtDotHat = (z1 - z0 + dr0 - dr1) / Tini;

Aini = zeros(16, 6 + 2*numBs);
Aini(1:3, 1:3) = eye(3);
Aini(4:6, 1:3) = (1/Tini)*eye(3);
Aini(4:6, 4:6) = (-1/Tini)*eye(3);

for i = 1:numBs
    biasIdx  = 6 + 2*i - 1;
    driftIdx = 6 + 2*i;

    Aini(biasIdx, 1:3)   = -hr1(i,:);
    Aini(biasIdx, 6+i)   = 1;

    Aini(driftIdx, 1:3)          = (-hr1(i,:))/Tini;
    Aini(driftIdx, 4:6)          = ( hr0(i,:))/Tini;
    Aini(driftIdx, 6+i)          =  1/Tini;
    Aini(driftIdx, 6+numBs+i)    = -1/Tini;
end

P = Aini * SigmaIni * Aini.';
P = (P + P.')/2;
[eV,eD] = eig(P);
eD = max(eD, 1e-6*eye(size(eD)));
P = eV*eD*eV.';

%% ------------------------------ Filter loop ------------------------------
sigmaPosHist = zeros(3, numEpoch);
sigmaPosHist(:,1) = [sqrt(P(1,1)); sqrt(P(2,2)); sqrt(P(3,3))];

for k = 2:numEpoch
    xPrior = F * xEst(:,k-1);
    PPrior = F * P * F' + Q;

    H = zeros(11,16);
    h = zeros(11,1);
    z = zeros(11,1);

    rPrior = xPrior(1:3);
    vPrior = xPrior(4:6);

    for i = 1:5
        delta = rPrior - basePosUtm(i,:)';
        dist  = norm(delta);
        u     = delta / dist;

        % Pseudorange
        H(i,1:3) = u';
        H(i,6+i*2-1) = 1;
        h(i) = dist + xPrior(6+i*2-1);

        % Pseudorange-rate
        vDotDelta = dot(vPrior, delta);
        H(i+5,1:3) = vPrior'/dist - (vDotDelta/(dist^3))*delta';
        H(i+5,4:6) = u';
        H(i+5,6+i*2) = 1;
        h(i+5) = dot(vPrior, u) + xPrior(6+i*2);
    end

    % Height constraint (SBG/RTK altitude)
    H(11,3) = 1;
    h(11)   = xPrior(3);

    z(1:5)   = measureCodeP(:,k);
    z(6:10)  = measureRhoRate(:,k);
    z(11)    = rxPosU(k);

    S = H * PPrior * H' + R;
    K = PPrior * H' / S;

    xEst(:,k) = xPrior + K * (z - h);
    P = (eye(16) - K * H) * PPrior;

    sigmaPosHist(:,k) = [sqrt(P(1,1)); sqrt(P(2,2)); sqrt(P(3,3))];
end

%% ------------------------------ 3-sigma envelope plot ------------------------------
timeMin = trajTimeSec(1:numEpoch) / 60;

posErrE = xEst(1,:).' - rxPosE(1:numEpoch).';
posErrN = xEst(2,:).' - rxPosN(1:numEpoch).';
posErrU = xEst(3,:).' - rxPosU(1:numEpoch).';

sigE = sigmaPosHist(1,:).';
sigN = sigmaPosHist(2,:).';
sigU = sigmaPosHist(3,:).';

b3E = 3*sigE; b3N = 3*sigN; b3U = 3*sigU;

ignoreSteps = 200;
if numEpoch > ignoreSteps
    scaleIdx = ignoreSteps:numEpoch;
else
    scaleIdx = 1:numEpoch;
end

limE = max(abs([b3E(scaleIdx); posErrE(scaleIdx)]));
limN = max(abs([b3N(scaleIdx); posErrN(scaleIdx)]));
limU = max(abs([b3U(scaleIdx); posErrU(scaleIdx)]));
limE = max(limE, 0.1);
limN = max(limN, 0.1);
limU = max(limU, 0.1);

fillColor  = [1, 0.85, 0.85];
lineColor  = [0, 0.4470, 0.7410];
boundColor = [0.8500, 0.3250, 0.0980];

figure('Units','pixels','Position',[150,150,800,700],'Color','w');

tFill = [timeMin; flipud(timeMin)];
subplot(3,1,1); hold on; box on; grid on;
fill(tFill, [b3E; flipud(-b3E)], fillColor, 'EdgeColor','none', 'FaceAlpha',0.5);
plot(timeMin,  b3E, '--', 'Color', boundColor, 'LineWidth', 1);
plot(timeMin, -b3E, '--', 'Color', boundColor, 'LineWidth', 1);
plot(timeMin, posErrE, '-', 'Color', lineColor, 'LineWidth', 1.2);
ylabel('East Error (m)'); set(gca,'FontSize',10);
ylim([-limE*1.5, limE*1.5]); xlim([0, timeMin(end)]);

subplot(3,1,2); hold on; box on; grid on;
fill(tFill, [b3N; flipud(-b3N)], fillColor, 'EdgeColor','none', 'FaceAlpha',0.5);
plot(timeMin,  b3N, '--', 'Color', boundColor, 'LineWidth', 1);
plot(timeMin, -b3N, '--', 'Color', boundColor, 'LineWidth', 1);
plot(timeMin, posErrN, '-', 'Color', lineColor, 'LineWidth', 1.2);
ylabel('North Error (m)'); set(gca,'FontSize',10);
ylim([-limN*1.5, limN*1.5]); xlim([0, timeMin(end)]);

subplot(3,1,3); hold on; box on; grid on;
fill(tFill, [b3U; flipud(-b3U)], fillColor, 'EdgeColor','none', 'FaceAlpha',0.5);
plot(timeMin,  b3U, '--', 'Color', boundColor, 'LineWidth', 1);
plot(timeMin, -b3U, '--', 'Color', boundColor, 'LineWidth', 1);
plot(timeMin, posErrU, '-', 'Color', lineColor, 'LineWidth', 1.2);
ylabel('Up Error (m)'); xlabel('Time (min)'); set(gca,'FontSize',10);
ylim([-limU*1.5, limU*1.5]); xlim([0, timeMin(end)]);
