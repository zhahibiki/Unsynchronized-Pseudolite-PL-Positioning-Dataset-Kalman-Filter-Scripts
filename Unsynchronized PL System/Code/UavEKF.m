clc; clear; close all; format long;
addpath('common');

%% ====================== 0) PATHS ======================
dataDir   = '..\observations\';
fileBase  = fullfile(dataDir, 'uav_base_stations.txt');         % LLH in txt (lat,lon,h)
fileObs   = fullfile(dataDir, 'uav_observations.txt');  

fileGT_Raw ='..\trajectory\UAV.txt';  

%% ====================== 1) CONST ======================
lightSpeed = 299792458;

%% ====================== 2) LOAD BASE (TXT: LLH) ======================
fprintf('Loading base stations (TXT: LLH)...\n');
if ~isfile(fileBase), error('Base file not found: %s', fileBase); end
baseTbl = readtable(fileBase, 'Delimiter', '\t');
baseLLH = [baseTbl{:,2}, baseTbl{:,3}, baseTbl{:,4}];%[lat(deg), lon(deg), h(m)]
% Convert once to UTM (same as your main script)
basePos = geo2utm_batch(baseLLH, 5);      % 5x3
numTx = size(basePos,1);

%% ====================== 3) LOAD OBS (PUBLIC TXT) ======================
fprintf('Loading observations (TXT: code PR)...\n');
obsTbl = readtable(fileObs, 'Delimiter', '\t');

obsTime_abs = obsTbl{:,1};                   % [s] 
measCode    = obsTbl{:,2:(1+numTx)};         % Nx5
numEpochs = length(obsTime_abs);
dt_real = mean(diff(obsTime_abs));

fprintf('Obs epochs: %d, dt≈%.4f s\n', numEpochs, dt_real);

%% ====================== 4) LOAD RTK TRAJ (SAME AS MAIN SCRIPT) ======================
fprintf('Loading RTK ground truth...\n');
opts = detectImportOptions(fileGT_Raw, 'FileType', 'text');
opts.DataLines = [77600, inf];
data = readtable(fileGT_Raw, opts);

requiredVars = {'UTCTime','Latitude','Longitude','HeightEllipsoid', ...
                'EastVelocity','NorthVelocity','DownVelocity'};

sampleSbg = 20*10;  % keep original
timeind = data.UTCTime(1:sampleSbg:end);
lat  = data.Latitude(1:sampleSbg:end);
lon  = data.Longitude(1:sampleSbg:end);
h_ell = data.HeightEllipsoid(1:sampleSbg:end);

rx.VelX = data.EastVelocity(1:sampleSbg:end);
rx.VelY = data.NorthVelocity(1:sampleSbg:end);
rx.VelZ = -data.DownVelocity(1:sampleSbg:end);

start_time1 = timeind(1);
trueT_rel   = seconds(timeind(:) - start_time1);  
times_abs   = seconds(timeind(:));                 

mat_pos_lbh = [lat, lon, h_ell];
mat_pos = geo2utm_batch(mat_pos_lbh, 5)';       
rx.PosX = mat_pos(1,:);
rx.PosY = mat_pos(2,:);
rx.PosZ = mat_pos(3,:);

sampleLen = numel(rx.PosX);

%% ====================== 5) ALIGN GT TO OBS BY ABS TIME ======================

idxGT = zeros(numEpochs,1);
for k = 1:numEpochs
    [~, idxGT(k)] = min(abs(times_abs - obsTime_abs(k)));
end

truePos = [rx.PosX(idxGT)', rx.PosY(idxGT)', rx.PosZ(idxGT)'];
trueVel = [rx.VelX(idxGT),  rx.VelY(idxGT),  rx.VelZ(idxGT)];

%% ====================== 6) INIT (BIAS/DRIFT FROM FIRST 2 OBS) ======================
fprintf('Initializing clock states from TXT epoch1-2...\n');

pos_init = truePos(1,:);
bias_init = zeros(numTx,1);
for i = 1:numTx
    dist = norm(pos_init - basePos(i,:));
    bias_init(i) = measCode(1,i) - dist;
end

dt_init = obsTime_abs(2) - obsTime_abs(1);
drift_init = (measCode(2,:)' - measCode(1,:)') / dt_init;

%% ====================== 7) EKF SETUP  ======================
% State: [pos(3), vel(3), bias(5), drift(5)] => 16

% ---- F ----
F_PV  = [1 0 0 dt_real 0 0;
         0 1 0 0 dt_real 0;
         0 0 1 0 0 dt_real;
         0 0 0 1 0 0;
         0 0 0 0 1 0;
         0 0 0 0 0 1];
F_clk = [1 dt_real; 0 1];
F = blkdiag(F_PV, F_clk, F_clk, F_clk, F_clk, F_clk);

% ---- Q (same values as your main script) ----
paraRub.h0  = 9.7/45*1e-21;  paraRub.h_2  = 1/45*1e-22;
paraTCXO.h0 = 9.5/45*1e-21;  paraTCXO.h_2 = 7/45*1e-20;

Sf_s = 2*pi^2*paraRub.h_2;   St_s = paraRub.h0/2;
Sf_r = 2*pi^2*paraTCXO.h_2;  St_r = paraTCXO.h0/2;

Q_s = [St_s*dt_real + Sf_s*dt_real^3/3, Sf_s*dt_real^2/2;
       Sf_s*dt_real^2/2,               Sf_s*dt_real];

Q_r = [St_r*dt_real + Sf_r*dt_real^3/3, Sf_r*dt_real^2/2;
       Sf_r*dt_real^2/2,               Sf_r*dt_real];

Q_clk = blkdiag(Q_s, Q_s, Q_s, Q_s, Q_s) + repmat(Q_r, 5, 5);
Q_clk = Q_clk * lightSpeed^2;

qx = 2.2^2; qy = 2.2^2; qz = 0.5^2;
Q_pv = [qx*dt_real^3/3 0 0 qx*dt_real^2/2 0 0;
        0 qy*dt_real^3/3 0 0 qy*dt_real^2/2 0;
        0 0 qz*dt_real^3/3 0 0 qz*dt_real^2/2;
        qx*dt_real^2/2 0 0 qx*dt_real 0 0;
        0 qy*dt_real^2/2 0 0 qy*dt_real 0;
        0 0 qz*dt_real^2/2 0 0 qz*dt_real];

Q = blkdiag(Q_pv, Q_clk);
R = diag([0.18^2*ones(1,5), 0.2^2]);

%% ====================== 8) ML-STYLE INITIAL COV (COPY YOUR LOGIC) ======================

k0 = 1;
k1 = k0 + 2;
T_ini = obsTime_abs(k1) - obsTime_abs(k0);

zrr0 = truePos(k0,:)';
zrr1 = truePos(k1,:)';

z0 = measCode(k0,:)';
z1 = measCode(k1,:)';

N_ml = numTx;
dr1 = zeros(N_ml,1); dr0 = zeros(N_ml,1);
hr1 = zeros(N_ml,3); hr0 = zeros(N_ml,3);
for n = 1:N_ml
    rsn = basePos(n,:).';
    dr1(n) = norm(zrr1 - rsn);
    dr0(n) = norm(zrr0 - rsn);
    hr1(n,:) = ((zrr1 - rsn)/dr1(n)).';
    hr0(n,:) = ((zrr0 - rsn)/dr0(n)).';
end

Sigma_ini = blkdiag(diag(1^2*ones(1,3)), diag(1^2*ones(1,3)), ...
                    diag(0.3^2*ones(1,N_ml)), diag(0.3^2*ones(1,N_ml)));

Aini = zeros(16, 6 + 2*N_ml);
Aini(1:3, 1:3) = eye(3);
Aini(4:6, 1:3) = (1/T_ini)*eye(3);
Aini(4:6, 4:6) = (-1/T_ini)*eye(3);

for n = 1:N_ml
    bias_idx  = 6 + 2*n - 1;
    drift_idx = 6 + 2*n;

    Aini(bias_idx, 1:3) = -hr1(n,:);
    Aini(bias_idx, 6 + n) = 1;

    Aini(drift_idx, 1:3) = (-hr1(n,:))/T_ini;
    Aini(drift_idx, 4:6) = ( hr0(n,:))/T_ini;
    Aini(drift_idx, 6 + n)         =  1/T_ini;
    Aini(drift_idx, 6 + N_ml + n)  = -1/T_ini;
end

P_ML_ini = Aini * Sigma_ini * Aini.';
P_ML_ini = (P_ML_ini + P_ML_ini.')/2;
[eV,eD] = eig(P_ML_ini);
eD = max(eD, 1e-6*eye(size(eD)));
P = eV*eD*eV.';

%% ====================== 9) INITIAL STATE X  ======================
X = zeros(16,1);
X(1:3) = truePos(1,:)';
X(4:6) = trueVel(1,:)';

for i = 1:numTx
    idx_b = 6 + (i-1)*2 + 1;
    idx_d = 6 + (i-1)*2 + 2;
    X(idx_b) = bias_init(i);
    X(idx_d) = drift_init(i);
end

%% ====================== 10) EKF LOOP  ======================
fprintf('Running EKF from TXT epoch 1...\n');

X_hist = zeros(numEpochs,16);
sigma_pos = zeros(numEpochs,3);

nu_log = zeros(6, numEpochs);
S_log  = zeros(6, 6, numEpochs);

X_hist(1,:) = X';
sigma_pos(1,:) = sqrt(diag(P(1:3,1:3)))';

for k = 2:numEpochs
    % predict
    X_prior = F * X;
    P_prior = F * P * F' + Q;

    % measurement: 5 code + 1 height constraint (use aligned GT height exactly like your old script)
    H = zeros(6, 16);
    h = zeros(6, 1);
    z = zeros(6, 1);

    pos_pred = X_prior(1:3);

    for i = 1:numTx
        delta = pos_pred - basePos(i,:)';
        dist  = norm(delta);
        u = delta / dist;

        idx_b = 6 + (i-1)*2 + 1;

        H(i,1:3) = u';
        H(i,idx_b) = 1;

        h(i) = dist + X_prior(idx_b);
        z(i) = measCode(k,i);
    end

    % height constraint
    H(6,3) = 1;
    h(6)   = pos_pred(3);
    z(6)   = truePos(k,3);

    % update
    S = H*P_prior*H' + R;
    K = P_prior*H'/S;

    nu = z - h;

    X = X_prior + K*nu;
    P = (eye(16) - K*H)*P_prior;

    X_hist(k,:) = X';
    sigma_pos(k,:) = sqrt(diag(P(1:3,1:3)))';

    nu_log(:,k)  = nu;
    S_log(:,:,k) = S;
end

fprintf('EKF finished.\n');

%% ====================== 11) ERROR / RMSE ======================
estPos = X_hist(:,1:3);
errPos = estPos - truePos;

err_E = errPos(:,1);
err_N = errPos(:,2);
err_U = errPos(:,3);
err_3d = sqrt(err_E.^2 + err_N.^2 + err_U.^2);

rmse_E = sqrt(mean(err_E.^2));
rmse_N = sqrt(mean(err_N.^2));
rmse_U = sqrt(mean(err_U.^2));
rmse_3D = sqrt(mean(err_3d.^2));

fprintf('RMSE (m): E=%.3f, N=%.3f, U=%.3f, 3D=%.3f\n', rmse_E, rmse_N, rmse_U, rmse_3D);

%% ====================== 12) PLOTS======================
tmin = (obsTime_abs - obsTime_abs(1))/60;

% --- use sigma from EKF ---
sigma_E = sigma_pos(:,1);
sigma_N = sigma_pos(:,2);
sigma_U = sigma_pos(:,3);

b3_E = 3 * sigma_E;
b3_N = 3 * sigma_N;
b3_U = 3 * sigma_U;

% optional: ignore first N points when auto-scaling y-limits
ignore_steps = 200;
len_plot = numEpochs;

if len_plot > ignore_steps
    scale_idx = ignore_steps:len_plot;
else
    scale_idx = 1:len_plot;
end

lim_E = max(abs([b3_E(scale_idx); err_E(scale_idx)]));
lim_N = max(abs([b3_N(scale_idx); err_N(scale_idx)]));
lim_U = max(abs([b3_U(scale_idx); err_U(scale_idx)]));

lim_E = max(lim_E, 0.1);
lim_N = max(lim_N, 0.1);
lim_U = max(lim_U, 0.1);

% ---- Plot: 3 panels with 3-sigma envelope ----
fill_color  = [1, 0.85, 0.85];     % light red background
bound_color = [0.85, 0.325, 0.098];% red dashed bound
line_color  = [0, 0.4470, 0.7410]; % MATLAB default blue

figure('Units','pixels','Position',[150, 150, 800, 700], 'Color','w');

% Common fill time polygon
t_fill = [tmin(:); flipud(tmin(:))];

% -------- East --------
subplot(3,1,1); hold on; box on; grid on;
y_fill_E = [b3_E(:); flipud(-b3_E(:))];
fill(t_fill, y_fill_E, fill_color, 'EdgeColor','none', 'FaceAlpha',0.5, ...
    'DisplayName','3\sigma Envelope');
plot(tmin,  b3_E, '--', 'Color', bound_color, 'LineWidth', 1, 'HandleVisibility','off');
plot(tmin, -b3_E, '--', 'Color', bound_color, 'LineWidth', 1, 'HandleVisibility','off');
plot(tmin, err_E, '-', 'Color', line_color, 'LineWidth', 1.2, ...
    'DisplayName','Estimation Error');
ylabel('East Error (m)', 'FontSize', 11);
legend('Location','northeast', 'NumColumns', 2);
set(gca,'FontSize',10);
xlim([tmin(1), tmin(end)]);
ylim([-lim_E*1.5, lim_E*1.5]);

% -------- North --------
subplot(3,1,2); hold on; box on; grid on;
y_fill_N = [b3_N(:); flipud(-b3_N(:))];
fill(t_fill, y_fill_N, fill_color, 'EdgeColor','none', 'FaceAlpha',0.5, ...
    'DisplayName','3\sigma Envelope');
plot(tmin,  b3_N, '--', 'Color', bound_color, 'LineWidth', 1, 'HandleVisibility','off');
plot(tmin, -b3_N, '--', 'Color', bound_color, 'LineWidth', 1, 'HandleVisibility','off');
plot(tmin, err_N, '-', 'Color', line_color, 'LineWidth', 1.2, ...
    'DisplayName','Estimation Error');
ylabel('North Error (m)', 'FontSize', 11);
legend('Location','northeast', 'NumColumns', 2);
set(gca,'FontSize',10);
xlim([tmin(1), tmin(end)]);
ylim([-lim_N*1.5, lim_N*1.5]);

% -------- Up --------
subplot(3,1,3); hold on; box on; grid on;
y_fill_U = [b3_U(:); flipud(-b3_U(:))];
fill(t_fill, y_fill_U, fill_color, 'EdgeColor','none', 'FaceAlpha',0.5, ...
    'DisplayName','3\sigma Envelope');
plot(tmin,  b3_U, '--', 'Color', bound_color, 'LineWidth', 1, 'HandleVisibility','off');
plot(tmin, -b3_U, '--', 'Color', bound_color, 'LineWidth', 1, 'HandleVisibility','off');
plot(tmin, err_U, '-', 'Color', line_color, 'LineWidth', 1.2, ...
    'DisplayName','Estimation Error');
ylabel('Up Error (m)', 'FontSize', 11);
xlabel('Time (min)', 'FontSize', 11);
legend('Location','northeast', 'NumColumns', 2);
set(gca,'FontSize',10);
xlim([tmin(1), tmin(end)]);
ylim([-lim_U*1.5, lim_U*1.5]);