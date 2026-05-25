clear; close all; clc;

repoRoot = fileparts(fileparts(mfilename('fullpath')));
dataFile = fullfile(repoRoot, 'data', 'transformer_dm_impedance.csv');
resultDir = fullfile(repoRoot, 'results');
if ~exist(resultDir, 'dir')
    mkdir(resultDir);
end

if exist('lsqnonlin', 'file') ~= 2
    error('fit_transformer_dm requires MATLAB Optimization Toolbox: lsqnonlin was not found.');
end

data = readtable(dataFile);
freq_meas = data.frequency_hz;
phase_meas_deg = data.phase_deg;
Z_mag_meas = data.magnitude_ohm;

f_min = 1e3;
f_max = 30e6;
mask = (freq_meas >= f_min) & (freq_meas <= f_max);
f = freq_meas(mask);
w = 2*pi*f;
mag = Z_mag_meas(mask);
phs_deg = phase_meas_deg(mask);
Z_meas = mag .* exp(1i*deg2rad(phs_deg));

Z_log = log10(mag);
Z_smooth = smoothdata(Z_log, 'gaussian', 15);

f_main_mask = (f >= 5e3) & (f <= 40e3);
[~, max_idx_main_local] = max(mag(f_main_mask));
f_main_candidates = f(f_main_mask);
f_peak1 = f_main_candidates(max_idx_main_local);
peak1_global_idx = find(f == f_peak1, 1);

low_f_idx = find(f >= 1.5e3, 1);
L_eq_total = mag(low_f_idx) / w(low_f_idx);

[~, locs, ~, prom] = findpeaks(Z_smooth, ...
    'MinPeakProminence', 0.15, 'MinPeakDistance', 20);
valid_mask = f(locs) > (f_peak1 * 1.5);
locs = locs(valid_mask);
prom = prom(valid_mask);

target_other_branches = 11;
if length(locs) > target_other_branches
    [~, sort_idx] = sort(prom, 'descend');
    locs = sort(locs(sort_idx(1:target_other_branches)));
end

all_locs = [peak1_global_idx; locs(:)];
num_branches = numel(all_locs);

p0 = zeros(1, 3*num_branches);
lb = zeros(1, 3*num_branches);
ub = zeros(1, 3*num_branches);

for i = 1:num_branches
    idx_pt = all_locs(i);
    fr = f(idx_pt);
    Zr = mag(idx_pt);
    wr = 2*pi*fr;
    base = (i-1)*3;

    if i == 1
        R_guess = Zr;
        L_guess = L_eq_total * 0.95;
        C_guess = 1 / (wr^2 * L_guess);
        lb(base+1) = log10(1000);
        ub(base+1) = log10(500e3);
    else
        R_guess = Zr;
        C_guess = 300e-12;
        L_guess = 1 / (wr^2 * C_guess);
        lb(base+1) = log10(0.1);
        ub(base+1) = log10(100e3);
    end

    p0(base+1) = log10(R_guess);
    p0(base+2) = log10(L_guess);
    p0(base+3) = log10(C_guess);
    lb(base+2) = log10(1e-8);
    ub(base+2) = log10(100e-3);
    lb(base+3) = log10(5e-12);
    ub(base+3) = log10(50000e-12);
end

options = optimoptions('lsqnonlin', ...
    'Display', 'iter', ...
    'Algorithm', 'trust-region-reflective', ...
    'MaxFunctionEvaluations', 30000, ...
    'FunctionTolerance', 1e-9, ...
    'StepTolerance', 1e-9);

fprintf('\n[Transformer DM] Foster-I fitting with %d parallel RLC tanks\n', num_branches);
fun = @(p) dm_cost_function(p, w, f, Z_meas, peak1_global_idx);
p_opt = lsqnonlin(fun, p0, lb, ub, options);

branch = (1:num_branches).';
R_ohm = zeros(num_branches,1);
L_H = zeros(num_branches,1);
C_F = zeros(num_branches,1);
f_res_Hz = zeros(num_branches,1);

Z_fit = zeros(size(w));
for k = 1:num_branches
    base = (k-1)*3;
    R_ohm(k) = 10^p_opt(base+1);
    L_H(k) = 10^p_opt(base+2);
    C_F(k) = 10^p_opt(base+3);
    f_res_Hz(k) = 1/(2*pi*sqrt(L_H(k)*C_F(k)));

    Y_tank = 1/R_ohm(k) + 1./(1i*w*L_H(k)) + 1i*w*C_F(k);
    Z_fit = Z_fit + 1./Y_tank;
end

params = table(branch, R_ohm, L_H, C_F, f_res_Hz);
writetable(params, fullfile(resultDir, 'transformer_dm_fitted_parameters.csv'));

fig = figure('Name', 'Transformer DM Model Verification', 'Color', 'w');
subplot(2,1,1);
loglog(f, mag, 'b-', f, abs(Z_fit), 'r--', 'LineWidth', 1.5);
grid on; grid minor;
xlabel('Frequency (Hz)');
ylabel('|Z_{dm}| (Ohm)');
legend('Measured','Fitted','Location','best');
title('Transformer DM magnitude');

subplot(2,1,2);
semilogx(f, phs_deg, 'b-', f, rad2deg(angle(Z_fit)), 'r--', 'LineWidth', 1.5);
grid on; grid minor;
xlabel('Frequency (Hz)');
ylabel('Phase (deg)');
legend('Measured','Fitted','Location','best');
title('Transformer DM phase');
saveas(fig, fullfile(resultDir, 'transformer_dm_verification.png'));

fprintf('[Transformer DM] Results written to %s\n', resultDir);

function residuals = dm_cost_function(p_log, w, f, Z_meas, peak_idx)
    p = 10.^p_log;
    num_branches = length(p)/3;
    Z_total = zeros(size(w));

    for k = 1:num_branches
        base = (k-1)*3;
        R = p(base+1);
        L = p(base+2);
        C = p(base+3);
        Y_tank = 1/R + 1./(1i*w*L) + 1i*w*C;
        Z_total = Z_total + 1./Y_tank;
    end

    diff_mag = log10(abs(Z_total)) - log10(abs(Z_meas));
    diff_angle = angle(Z_total) - angle(Z_meas);

    W = ones(size(f));
    W((f >= 5e3) & (f <= 500e3)) = 2.0;
    peak_range = max(1, peak_idx-15) : min(length(f), peak_idx+15);
    W(peak_range) = 10.0;
    W(1:50) = 5.0;

    residuals = [diff_mag .* W * 1.5; diff_angle .* W * 0.5];
    residuals(isnan(residuals)) = 0;
end

