clear; close all; clc;

repoRoot = fileparts(fileparts(mfilename('fullpath')));
dataFile = fullfile(repoRoot, 'data', 'transformer_cm_impedance.csv');
resultDir = fullfile(repoRoot, 'results');
if ~exist(resultDir, 'dir')
    mkdir(resultDir);
end

if exist('lsqnonlin', 'file') ~= 2
    error('fit_transformer_cm requires MATLAB Optimization Toolbox: lsqnonlin was not found.');
end

data = readtable(dataFile);
freq_meas = data.frequency_hz;
phase_meas_deg = data.phase_deg;
Z_mag_meas = data.magnitude_ohm * 2;

f_min = 1e3;
f_max = 30e6;
mask = (freq_meas >= f_min) & (freq_meas <= f_max);
f = freq_meas(mask);
w = 2*pi*f;
mag = Z_mag_meas(mask);
phs_deg = phase_meas_deg(mask);
Z_meas = mag .* exp(1i*deg2rad(phs_deg));

f_linear_idx = f < 50e3;
if sum(f_linear_idx) > 5
    p_fit = polyfit(log(w(f_linear_idx)), log(mag(f_linear_idx)), 1);
    C0_est = exp(-p_fit(2));
else
    C0_est = 100e-12;
end

Z_log = log10(mag);
Z_smooth = smoothdata(Z_log, 'gaussian', 25);
[~, locs, ~, prom] = findpeaks(-Z_smooth, ...
    'MinPeakProminence', 0.25, 'MinPeakDistance', 40);
valid_mask = f(locs) > 20e3;
locs = locs(valid_mask);
prom = prom(valid_mask);

target_branches = 12;
if length(locs) > target_branches
    [~, sort_idx] = sort(prom, 'descend');
    locs = sort(locs(sort_idx(1:target_branches)));
end

num_branches = length(locs);
p0 = zeros(1, 1 + 3*num_branches);
p0(1) = log10(C0_est);

for i = 1:num_branches
    idx_pt = locs(i);
    fr = f(idx_pt);
    Zr = mag(idx_pt);
    wr = 2*pi*fr;

    R_guess = min(Zr, 1000);
    C_guess = 200e-12;
    L_guess = 1 / (wr^2 * C_guess);

    base = 1 + (i-1)*3;
    p0(base+1) = log10(R_guess);
    p0(base+2) = log10(L_guess);
    p0(base+3) = log10(C_guess);
end

lb = zeros(size(p0));
ub = zeros(size(p0));
lb(1) = log10(50e-12);
ub(1) = log10(5000e-12);

for i = 1:num_branches
    base = 1 + (i-1)*3;
    lb(base+1) = log10(0.01);
    ub(base+1) = log10(2000);
    lb(base+2) = log10(1e-8);
    ub(base+2) = log10(100e-3);
    lb(base+3) = log10(10e-12);
    ub(base+3) = log10(20000e-12);
end

options = optimoptions('lsqnonlin', ...
    'Display', 'iter', ...
    'Algorithm', 'trust-region-reflective', ...
    'MaxFunctionEvaluations', 20000, ...
    'FunctionTolerance', 1e-9, ...
    'StepTolerance', 1e-9);

fprintf('\n[Transformer CM] Foster-II fitting with C0 and %d series RLC branches\n', num_branches);
fun = @(p) cm_cost_function(p, w, Z_meas);
p_opt = lsqnonlin(fun, p0, lb, ub, options);

C0_F = 10^p_opt(1);
branch = (1:num_branches).';
R_ohm = zeros(num_branches,1);
L_H = zeros(num_branches,1);
C_F = zeros(num_branches,1);
f_res_Hz = zeros(num_branches,1);

Y_fit = 1i*w*C0_F;
for k = 1:num_branches
    base = 1 + (k-1)*3;
    R_ohm(k) = 10^p_opt(base+1);
    L_H(k) = 10^p_opt(base+2);
    C_F(k) = 10^p_opt(base+3);
    f_res_Hz(k) = 1/(2*pi*sqrt(L_H(k)*C_F(k)));

    Z_branch = R_ohm(k) + 1i*w*L_H(k) + 1./(1i*w*C_F(k));
    Y_fit = Y_fit + 1./Z_branch;
end
Z_fit = 1 ./ Y_fit;

params = table(branch, R_ohm, L_H, C_F, f_res_Hz);
writetable(params, fullfile(resultDir, 'transformer_cm_fitted_parameters.csv'));
writetable(table(C0_F), fullfile(resultDir, 'transformer_cm_C0.csv'));

fig = figure('Name', 'Transformer CM Model Verification', 'Color', 'w');
subplot(2,1,1);
loglog(f, mag, 'b-', f, abs(Z_fit), 'r--', 'LineWidth', 1.5);
grid on; grid minor;
xlabel('Frequency (Hz)');
ylabel('|Z_{cm}| (Ohm)');
legend('Measured','Fitted','Location','best');
title(sprintf('Transformer CM magnitude, C0 = %.2f pF', C0_F*1e12));

subplot(2,1,2);
semilogx(f, phs_deg, 'b-', f, rad2deg(angle(Z_fit)), 'r--', 'LineWidth', 1.5);
grid on; grid minor;
xlabel('Frequency (Hz)');
ylabel('Phase (deg)');
legend('Measured','Fitted','Location','best');
title('Transformer CM phase');
saveas(fig, fullfile(resultDir, 'transformer_cm_verification.png'));

fprintf('[Transformer CM] Results written to %s\n', resultDir);

function residuals = cm_cost_function(p_log, w, Z_meas)
    p = 10.^p_log;
    C0 = p(1);
    num_branches = (length(p)-1)/3;

    Y_total = 1i*w*C0;
    for k = 1:num_branches
        base = 1 + (k-1)*3;
        R = p(base+1);
        L = p(base+2);
        C = p(base+3);
        Z_branch = R + 1i*w*L + 1./(1i*w*C);
        Y_total = Y_total + 1./Z_branch;
    end

    Z_total = 1 ./ Y_total;
    diff_mag = log10(abs(Z_total)) - log10(abs(Z_meas));
    diff_angle = angle(Z_total) - angle(Z_meas);
    residuals = [diff_mag * 1.2; diff_angle * 0.8];
    residuals(isnan(residuals)) = 0;
end

