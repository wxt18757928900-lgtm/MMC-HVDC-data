clear; close all; clc;

repoRoot = fileparts(fileparts(mfilename('fullpath')));
dataFile = fullfile(repoRoot, 'data', 'arm_reactor_impedance.csv');
resultDir = fullfile(repoRoot, 'results');
if ~exist(resultDir, 'dir')
    mkdir(resultDir);
end

if exist('rationalfit', 'file') ~= 2
    error('fit_arm_reactor requires MATLAB RF Toolbox: rationalfit was not found.');
end

data = readtable(dataFile);
freq = data.frequency_hz;
phase_deg = data.phase_deg;
mag = data.magnitude_ohm;
w = 2*pi*freq;

Z_meas = mag .* exp(1i*deg2rad(phase_deg));
Y_meas = 1 ./ Z_meas;
weights = 1 ./ abs(Y_meas);

target_poles = 20;
best_fit = [];
final_error_db = NaN;

fprintf('\n[Arm reactor] Admittance-domain vector fitting\n');
for n = 4:2:target_poles
    warning('off', 'RF:rationalfit:ErrorToleranceNotMet');
    fit_temp = rationalfit(freq, Y_meas, 'NPoles', n, ...
        'Weight', weights, 'Tolerance', -40);
    warning('on', 'RF:rationalfit:ErrorToleranceNotMet');

    Y_fit = freqresp(fit_temp, freq);
    avg_error_db = 20*log10(mean(abs(Y_meas - Y_fit) ./ abs(Y_meas)));
    fprintf('  poles = %2d, average relative error = %.2f dB\n', n, avg_error_db);

    if n == target_poles
        best_fit = fit_temp;
        final_error_db = avg_error_db;
    end
end

poles = best_fit.A;
residues = best_fit.C;
D = best_fit.D;
E = best_fit.E;

raw_params = [];
count = 1;

if abs(real(E)) > 1e-18
    raw_params(count,:) = [1, 0, 0, real(E), 0];
    count = count + 1;
end

if abs(real(D)) > 1e-9
    raw_params(count,:) = [2, 1/real(D), 0, 0, 0];
    count = count + 1;
end

processed = [];
for k = 1:length(poles)
    if ismember(k, processed)
        continue;
    end

    p = poles(k);
    r = residues(k);

    if abs(imag(p)) < 1e-4 * abs(real(p))
        L = 1/real(r);
        R = -real(p)/real(r);
        if L > 0 && R > 0
            raw_params(count,:) = [3, R, L, 0, 0];
            count = count + 1;
        end
        processed(end+1) = k;
    elseif k < length(poles)
        alpha = -real(p);
        beta = abs(imag(p));
        r_re = real(r);
        L = 1/(2*r_re);
        R = 2*alpha*L;
        C = 1/((alpha^2 + beta^2)*L);
        f0 = 1/(2*pi*sqrt(L*C));
        if L > 0 && R > 0 && C > 0
            raw_params(count,:) = [4, R, L, C, f0];
            count = count + 1;
        end
        processed(end+1) = k;
        processed(end+1) = k + 1;
    end
end

final_circuit = [];
rl_idxs = find(raw_params(:,1) == 3);
[~, max_l_loc] = max(raw_params(rl_idxs, 3));
main_L_idx = rl_idxs(max_l_loc);
main_L_val = raw_params(main_L_idx, 3);

for i = 1:size(raw_params, 1)
    type = raw_params(i, 1);
    R = raw_params(i, 2);
    L = raw_params(i, 3);
    keep_flag = true;

    if type == 3 && abs(L - main_L_val) > 1e-5
        if R < 200 && L < 1e-6
            keep_flag = false;
        end
        if R > 2000 && R < 3000 && L < 1e-3
            keep_flag = false;
        end
    end

    if keep_flag
        final_circuit = [final_circuit; raw_params(i,:)];
    end
end

R_damp = 200000;
final_circuit = [final_circuit; 2, R_damp, 0, 0, 0];

branch_type = strings(size(final_circuit,1), 1);
for i = 1:size(final_circuit, 1)
    switch final_circuit(i,1)
        case 1
            branch_type(i) = "shunt_capacitance";
        case 2
            if abs(final_circuit(i,2) - R_damp) < 1e-9
                branch_type(i) = "global_damping_resistor";
            else
                branch_type(i) = "shunt_resistor";
            end
        case 3
            branch_type(i) = "series_RL_branch";
        case 4
            branch_type(i) = "series_RLC_branch";
    end
end

params = table(branch_type, final_circuit(:,2), final_circuit(:,3), ...
    final_circuit(:,4), final_circuit(:,5), ...
    'VariableNames', {'branch_type','R_ohm','L_H','C_F','f_res_Hz'});
writetable(params, fullfile(resultDir, 'arm_reactor_fitted_parameters.csv'));

Y_fit_final = zeros(size(freq));
for i = 1:size(final_circuit, 1)
    type = final_circuit(i, 1);
    R = final_circuit(i, 2);
    L = final_circuit(i, 3);
    C = final_circuit(i, 4);

    if type == 1
        Y_branch = 1i*w*C;
    elseif type == 2
        Y_branch = 1/R;
    elseif type == 3
        Y_branch = 1 ./ (R + 1i*w*L);
    elseif type == 4
        Y_branch = 1 ./ (R + 1i*w*L + 1./(1i*w*C));
    end
    Y_fit_final = Y_fit_final + Y_branch;
end
Z_fit = 1 ./ Y_fit_final;

fig = figure('Name', 'Arm Reactor Model Verification', 'Color', 'w');
subplot(2,1,1);
loglog(freq, mag, 'b-', freq, abs(Z_fit), 'r--', 'LineWidth', 1.5);
grid on; grid minor;
xlabel('Frequency (Hz)');
ylabel('|Z| (Ohm)');
legend('Measured','Fitted','Location','best');
title(sprintf('Arm reactor magnitude, average admittance error %.2f dB', final_error_db));

subplot(2,1,2);
semilogx(freq, phase_deg, 'b-', freq, rad2deg(angle(Z_fit)), 'r--', 'LineWidth', 1.5);
grid on; grid minor;
xlabel('Frequency (Hz)');
ylabel('Phase (deg)');
legend('Measured','Fitted','Location','best');
title('Arm reactor phase');
saveas(fig, fullfile(resultDir, 'arm_reactor_verification.png'));

fprintf('[Arm reactor] Results written to %s\n', resultDir);

