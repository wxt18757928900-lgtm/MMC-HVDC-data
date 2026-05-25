clear; close all; clc;

scriptDir = fileparts(mfilename('fullpath'));

run(fullfile(scriptDir, 'fit_arm_reactor.m'));
run(fullfile(scriptDir, 'fit_transformer_dm.m'));
run(fullfile(scriptDir, 'fit_transformer_cm.m'));

fprintf('\nAll fitting scripts finished. Results are saved in the results folder.\n');

