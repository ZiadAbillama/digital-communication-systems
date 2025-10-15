% EECE 442 - Project Part 1 Main Script
clc; clear; close all;

% === SECTION 1: Sampling and Reconstruction ===
run('experiments/sampling_testing.m');
run('experiments/reconstruction_testing.m');

% === SECTION 2: Fourier Series Sampling ===
run('modules/sampling/ffs_test.m');

% === SECTION 3: Quantization ===
% run('experiments/quantization_testing.m');

% === SECTION 4: Source Coding ===
% run('experiments/source_coding_testing.m');
