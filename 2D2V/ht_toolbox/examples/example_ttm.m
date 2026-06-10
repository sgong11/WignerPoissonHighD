% EXAMPLE_TTM.M
% Demonstrates how to perform n-mode multiplication (TTM)
% on a full 3D tensor using ht.ttm.

addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..')));

% === Create a moderate-size 3D tensor ===
X = randn(10, 20, 30);  % size: [I1, I2, I3]

% === Mode-2 multiplication: project 2nd dimension (20 → 8) ===
A2 = randn(8, 20);                % shape: [8 × 20]
Y2 = ht.ttm(X, A2, 2);            % result: [10 × 8 × 30]

disp('--- TTM along mode-2 ---');
disp(['Original size: ', mat2str(size(X))]);
disp(['Projected size (mode-2): ', mat2str(size(Y2))]);

% === Mode-3 multiplication: project 3rd dimension (30 → 5) ===
A3 = randn(5, 30);                % shape: [5 × 30]
Y3 = ht.ttm(X, A3, 3);            % result: [10 × 20 × 5]

disp('--- TTM along mode-3 ---');
disp(['Projected size (mode-3): ', mat2str(size(Y3))]);

% === Mode-1 multiplication: project 1st dimension (10 → 6) ===
A1 = randn(6, 10);                % shape: [6 × 10]
Y1 = ht.ttm(X, A1, 1);            % result: [6 × 20 × 30]

disp('--- TTM along mode-1 ---');
disp(['Projected size (mode-1): ', mat2str(size(Y1))]);
