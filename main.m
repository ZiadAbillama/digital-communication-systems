%% ==========================================================
%  EECE 442 Project - Part 1 (End-to-End)
%  Sampling, Quantization (Uniform & Lloyd–Max),
%  Source Coding (Entropy & Huffman) + Saved Figures
%  ----------------------------------------------------------
%  This script matches your current structure:
%    sampling/      -> sample.m, reconstruct.m, ffs.m
%    quantization/  -> quan.m, lloyd_max.m
%    source_coding/ -> entropy.m, huffman_encode.m, huffman_decode.m
%    tests/         -> your test scripts (not required here)
%
%  Figures are saved to results/figures/.
%  Menu lets you run each module separately.
%  ----------------------------------------------------------
clc; clear; close all;

% --- Path setup
addpath('sampling', 'quantization', 'source_coding', 'tests');

% --- Results dir
resultsDir = fullfile(pwd, 'results', 'figures');
if ~exist(resultsDir, 'dir'), mkdir(resultsDir); end

disp('=== EECE 442 - Part 1: Main Menu ===');

while true
    fprintf('\nSelect a section to run:\n');
    fprintf('  1) Sampling & Reconstruction (incl. FS error plots)\n');
    fprintf('  2) Quantization (Uniform vs Lloyd–Max + MSE vs M)\n');
    fprintf('  3) Source Coding (Huffman integrated with quantizers)\n');
    fprintf('  4) Exit\n');
    choice = input('Enter choice (1–4): ');

    switch choice
        %% =========================================
        % 1) Sampling & Reconstruction
        % =========================================
        case 1
            disp('--- Sampling & Reconstruction ---');

            % Message: cosine at f0 Hz (as in requirements)
            f0 = 10;                           % Hz
            t  = linspace(0, 1, 20000);        % "continuous" time
            xt = cos(2*pi*f0*t);

            % Nyquist freq for cosine is fN = f0 (Nyquist rate = 2*f0)
            fN = f0;

            % Two sampling rates: 0.5 fN (aliasing) and 2 fN (Nyquist rate)
            fs_low  = 0.5 * fN;                % below Nyquist -> aliasing
            fs_high = 2.0 * fN;                % at Nyquist rate

            % Sample (your signature: [t_sample, x_sample] = sample(t, xt, fs))
            [tL, xL] = sample(t, xt, fs_low);
            [tH, xH] = sample(t, xt, fs_high);

            % Reconstruct (your signature: xrcon = reconstruct(t, x_sample, fs))
            xL_rec = reconstruct(t, xL, fs_low);
            xH_rec = reconstruct(t, xH, fs_high);

            % Plot comparisons
            figure('Name','Sampling_Reconstruction','NumberTitle','off');
            subplot(2,1,1);
            plot(t, xt, 'b'); hold on;
            stem(tL, xL, 'r', 'filled');
            plot(t, xL_rec, 'k--', 'LineWidth', 1.1);
            title(sprintf('Below Nyquist: fs = %.2f Hz (0.5 f_N)', fs_low));
            legend('Original x(t)','Samples','Reconstruction'); grid on;

            subplot(2,1,2);
            plot(t, xt, 'b'); hold on;
            stem(tH, xH, 'm', 'filled');
            plot(t, xH_rec, 'k--', 'LineWidth', 1.1);
            title(sprintf('At Nyquist rate: fs = %.2f Hz (2 f_N)', fs_high));
            legend('Original x(t)','Samples','Reconstruction'); grid on;

            saveas(gcf, fullfile(resultsDir, 'Sampling_Reconstruction.png'));

            % -------- Finite Fourier Series error plots (requirements d,e)
            % Use a square wave chunk to show Gibbs & error trends
            tFS = linspace(-1, 1, 4000);
            xFS = square(2*pi*2*tFS);  % 2 Hz
            T   = 2;                   % matches window (-1..1)
            n_values = [1,3,5,10,20,40,80];
            E_n = zeros(size(n_values));
            for i = 1:numel(n_values)
                [xhat, ~] = ffs(xFS, tFS, n_values(i), T);
                E_n(i) = trapz(tFS, abs(xFS - xhat).^2);
            end

            T_values = [0.5, 1, 2, 4, 8];  % vary T at fixed n
            n_fixed  = 20;
            E_T = zeros(size(T_values));
            for i = 1:numel(T_values)
                [xhat, ~] = ffs(xFS, tFS, n_fixed, T_values(i));
                E_T(i) = trapz(tFS, abs(xFS - xhat).^2);
            end

            figure('Name','FFS_Error','NumberTitle','off');
            subplot(1,2,1);
            semilogy(n_values, E_n, 'or-','LineWidth',1.2);
            xlabel('n'); ylabel('E(n,T)'); grid on;
            title(sprintf('Error vs n (T=%.1f)', T));

            subplot(1,2,2);
            semilogy(T_values, E_T, 'ob-','LineWidth',1.2);
            xlabel('T (sec)'); ylabel('E(n,T)'); grid on;
            title(sprintf('Error vs T (n=%d)', n_fixed));

            saveas(gcf, fullfile(resultsDir, 'FFS_Error_Analysis.png'));
            disp('✅ Sampling + FS error figures saved.');

        %% =========================================
        % 2) Quantization (Uniform vs Lloyd–Max)
        % =========================================
        case 2
            disp('--- Quantization (Uniform vs Lloyd–Max) ---');

            % Use a representative analog sequence (same as tests)
            t = linspace(0, 1, 2000);
            x = sin(2*pi*3*t) + 0.3*sin(2*pi*7*t);

            % Choose M for side-by-side comparison
            M = 8;

            % --- Uniform thresholds & levels
            xmin = min(x); xmax = max(x);
            thrU = linspace(xmin, xmax, M+1); thrU = thrU(2:end-1);
            lvlU = linspace(xmin+(xmax-xmin)/(2*M), xmax-(xmax-xmin)/(2*M), M);
            xqU  = quan(x, thrU, lvlU);

            % --- Lloyd–Max (your function returns [xq, thr, lvl])
            [xqL, thrL, lvlL] = lloyd_max(x, M);

            % MSEs
            mseU = mean((x - xqU).^2);
            mseL = mean((x - xqL).^2);

            fprintf('Uniform MSE (M=%d):   %.6e\n', M, mseU);
            fprintf('Lloyd–Max MSE (M=%d): %.6e\n', M, mseL);

            % Plot quantized vs original
            figure('Name','Quantization_Comparison','NumberTitle','off');
            subplot(2,1,1);
            plot(t, x, 'b'); hold on; stairs(t, xqU, 'r','LineWidth',1);
            title(sprintf('Uniform Quantization (M=%d), MSE=%.2e', M, mseU));
            legend('Original','Quantized'); grid on;

            subplot(2,1,2);
            plot(t, x, 'b'); hold on; stairs(t, xqL, 'm','LineWidth',1);
            title(sprintf('Lloyd–Max Quantization (M=%d), MSE=%.2e', M, mseL));
            legend('Original','Quantized'); grid on;

            saveas(gcf, fullfile(resultsDir, sprintf('Quantization_Comparison_M%d.png', M)));

            % --- MSE vs M (both Uniform and Lloyd–Max)
            Mvals = [2 4 8 16 32 64 128];
            mseU_vec = zeros(size(Mvals));
            mseL_vec = zeros(size(Mvals));
            for i = 1:numel(Mvals)
                Mi = Mvals(i);
                % uniform
                thr = linspace(xmin, xmax, Mi+1); thr = thr(2:end-1);
                lvl = linspace(xmin+(xmax-xmin)/(2*Mi), xmax-(xmax-xmin)/(2*Mi), Mi);
                xq  = quan(x, thr, lvl);
                mseU_vec(i) = mean((x - xq).^2);
                % lloyd
                [xqOpt, ~, ~] = lloyd_max(x, Mi);
                mseL_vec(i) = mean((x - xqOpt).^2);
            end

            figure('Name','MSE_vs_M','NumberTitle','off');
            semilogy(Mvals, mseU_vec, 'r-o','LineWidth',1.1); hold on;
            semilogy(Mvals, mseL_vec, 'm-s','LineWidth',1.1);
            grid on; xlabel('Number of levels (M)'); ylabel('MSE');
            title('Quantization Error vs Number of Levels');
            legend('Uniform','Lloyd–Max','Location','southwest');

            saveas(gcf, fullfile(resultsDir, 'Quantization_MSE_vs_M.png'));
            disp('✅ Quantization figures saved.');

        %% =========================================
        % 3) Source Coding (Integrated with Quantizers)
        % =========================================
        case 3
            disp('--- Source Coding (Huffman) w/ Quantizer Outputs ---');

            % Use same analog sequence as quantization section
            t = linspace(0, 1, 2000);
            x = sin(2*pi*3*t) + 0.3*sin(2*pi*7*t);

            % Compare a few M values
            Mvals = [4 8 16];

            % Containers for a small summary table in console
            summary = [];  % [M, H_uni, Lbar_uni, Eff_uni, H_ll, Lbar_ll, Eff_ll, MSE_uni, MSE_ll]

            for Mi = Mvals
                % ----- Uniform quantization
                xmin = min(x); xmax = max(x);
                thr = linspace(xmin, xmax, Mi+1); thr = thr(2:end-1);
                lvl = linspace(xmin+(xmax-xmin)/(2*Mi), xmax-(xmax-xmin)/(2*Mi), Mi);
                xqU = quan(x, thr, lvl);
                mseU = mean((x - xqU).^2);

                % ----- Lloyd–Max
                [xqL, ~, ~] = lloyd_max(x, Mi);
                mseL = mean((x - xqL).^2);

                % ----- Huffman (your huffman_encode returns [dict, avglen])
                % Uniform
                [dictU, LbarU] = huffman_encode(xqU);
                HU   = entropy(xqU);
                encU = huffmanenco(num2cell(xqU), dictU);
                decU = huffman_decode(encU, dictU);
                losslessU = isequal(xqU(:), decU(:));
                % Lloyd–Max
                [dictL, LbarL] = huffman_encode(xqL);
                HL   = entropy(xqL);
                encL = huffmanenco(num2cell(xqL), dictL);
                decL = huffman_decode(encL, dictL);
                losslessL = isequal(xqL(:), decL(:));

                EffU = 100 * (HU / LbarU);
                EffL = 100 * (HL / LbarL);

                fprintf('\nM = %d\n', Mi);
                fprintf('  Uniform:   H=%.4f, Lbar=%.4f, Eff=%.2f%%, MSE=%.3e, Lossless=%d\n', ...
                    HU, LbarU, EffU, mseU, losslessU);
                fprintf('  Lloyd–Max: H=%.4f, Lbar=%.4f, Eff=%.2f%%, MSE=%.3e, Lossless=%d\n', ...
                    HL, LbarL, EffL, mseL, losslessL);

                summary = [summary; Mi, HU, LbarU, EffU, HL, LbarL, EffL, mseU, mseL]; %#ok<AGROW>
            end

            % Save a small summary figure (table-like)
            fig = figure('Name','RateDistortionSummary','NumberTitle','off','Visible','on');
            uitable('Data', summary, ...
                'ColumnName', {'M','H_U','Lbar_U','Eff_U(%)','H_LM','Lbar_LM','Eff_LM(%)','MSE_U','MSE_LM'}, ...
                'Units','normalized','Position',[0 0 1 1]);
            title('Rate/Distortion/Efficiency Summary');
            fpath = fullfile(resultsDir, 'Huffman_Rate_Distortion_Summary.png');
            saveas(fig, fpath);

            % Also save a histogram + level-placement figure for one M
            Mi = 8;  % pick an illustrative M
            thr = linspace(min(x), max(x), Mi+1); thr = thr(2:end-1);
            lvl = linspace(min(x)+(max(x)-min(x))/(2*Mi), max(x)-(max(x)-min(x))/(2*Mi), Mi);
            xqU = quan(x, thr, lvl);
            [xqL, ~, ~] = lloyd_max(x, Mi);

            figure('Name','LevelPlacement','NumberTitle','off');
            histogram(x, 100, 'Normalization','pdf','FaceAlpha',0.35); hold on;
            y0 = ylim; y0 = y0(2)*0.05;
            stem(lvl, y0*ones(size(lvl)), 'r','LineWidth',1.5);
            % For Lloyd levels, estimate from unique recon levels (mode not required here)
            lvlLM = unique(xqL); lvlLM = sort(lvlLM(:).');
            stem(lvlLM, (y0*1.2)*ones(size(lvlLM)), 'm','LineWidth',1.5);
            legend('Signal PDF','Uniform Levels','Lloyd–Max Levels','Location','best');
            title(sprintf('Level Placement (M=%d)', Mi)); grid on;

            saveas(gcf, fullfile(resultsDir, sprintf('Level_Placement_M%d.png', Mi)));
            disp('✅ Source coding metrics + figures saved.');

        %% =========================================
        % 4) Exit
        % =========================================
        case 4
            disp('Exiting.');
            break;

        otherwise
            disp('Invalid selection. Please enter 1–4.');
    end
end

disp('=== Done. All figures saved under results/figures/ ===');
