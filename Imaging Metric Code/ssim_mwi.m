% ssim_mwi.m
% Standalone SSIM for microwave imaging from .s6p files (RF Toolbox + IPT).
% Supports two antenna arrangements and several comparison views.
%
% ARR#1 (2x3):  Array_E_CutModel_vivaldi.s6p, Array_F_CutModel_vivaldi.s6p, antenna_locations.csv
% ARR#2 (1x6):  Array2_E_CutModel_vivaldi.s6p, Array2_F_CutModel_vivaldi.s6p, antenna2_locations.csv
%
% Run:
%   cd ~/Desktop/MERIT
%   ssim_mwi
%
% Outputs (plots + CSVs):
%   - SSIM vs z (Array1 vs Array2 overlay) + difference curve
%   - Volume-mean SSIM bar chart (Array1 vs Array2)
%   - Cross-array SSIM vs z: Full↔Full, Empty↔Empty, Diff↔Diff
%   - Tables in ssim_results_pair.csv and ssim_results_crossarray.csv

clear; clc;

%% ------------------ CONFIG ------------------
% Imaging domain (ellipse "cut" model)
a = 0.1064;        % X semi-axis (m)
b = 0.1705;        % Y semi-axis (m)
height = 0.15;     % Z extent (m)
resolution = 2e-3; % voxel spacing (m)
cutY = 0.140;      % crop in Y (m)

% Medium
epsilon_r = 19.72; % real relative permittivity
c0 = 299792458;    % m/s
v  = c0 / sqrt(epsilon_r);

% Optional band-limit: set [] to use all
band_Hz = [];      % e.g., [1.5e9 3.0e9]

% S-parameter magnitude floor (in dB). Set [] to disable.
threshold_dB = -50;

% Z-slices to evaluate SSIM on
z_vals = 0:0.01:height;

% Outputs
out_dir = pwd;
pair_csv  = fullfile(out_dir, 'ssim_results_pair.csv');
cross_csv = fullfile(out_dir, 'ssim_results_crossarray.csv');

%% ------------------ Check toolboxes ------------------
if isempty(which('sparameters')), error('RF Toolbox required.'); end
if isempty(which('ssim')), error('Image Processing Toolbox required.'); end

%% ------------------ Locate files ------------------
Folder = pwd;

% Array #1 (2x3)
fE1 = fullfile(Folder,'Array_E_CutModel_vivaldi.s6p');
fF1 = fullfile(Folder,'Array_F_CutModel_vivaldi.s6p');
ants1 = fullfile(Folder,'antenna_locations.csv');

% Array #2 (1x6)
fE2 = fullfile(Folder,'Array2_E_CutModel_vivaldi.s6p');
fF2 = fullfile(Folder,'Array2_F_CutModel_vivaldi.s6p');
ants2 = fullfile(Folder,'antenna2_locations.csv');

haveArr1 = isfile(fE1) && isfile(fF1) && isfile(ants1);
haveArr2 = isfile(fE2) && isfile(fF2) && isfile(ants2);
if ~haveArr1 && ~haveArr2
    error('No complete array dataset found. Need either ARR#1 or ARR#2 files present.');
end

%% ------------------ Build imaging grid ------------------------
[xv, yv, zv] = deal(0:resolution:a, -cutY:resolution:cutY, 0:resolution:height);
[nx, ny, nz] = deal(numel(xv), numel(yv), numel(zv));
[X, Y, Z]    = meshgrid(xv, yv, zv);
mask = (X./a).^2 + (Y./b).^2 <= 1;
pts  = [X(mask), Y(mask), Z(mask)];
[ix_all, iy_all, iz_all] = ind2sub([nx,ny,nz], find(mask));
Npts = size(pts,1);

%% ------------------ Helpers ------------------
% Read .s6p to [Nf x N^2] and channel list
function [freqHz, scan, channels, Nports] = read_snp_as_matrix(fname, threshold_dB, band_Hz)
    S = sparameters(fname);
    freqHz = S.Frequencies(:);
    S3 = S.Parameters;         % [N x N x Nf]
    N  = size(S3,1);
    if ~isempty(threshold_dB)
        SdB = 20*log10(abs(S3));
        S3(SdB < threshold_dB) = 0;
    end
    % Band-limit
    keep = true(size(freqHz));
    if ~isempty(band_Hz)
        keep = freqHz >= band_Hz(1) & freqHz <= band_Hz(2);
    end
    freqHz = freqHz(keep); S3 = S3(:,:,keep);
    % Flatten
    Nf = numel(freqHz); Nch = N*N;
    scan = zeros(Nf, Nch);
    channels = zeros(Nch,2);
    idx = 1;
    for r = 1:N
        for c = 1:N
            scan(:,idx)   = squeeze(S3(r,c,:));
            channels(idx,:) = [r c];
            idx = idx+1;
        end
    end
    Nports = N;
end

% DAS reconstruction
function img_vals = recon_das(scan, freqHz, channels, ant_pos, pts, v)
    k = 2*pi*freqHz(:)/v;          % [Nf x 1]
    Nch = size(channels,1);
    img_vals = zeros(size(pts,1),1);
    for ch = 1:Nch
        tx = channels(ch,1); rx = channels(ch,2);
        rt = ant_pos(tx,:);  rr = ant_pos(rx,:);
        dt = sqrt(sum((pts - rt).^2, 2));
        dr = sqrt(sum((pts - rr).^2, 2));
        L  = dt + dr;                             % [Npts x 1]
        phase = exp(-1j * (k * L.'));             % [Nf x Npts]
        img_vals = img_vals + (phase.' * scan(:,ch));  % [Npts x 1]
    end
    img_vals = abs(img_vals);
end

% Scatter values back to grid (NaN outside ellipse)
function V = vals_to_grid(vals, nx, ny, nz, ix, iy, iz)
    V = nan(nx,ny,nz);
    for n = 1:numel(ix)
        V(ix(n), iy(n), iz(n)) = vals(n);
    end
end

% Per-slice SSIM with zero-mean, unit-std normalization (paper-style)
function [ssim_per_z, vol_mean, tableZ] = ssim_volume(Agrid, Bgrid, xv, yv, zv, z_vals, global_max)
    if nargin<7 || isempty(global_max)
        global_max = max([max(Agrid,[],'all','omitnan'), max(Bgrid,[],'all','omitnan')]);
        if isempty(global_max) || global_max==0, global_max = 1; end
    end
    [~, iz_req] = arrayfun(@(z0) min(abs(zv - z0)), z_vals);
    ssim_per_z = zeros(numel(z_vals),1);
    for kk = 1:numel(z_vals)
        iz = iz_req(kk);
        A = squeeze(Agrid(:,:,iz)).';
        B = squeeze(Bgrid(:,:,iz)).';
        A = A/global_max; B = B/global_max;
        A(~isfinite(A)) = 0; B(~isfinite(B)) = 0;
        % Zero-mean, unit-std per slice (matches the passage)
        Am = mean(A,'all'); Bm = mean(B,'all');
        As = std(A,0,'all'); if As==0, As=1; end
        Bs = std(B,0,'all'); if Bs==0, Bs=1; end
        A = (A - Am)/As; B = (B - Bm)/Bs;
        ssim_per_z(kk) = ssim(A,B);
    end
    vol_mean = mean(ssim_per_z);
    tableZ = table(z_vals(:), ssim_per_z, 'VariableNames', {'z_m','SSIM'});
end

% Utility: per-slice cross-array SSIM curve + mean
function [ssim_z, ssim_mean] = crossarray_ssim(VA, VB, xv, yv, zv, z_vals)
    gmax = max([max(VA,[],'all','omitnan'), max(VB,[],'all','omitnan')]); if gmax==0, gmax=1; end
    [ssim_z, ssim_mean, ~] = ssim_volume(VA, VB, xv, yv, zv, z_vals, gmax);
end

%% ------------------ ARR#1: reconstruct & SSIM (Full vs Empty) ------------
arr1 = struct();
if haveArr1
    ant1 = readmatrix(ants1);  assert(size(ant1,2) >= 3, 'antenna_locations.csv must have >=3 columns.');
    [fE1_Hz, scanE1, ch1, Nports1] = read_snp_as_matrix(fE1, threshold_dB, band_Hz);
    [fF1_Hz, scanF1, ch1b, ~]       = read_snp_as_matrix(fF1, threshold_dB, band_Hz);
    assert(isequal(fE1_Hz,fF1_Hz) && isequal(ch1,ch1b), 'ARR1 Empty/Full grids or channels differ.');
    arr1.freq = fE1_Hz; arr1.channels = ch1; arr1.ant = ant1(1:Nports1,1:3);
    valsE1 = recon_das(scanE1, arr1.freq, arr1.channels, arr1.ant, pts, v);
    valsF1 = recon_das(scanF1, arr1.freq, arr1.channels, arr1.ant, pts, v);
    arr1.VE = vals_to_grid(valsE1, nx, ny, nz, ix_all, iy_all, iz_all);
    arr1.VF = vals_to_grid(valsF1, nx, ny, nz, ix_all, iy_all, iz_all);
    gmax1 = max([max(valsE1), max(valsF1)]); if gmax1==0, gmax1=1; end
    [arr1.ssim_z, arr1.ssim_mean, arr1.tblZ] = ssim_volume(arr1.VF, arr1.VE, xv, yv, zv, z_vals, gmax1);
end

%% ------------------ ARR#2: reconstruct & SSIM (Full vs Empty) ------------
arr2 = struct();
if haveArr2
    ant2 = readmatrix(ants2);  assert(size(ant2,2) >= 3, 'antenna2_locations.csv must have >=3 columns.');
    [fE2_Hz, scanE2, ch2, Nports2] = read_snp_as_matrix(fE2, threshold_dB, band_Hz);
    [fF2_Hz, scanF2, ch2b, ~]       = read_snp_as_matrix(fF2, threshold_dB, band_Hz);
    assert(isequal(fE2_Hz,fF2_Hz) && isequal(ch2,ch2b), 'ARR2 Empty/Full grids or channels differ.');
    arr2.freq = fE2_Hz; arr2.channels = ch2; arr2.ant = ant2(1:Nports2,1:3);
    valsE2 = recon_das(scanE2, arr2.freq, arr2.channels, arr2.ant, pts, v);
    valsF2 = recon_das(scanF2, arr2.freq, arr2.channels, arr2.ant, pts, v);
    arr2.VE = vals_to_grid(valsE2, nx, ny, nz, ix_all, iy_all, iz_all);
    arr2.VF = vals_to_grid(valsF2, nx, ny, nz, ix_all, iy_all, iz_all);
    gmax2 = max([max(valsE2), max(valsF2)]); if gmax2==0, gmax2=1; end
    [arr2.ssim_z, arr2.ssim_mean, arr2.tblZ] = ssim_volume(arr2.VF, arr2.VE, xv, yv, zv, z_vals, gmax2);
end

%% ------------------ FIGURE A: SSIM vs z overlay + difference -------------
if haveArr1 && haveArr2
    figure('Name','Per-slice SSIM (Full vs Empty): Array1 vs Array2');
    plot(z_vals, arr1.ssim_z, 'o-','LineWidth',1.5); hold on;
    plot(z_vals, arr2.ssim_z, 's-','LineWidth',1.5);
    plot(z_vals, arr1.ssim_z - arr2.ssim_z, 'k--','LineWidth',1.2);
    grid on; xlabel('z (m)'); ylabel('SSIM (Full vs Empty)');
    legend({'Array1 (2x3)','Array2 (1x6)','Array1 − Array2'}, 'Location','best');
    title('Per-slice SSIM comparison by array');
elseif haveArr1
    figure('Name','Per-slice SSIM (Full vs Empty): Array1 only');
    plot(z_vals, arr1.ssim_z, 'o-','LineWidth',1.5); grid on;
    xlabel('z (m)'); ylabel('SSIM (Full vs Empty)');
    title(sprintf('Array1 per-slice SSIM; volume-mean = %.3f', arr1.ssim_mean));
else
    figure('Name','Per-slice SSIM (Full vs Empty): Array2 only');
    plot(z_vals, arr2.ssim_z, 's-','LineWidth',1.5); grid on;
    xlabel('z (m)'); ylabel('SSIM (Full vs Empty)');
    title(sprintf('Array2 per-slice SSIM; volume-mean = %.3f', arr2.ssim_mean));
end

%% ------------------ FIGURE B: Volume-mean SSIM bar chart -----------------
labels = strings(0,1); means = [];
if haveArr1, labels(end+1,1)="Array1 (2x3)"; means(end+1,1)=arr1.ssim_mean; end
if haveArr2, labels(end+1,1)="Array2 (1x6)"; means(end+1,1)=arr2.ssim_mean; end
figure('Name','Volume-mean SSIM (Full vs Empty)');
bar(means);
set(gca,'XTick',1:numel(means),'XTickLabel',labels,'XTickLabelRotation',10);
ylabel('Volume-mean SSIM'); grid on;
title('Volume-mean SSIM (Full vs Empty) by array');

% Write pair-mode CSV (per-slice + volume-mean for each available array)
Tpair = table();

if haveArr1
    T1 = arr1.tblZ;                         % columns: z_m, SSIM
    n1 = size(T1,1);
    ArrayCol = repmat({'Array1 (2x3)'}, n1, 1);   % use cellstr for wide compatibility
    MeanCol  = repmat(arr1.ssim_mean, n1, 1);
    T1 = [ table(ArrayCol, 'VariableNames', {'Array'}) , ...
           T1 , ...
           table(MeanCol,  'VariableNames', {'VolumeMeanSSIM'}) ];
    Tpair = [Tpair; T1];
end

if haveArr2
    T2 = arr2.tblZ;                         % columns: z_m, SSIM
    n2 = size(T2,1);
    ArrayCol = repmat({'Array2 (1x6)'}, n2, 1);
    MeanCol  = repmat(arr2.ssim_mean, n2, 1);
    T2 = [ table(ArrayCol, 'VariableNames', {'Array'}) , ...
           T2 , ...
           table(MeanCol,  'VariableNames', {'VolumeMeanSSIM'}) ];
    Tpair = [Tpair; T2];
end

if ~isempty(Tpair)
    writetable(Tpair, pair_csv);
    fprintf('Wrote pair-mode SSIM table: %s\n', pair_csv);
end

%% ------------------ CROSS-ARRAY COMPARISONS ------------------------------
% Cross-array SSIM vs z for: Full↔Full, Empty↔Empty, Diff↔Diff
if haveArr1 && haveArr2
    % Full vs Full
    [ssimFF_z, ssimFF_mean] = crossarray_ssim(arr1.VF, arr2.VF, xv, yv, zv, z_vals);
    % Empty vs Empty
    [ssimEE_z, ssimEE_mean] = crossarray_ssim(arr1.VE, arr2.VE, xv, yv, zv, z_vals);
    % Differential volumes (Full-Empty)
    Vdiff1 = arr1.VF - arr1.VE;  % same grid, so direct subtraction
    Vdiff2 = arr2.VF - arr2.VE;
    [ssimDD_z, ssimDD_mean] = crossarray_ssim(Vdiff1, Vdiff2, xv, yv, zv, z_vals);

    % =======================
    % FIGURE C: SSIM curves
    % =======================
    figure('Name','Cross-array SSIM vs z');
    plot(z_vals, ssimFF_z, 'o-','LineWidth',1.5); hold on;
    plot(z_vals, ssimEE_z, 's-','LineWidth',1.5);
    plot(z_vals, ssimDD_z, 'd-','LineWidth',1.5);
    grid on; xlabel('z (m)'); ylabel('SSIM (Array1 ↔ Array2)');
    legend({'Full↔Full','Empty↔Empty','(Full−Empty)↔(Full−Empty)'}, 'Location','best');
    title(sprintf('Cross-array SSIM; means = [%.3f, %.3f, %.3f]', ssimFF_mean, ssimEE_mean, ssimDD_mean));

    % Helper: extract XY slice (transpose so X horizontal, Y vertical)
    getXY = @(V, iz) (squeeze(V(:,:,iz)).');

    % How many worst slices to show
    Nshow = min(6, numel(z_vals));

    % ============================================================
    % FIGURE D (revised): Worst cross-array FULL differences
    % GIF-like normalization: one global max across BOTH arrays
    % ============================================================
    [~, orderF] = sort(ssimFF_z, 'ascend');
    worstIdxF = orderF(1:Nshow);

    overall_max_full = max([max(arr1.VF,[],'all','omitnan'), max(arr2.VF,[],'all','omitnan')]);
    if ~isfinite(overall_max_full) || overall_max_full==0, overall_max_full = 1; end

    figure('Name','Worst z-slices where Full (Arr1 vs Arr2) differ most (global-normalized)');
    tF = tiledlayout(2,Nshow,'TileSpacing','compact','Padding','compact');
    for jj = 1:Nshow
        iz = find(abs(zv - z_vals(worstIdxF(jj))) < 1e-9, 1);
        A = getXY(arr1.VF, iz); B = getXY(arr2.VF, iz);
        A(~isfinite(A)) = 0; B(~isfinite(B)) = 0;

        nexttile; imagesc(xv, yv, A/overall_max_full); axis image; set(gca,'YDir','normal');
        title(sprintf('Arr1 Full @ z=%.2f m', z_vals(worstIdxF(jj))));
        xlabel('X (m)'); ylabel('Y (m)'); colorbar; caxis([0 1]);

        nexttile; imagesc(xv, yv, B/overall_max_full); axis image; set(gca,'YDir','normal');
        title(sprintf('Arr2 Full @ z=%.2f m', z_vals(worstIdxF(jj))));
        xlabel('X (m)'); ylabel('Y (m)'); colorbar; caxis([0 1]);
    end
    title(tF,'Worst cross-array Full differences');

    % =============================================================
    % FIGURE E: Worst cross-array EMPTY differences
    % GIF-like normalization: one global max across BOTH arrays
    % =============================================================
    [~, orderE] = sort(ssimEE_z, 'ascend');
    worstIdxE = orderE(1:Nshow);

    overall_max_empty = max([max(arr1.VE,[],'all','omitnan'), max(arr2.VE,[],'all','omitnan')]);
    if ~isfinite(overall_max_empty) || overall_max_empty==0, overall_max_empty = 1; end

    figure('Name','Worst z-slices where Empty (Arr1 vs Arr2) differ most (global-normalized)');
    tE = tiledlayout(2,Nshow,'TileSpacing','compact','Padding','compact');
    for jj = 1:Nshow
        iz = find(abs(zv - z_vals(worstIdxE(jj))) < 1e-9, 1);
        A = getXY(arr1.VE, iz); B = getXY(arr2.VE, iz);
        A(~isfinite(A)) = 0; B(~isfinite(B)) = 0;

        nexttile; imagesc(xv, yv, A/overall_max_empty); axis image; set(gca,'YDir','normal');
        title(sprintf('Arr1 Empty @ z=%.2f m', z_vals(worstIdxE(jj))));
        xlabel('X (m)'); ylabel('Y (m)'); colorbar; caxis([0 1]);

        nexttile; imagesc(xv, yv, B/overall_max_empty); axis image; set(gca,'YDir','normal');
        title(sprintf('Arr2 Empty @ z=%.2f m', z_vals(worstIdxE(jj))));
        xlabel('X (m)'); ylabel('Y (m)'); colorbar; caxis([0 1]);
    end
    title(tE,'Worst cross-array Empty differences');

    % =======================================================================
    % FIGURE F: Worst cross-array (Full−Empty) differences — SIGNED display
    % Symmetric color limits around zero using robust 99th percentile
    % =======================================================================
    [~, orderD] = sort(ssimDD_z, 'ascend');
    worstIdxD = orderD(1:Nshow);

    % Robust symmetric range across BOTH arrays and ALL worst slices
    diffs_abs = [];
    for jj = 1:Nshow
        iz = find(abs(zv - z_vals(worstIdxD(jj))) < 1e-9, 1);
        d1 = getXY(Vdiff1, iz); d2 = getXY(Vdiff2, iz);
        diffs_abs = [diffs_abs; abs(d1(:)); abs(d2(:))]; %#ok<AGROW>
    end
    M = prctile(diffs_abs(isfinite(diffs_abs)), 99);
    if ~isfinite(M) || M<=0, M = 1; end
    clim_sym = [-M, +M];

    figure('Name','Worst z-slices where (Full−Empty) Arr1 vs Arr2 differ most (signed)');
    tD = tiledlayout(2,Nshow,'TileSpacing','compact','Padding','compact');
    for jj = 1:Nshow
        iz = find(abs(zv - z_vals(worstIdxD(jj))) < 1e-9, 1);
        A = getXY(Vdiff1, iz); B = getXY(Vdiff2, iz);
        A(~isfinite(A)) = 0; B(~isfinite(B)) = 0;

        nexttile; imagesc(xv, yv, A); axis image; set(gca,'YDir','normal');
        title(sprintf('Arr1 (Full−Empty) @ z=%.2f m', z_vals(worstIdxD(jj))));
        xlabel('X (m)'); ylabel('Y (m)'); colorbar; caxis(clim_sym);

        nexttile; imagesc(xv, yv, B); axis image; set(gca,'YDir','normal');
        title(sprintf('Arr2 (Full−Empty) @ z=%.2f m', z_vals(worstIdxD(jj))));
        xlabel('X (m)'); ylabel('Y (m)'); colorbar; caxis(clim_sym);
    end
    title(tD,'Worst cross-array (Full−Empty) differences');

    % ---------------------
    % CROSS CSV table
    % ---------------------
    Tcross = table(z_vals(:), ssimFF_z(:), ssimEE_z(:), ssimDD_z(:), ...
                   'VariableNames', {'z_m','SSIM_Full_Full','SSIM_Empty_Empty','SSIM_Diff_Diff'});
    nrows = size(Tcross, 1);
    Tcross = [ Tcross, ...
               table( repmat(ssimFF_mean, nrows, 1), ...
                      repmat(ssimEE_mean, nrows, 1), ...
                      repmat(ssimDD_mean, nrows, 1), ...
                      'VariableNames', {'VolumeMean_Full_Full','VolumeMean_Empty_Empty','VolumeMean_Diff_Diff'} ) ];
    writetable(Tcross, cross_csv);
    fprintf('Wrote cross-array SSIM table: %s\n', cross_csv);
end

    % -----------------------
    % Add a self-check curve
    % -----------------------
    % Choose one stable self-comparison (identical image vs itself).
    % Using Arr1 Full if available, else Arr2 Full (identical slices -> SSIM = 1).
    if haveArr1
        self_label = 'Self-check (Arr1 Full vs itself)';
    else
        self_label = 'Self-check (Arr2 Full vs itself)';
    end
    ssimSelf_z = ones(size(z_vals));  % exact-1 baseline for identical images

    % =======================
    % FIGURE C: SSIM curves
    % =======================
    figure('Name','Cross-array SSIM vs z');
    plot(z_vals, ssimFF_z, 'o-','LineWidth',1.5); hold on;
    plot(z_vals, ssimEE_z, 's-','LineWidth',1.5);
    plot(z_vals, ssimDD_z, 'd-','LineWidth',1.5);
    % 4th line: self-check (SSIM=1)
    plot(z_vals, ssimSelf_z, ':','LineWidth',1.5);
    grid on; xlabel('z (m)'); ylabel('SSIM (Array1 ↔ Array2)');
    legend({'Full↔Full','Empty↔Empty','(Full−Empty)↔(Full−Empty)', self_label}, 'Location','best');
    title(sprintf('Cross-array SSIM; means = [%.3f, %.3f, %.3f]', ssimFF_mean, ssimEE_mean, ssimDD_mean));

    % ===========================================
    % MERGED OVERVIEW: curves + volume-mean bars
    % ===========================================
    % Prepare bar data (reuse the volume-mean SSIM from pair mode)
    labels = strings(0,1); means = [];
    if haveArr1, labels(end+1,1) = "Array1 (2x3)"; means(end+1,1) = arr1.ssim_mean; end
    if haveArr2, labels(end+1,1) = "Array2 (1x6)"; means(end+1,1) = arr2.ssim_mean; end

    figMerged = figure('Name','Merged SSIM overview: curves + volume means');
    tl = tiledlayout(figMerged,1,3,'TileSpacing','compact','Padding','compact');

    % Left: the same cross-array SSIM vs z curves (with self-check), spanning 2 tiles
    ax1 = nexttile(tl,[1 2]);
    plot(z_vals, ssimFF_z, 'o-','LineWidth',1.5); hold on;
    plot(z_vals, ssimEE_z, 's-','LineWidth',1.5);
    plot(z_vals, ssimDD_z, 'd-','LineWidth',1.5);
    plot(z_vals, ssimSelf_z, ':','LineWidth',1.5);
    grid on; xlabel('z (m)'); ylabel('SSIM');
    legend({'Full↔Full','Empty↔Empty','(Full−Empty)↔(Full−Empty)', self_label}, 'Location','southwest');
    title('Cross-array SSIM vs depth (with self-check baseline)');

    % Right: compact volume-mean bar plot (pair-mode means for each array)
    ax2 = nexttile(tl);
    if ~isempty(means)
        bar(means);
        set(gca,'XTick',1:numel(means),'XTickLabel',labels,'XTickLabelRotation',10);
        ylabel('Volume-mean SSIM'); grid on;
        title('Volume-mean SSIM (Full vs Empty)');
        ylim([0 1]);  % SSIM in [0,1]
    else
        text(0.5,0.5,'No volume-mean data','HorizontalAlignment','center');
        axis off;
    end