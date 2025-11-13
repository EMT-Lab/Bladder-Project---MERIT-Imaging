% SSIM between FULL and EMPTY reconstructions, per z-slice, reading .s6p files.
% Standalone version: no MERIT.* calls needed.
% Run from:  cd ~/Desktop/MERIT ;  ssim_zslices_from_touchstone

clear; clc;

%% ------------------ Folder & file discovery ------------------
Folder   = pwd;
ants_file = fullfile(Folder,'antenna_locations.csv');

% Expected names (Empty = E, Full = F). If F missing, pick any other .s6p
s6p_empty = fullfile(Folder,'Array_E_CutModel_vivaldi.s6p');
s6p_full  = fullfile(Folder,'Array_F_CutModel_vivaldi.s6p');

if ~isfile(s6p_empty)
    error('Missing %s in %s (Empty case).', 'Array_E_CutModel_vivaldi.s6p', Folder);
end
if ~isfile(s6p_full)
    d = dir(fullfile(Folder,'*.s6p'));
    alt = setdiff(string(fullfile({d.folder},{d.name})), string(s6p_empty));
    if isempty(alt)
        error('Full .s6p not found. Expected Array_F_CutModel_vivaldi.s6p or another .s6p besides the E file.');
    end
    s6p_full = char(alt(1));
    fprintf('WARNING: using %s as FULL (Array_F_CutModel_vivaldi.s6p not found).\n', s6p_full);
end

%% ------------------ Imaging domain (ellipse "cut" model) ------------------
a = 0.1064;         % X semi-axis (m)
b = 0.1705;         % Y semi-axis (m)
height = 0.15;      % Z extent (m)
resolution = 2e-3;  % voxel spacing (m)
cutY = 0.140;       % crop in Y (m)

% Medium & slices
epsilon_r = 19.72;             % real relative permittivity
c0 = 299792458;                % m/s
v  = c0 / sqrt(epsilon_r);     % wave speed in medium
z_vals = 0:0.01:height;        % slices to evaluate

% Optional S-parameter magnitude floor (in dB)
threshold_dB = -50;            % set [] to disable

% Output GIFs
gif_ssim_maps    = fullfile(Folder,'ssim_maps.gif');
gif_full_slices  = fullfile(Folder,'full_slices.gif');
gif_empty_slices = fullfile(Folder,'empty_slices.gif');

%% ------------------ Load Touchstone with RF Toolbox ----------------------
try
    S_F = sparameters(s6p_full);
    S_E = sparameters(s6p_empty);
catch ME
    error(['Failed to read .s6p via RF Toolbox (sparameters). ' ...
           'Toolbox not found or file parse error.\n%s'], ME.message);
end

frequencies = S_F.Frequencies(:);                 % Hz
assert(isequal(size(S_E.Frequencies), size(frequencies)) && ...
       all(abs(S_E.Frequencies - frequencies) < 1e-9), ...
       'Full and Empty .s6p must share the same frequency grid.');

SF = S_F.Parameters;   % [6 x 6 x Nf] complex
SE = S_E.Parameters;

% Optional magnitude floor
if ~isempty(threshold_dB)
    SFdB = 20*log10(abs(SF));  SF(SFdB < threshold_dB) = 0;
    SEdB = 20*log10(abs(SE));  SE(SEdB < threshold_dB) = 0;
end

% Flatten to [Nf x 36] (row-major: S11,S12,...,S66) + channel list
num_ports = size(SF,1);    % 6
Nf        = numel(frequencies);
Nch       = num_ports^2;
scan_full  = zeros(Nf,Nch);
scan_empty = zeros(Nf,Nch);
channels   = zeros(Nch,2);
idx = 1;
for r = 1:num_ports
    for c = 1:num_ports
        scan_full(:,idx)  = squeeze(SF(r,c,:));
        scan_empty(:,idx) = squeeze(SE(r,c,:));
        channels(idx,:)   = [r c];   % [tx rx] port indices
        idx = idx + 1;
    end
end

%% ------------------ Antenna locations -------------------------
antenna_locations = readmatrix(ants_file);   % [x y z] in meters, one row per port (1..6)
assert(size(antenna_locations,2) >= 3, 'antenna_locations.csv must have >= 3 columns [x y z].');
assert(size(antenna_locations,1) >= num_ports, 'antenna_locations must have at least %d rows.', num_ports);

%% ------------------ Build imaging grid & indexers ------------------------
[xv, yv, zv] = deal(0:resolution:a, -cutY:resolution:cutY, 0:resolution:height);
[nx, ny, nz] = deal(numel(xv), numel(yv), numel(zv));
[X, Y, Z]    = meshgrid(xv, yv, zv);       % Note: meshgrid returns Y,X,Z arrangement for images; we'll index carefully

% Elliptical mask on each Z
mask = (X./a).^2 + (Y./b).^2 <= 1;

% Linearized list of voxel centers inside ellipse
pts = [X(mask), Y(mask), Z(mask)];          % [Npts x 3]

% Build a mapping from voxel coords to linear index in a 3D array
img_full_grid  = nan(nx, ny, nz);           % we'll fill [ix,iy,iz] later
img_empty_grid = nan(nx, ny, nz);

% Precompute index triples for masked points
[ix_all, iy_all, iz_all] = ind2sub([nx,ny,nz], find(mask));
idx_pts = find(mask); %#ok<NASGU> % (not needed but kept for clarity)

%% ------------------ Standalone DAS beamforming ---------------------------
% Frequency wavenumbers in medium
k = 2*pi*frequencies / v;      % [Nf x 1], rad/m

% Precompute antenna positions
ant_pos = antenna_locations(1:num_ports, 1:3);  % [Nports x 3]

% For each channel (tx,rx), compute total path length to each voxel:
% L(p) = ||r_voxel - r_tx|| + ||r_voxel - r_rx||
% Then coherent sum over frequencies with the phase term exp(-1j * k * L)
% and weight by the complex S-parameter value at each frequency.

Npts = size(pts,1);
img_full_vals  = zeros(Npts,1);
img_empty_vals = zeros(Npts,1);

for ch = 1:Nch
    tx = channels(ch,1);
    rx = channels(ch,2);
    rt = ant_pos(tx,:);       % [1x3]
    rr = ant_pos(rx,:);       % [1x3]

    % Distances to all voxels (vectorized)
    dt = sqrt(sum((pts - rt).^2, 2));   % [Npts x 1]
    dr = sqrt(sum((pts - rr).^2, 2));   % [Npts x 1]
    L  = dt + dr;                        % [Npts x 1]

    % Phase term for all (freq, voxel): exp(-j * k(f) * L(voxel))
    % Compute as outer product: k [Nf x 1] and L' [1 x Npts]
    phase = exp(-1j * (k * L.'));        % [Nf x Npts]

    % Apply S-parameters at this channel across frequency and sum coherently
    img_full_vals  = img_full_vals  + (phase.' * scan_full(:,ch));   % [Npts x 1]
    img_empty_vals = img_empty_vals + (phase.' * scan_empty(:,ch));  % [Npts x 1]
end

% Magnitude image
img_full_vals  = abs(img_full_vals);
img_empty_vals = abs(img_empty_vals);

% Place values back onto a 3D grid (only at masked indices)
for n = 1:numel(ix_all)
    img_full_grid(ix_all(n),  iy_all(n),  iz_all(n))  = img_full_vals(n);
    img_empty_grid(ix_all(n), iy_all(n), iz_all(n))   = img_empty_vals(n);
end

% Shared normalization for fair SSIM
global_max = max([max(img_full_vals), max(img_empty_vals)]);
if global_max == 0, global_max = 1; end

%% ------------------ SSIM per z-slice + GIFs -------------------
if isempty(which('ssim'))
    error('The "ssim" function requires the Image Processing Toolbox.');
end

ssim_vals = zeros(numel(z_vals),1);
firstFrame = true;

% Find the z indices closest to requested z_vals
[~, iz_req] = arrayfun(@(z0) min(abs(zv - z0)), z_vals);

for kk = 1:numel(z_vals)
    iz = iz_req(kk);
    z0 = zv(iz);

    % Extract slices in X-Y at fixed Z index
    A = squeeze(img_full_grid(:,:,iz)).';    % [ny x nx] -> transpose so x on horizontal axis, y vertical
    B = squeeze(img_empty_grid(:,:,iz)).';

    % Normalize to [0,1] using global max
    A = A / global_max;  B = B / global_max;

    % Replace NaNs (outside ellipse) with 0 for SSIM math; use alpha to visualize holes
    A_for_ssim = A;  B_for_ssim = B;
    A_for_ssim(~isfinite(A_for_ssim)) = 0;
    B_for_ssim(~isfinite(B_for_ssim)) = 0;

    [val, ssim_map] = ssim(A_for_ssim, B_for_ssim);
    ssim_vals(kk) = val;

    % ---------- SSIM map GIF ----------
    f1 = figure('Visible','off');
    imagesc(xv, yv, ssim_map); axis image; set(gca,'YDir','normal');
    title(sprintf('SSIM map at z = %.2f m (SSIM = %.3f)', z0, val));
    xlabel('X (m)'); ylabel('Y (m)'); colorbar; caxis([0 1]);
    frame = getframe(f1); [imind, cm] = rgb2ind(frame2im(frame), 256);
    if firstFrame, imwrite(imind, cm, gif_ssim_maps, 'gif', 'LoopCount', inf, 'DelayTime', 0.25);
    else,         imwrite(imind, cm, gif_ssim_maps, 'gif', 'WriteMode', 'append', 'DelayTime', 0.25); end
    close(f1);

    % ---------- FULL slice GIF ----------
    f2 = figure('Visible','off');
    imagesc(xv, yv, A); axis image; set(gca,'YDir','normal');
    title(sprintf('FULL slice (norm) at z = %.2f m', z0));
    xlabel('X (m)'); ylabel('Y (m)'); colorbar; caxis([0 1]);
    frame = getframe(f2); [imind, cm] = rgb2ind(frame2im(frame), 256);
    if firstFrame, imwrite(imind, cm, gif_full_slices, 'gif', 'LoopCount', inf, 'DelayTime', 0.25);
    else,         imwrite(imind, cm, gif_full_slices, 'gif', 'WriteMode', 'append', 'DelayTime', 0.25); end
    close(f2);

    % ---------- EMPTY slice GIF ----------
    f3 = figure('Visible','off');
    imagesc(xv, yv, B); axis image; set(gca,'YDir','normal');
    title(sprintf('EMPTY slice (norm) at z = %.2f m', z0));
    xlabel('X (m)'); ylabel('Y (m)'); colorbar; caxis([0 1]);
    frame = getframe(f3); [imind, cm] = rgb2ind(frame2im(frame), 256);
    if firstFrame, imwrite(imind, cm, gif_empty_slices, 'gif', 'LoopCount', inf, 'DelayTime', 0.25);
    else,         imwrite(imind, cm, gif_empty_slices, 'gif', 'WriteMode', 'append', 'DelayTime', 0.25); end
    close(f3);

    firstFrame = false;
end

%% ------------------ Summary: SSIM vs depth --------------------
figure;
plot(z_vals, ssim_vals, 'o-','LineWidth',1.5);
xlabel('z (m)'); ylabel('SSIM (FULL vs EMPTY)'); grid on;
title('Per-slice SSIM across depth');

fprintf('\nMean SSIM over slices: %.4f\n', mean(ssim_vals));
[bestSSIM, idxBest] = max(ssim_vals);
fprintf('Best SSIM: %.4f at z = %.3f m\n', bestSSIM, z_vals(idxBest));
