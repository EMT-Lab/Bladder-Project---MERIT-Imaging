% Pick Folder and Files
% Load the S-parameters
Folder = 'data_6ant_final';
SnP_E = 'Empty_6ant_1GHz_5GHz.s6p';
SnP_F = 'Full_6ant_1GHz_5GHz.s6p';

% Load the S-parameters
SparamsE = sparameters([Folder '/' SnP_E]);
SparamsF = sparameters([Folder '/' SnP_F]);

% Extract frequency (Hz) as a column vector (assuming both files have the same frequencies)
freq = SparamsE.Frequencies;
num_freq = length(freq);

% Extract S-parameter matrices
S_E = SparamsE.Parameters; % Complex S-parameters for file E
S_F = SparamsF.Parameters; % Complex S-parameters for file F

% signal filtered version of s-parameters
S_E_mag = 20 * log10(abs(S_E));
S_F_mag = 20 * log10(abs(S_F));

threshold_dB = -50; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
low_mag_mask_E = S_E_mag < threshold_dB;
low_mag_mask_F = S_F_mag < threshold_dB;

S_E(low_mag_mask_E) = 0;
S_F(low_mag_mask_F) = 0;

% Get the number of ports
num_ports = size(S_E, 1);
num_channels = num_ports^2; % Total channels (S11, S12, ..., Snn)

% Prepare channel names matrix (Nx2)
channels = zeros(num_channels, 2);
index = 1;
for row = 1:num_ports
    for col = 1:num_ports
        channels(index, 1) = row;
        channels(index, 2) = col;
        index = index + 1;
    end
end

% Prepare S-parameter data for E and F (real and imaginary combined)
S_combined_E = zeros(num_freq, num_channels);
S_combined_F = zeros(num_freq, num_channels);

index = 1;
for row = 1:num_ports
    for col = 1:num_ports
        S_ij_E = squeeze(S_E(row, col, :));  % Extract specific S-parameter for E
        S_ij_F = squeeze(S_F(row, col, :));  % Extract specific S-parameter for F
        S_combined_E(:, index) = S_ij_E;  % Store complex values for E
        S_combined_F(:, index) = S_ij_F;  % Store complex values for F
        index = index + 1;
    end
end

% Save frequencies as CSV
writematrix(freq, [Folder '/frequencies.csv']);

% Save channel names as CSV
writematrix(channels, [Folder '/channels.csv']);

% Save S-parameter data for both files (real and imaginary combined)
writematrix(S_combined_F - S_combined_E, [Folder '/scan.csv']);
writematrix(S_combined_F, [Folder '/scan1.csv']);
writematrix(S_combined_E, [Folder '/scan2.csv']);