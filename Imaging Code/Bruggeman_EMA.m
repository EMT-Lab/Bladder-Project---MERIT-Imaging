%% Two Mediums

MD = 1090; %muscle density
FD = 911; %fat density
FMR = 0.4; %fat to muscle mass ratio
Mixed_Volume = 0.00316; %volume of mixed fat and muscle area
MFMV = Mixed_Volume; %same as above. easier for equations

MM = MFMV*((FD+FMR*MD)/(MD*FD))^(-1); %muscle mass in mixed area
FM = FMR*MM; %fat mass in mixed area

MVF = MM/(MD*MFMV); %muscle volume fraction
FVF = FM/(FD*MFMV); %fat volume fraction

cm = MVF;
cd = FVF;
eo = 8.85418782*10^(-12);
muscle = load("Muscle_Perm.tab");
muscle_c = load("Muscle_Cond.tab");
fat = load("Fat_Perm.tab");
fat_c = load("Fat_Cond.tab");

freq = muscle(:,1);
ang_freq = freq.*2*pi;

muscle_perm = muscle(:,2);
muscle_imag = muscle_c(:,2)./(ang_freq.*eo); %cond loaded in and transformed to imaginary part of permittivity

fat_perm = fat(:,2);
fat_imag = fat_c(:,2)./(ang_freq.*eo);

em = muscle_perm - 1i.*muscle_imag;
ed = fat_perm - 1i.*fat_imag;

Hb = (3.*cd - 1).*ed + (3.*cm - 1).*em;

E_eff = (Hb + sqrt(Hb.^(2) + 8.*em.*ed))./4;
E_eff_real = real(E_eff);
E_eff_cond = imag(E_eff).*(-1.*ang_freq*eo);

E_perm_out = [freq,E_eff_real];
E_cond_out = [freq,E_eff_cond];

%% Four Mediums
% 
% % Define layer thicknesses (in meters)
% t_muscle = 0.01;  % Muscle layer thickness
% t_fat = 0.0135;     % Fat layer thickness
% t_skin = 0.002;   % Skin layer thickness
% t_mfm = 0.02025;     % Muscle-fat mixture layer thickness
% 
% % Define densities (kg/m³)
% MD = 1090;  % Muscle density
% FD = 911;   % Fat density
% SD = 1109;  % Skin density
% MFMD = 1032; % Muscle-Fat Mixture density (FtoMRatio = 0.4)
% 
% % Compute mass for each layer
% M_muscle = MD * t_muscle;
% M_fat = FD * t_fat;
% M_skin = SD * t_skin;
% M_mfm = MFMD * t_mfm;
% 
% % Compute total mass for each Bruggeman step
% M_sf = M_skin + M_fat;
% M_mmfm = M_muscle + M_mfm;
% M_total = M_sf + M_mmfm;
% 
% % Compute mass fractions
% MF_SF_skin = M_skin / M_sf;
% MF_SF_fat = M_fat / M_sf;
% MF_MMFM_muscle = M_muscle / M_mmfm;
% MF_MMFM_mfm = M_mfm / M_mmfm;
% MF_total_SF = M_sf / M_total;
% MF_total_MMFM = M_mmfm / M_total;
% 
% % Convert mass fractions to volume fractions
% VF_SF_skin = MF_SF_skin * (FD / SD); % Skin fraction in step 1
% VF_SF_fat = MF_SF_fat * (SD / FD);   % Fat fraction in step 1
% VF_MMFM_muscle = MF_MMFM_muscle * (MFMD / MD); % Muscle fraction in step 2
% VF_MMFM_mfm = MF_MMFM_mfm * (MD / MFMD);      % MFM fraction in step 2
% VF_total_SF = MF_total_SF * (M_mmfm / M_sf);  % Skin-fat mixture in step 3
% VF_total_MMFM = MF_total_MMFM * (M_sf / M_mmfm); % Muscle-MFM mixture in step 3
% 
% % Load permittivity and conductivity data
% eo = 8.85418782e-12;
% muscle = load("Muscle_Perm.tab");
% muscle_c = load("Muscle_Cond.tab");
% fat = load("Fat_Perm.tab");
% fat_c = load("Fat_Cond.tab");
% skin = load("Skin_Perm.tab");
% skin_c = load("Skin_Cond.tab");
% mfm = load("Pelvic_Mixed_Tissue_Perm.tab");
% mfm_c = load("Pelvic_Mixed_Tissue_Cond.tab");
% 
% freq = muscle(:,1);
% ang_freq = freq .* 2 * pi;
% 
% % Compute complex permittivities
% muscle_perm = muscle(:,2);
% muscle_imag = muscle_c(:,2) ./ (ang_freq .* eo);
% 
% fat_perm = fat(:,2);
% fat_imag = fat_c(:,2) ./ (ang_freq .* eo);
% 
% skin_perm = skin(:,2);
% skin_imag = skin_c(:,2) ./ (ang_freq .* eo);
% 
% mfm_perm = mfm(:,2);
% mfm_imag = mfm_c(:,2) ./ (ang_freq .* eo);
% 
% % Complex permittivities
% em = muscle_perm - 1i * muscle_imag;
% ed = fat_perm - 1i * fat_imag;
% es = skin_perm - 1i * skin_imag;
% emfm = mfm_perm - 1i * mfm_imag;
% 
% % Step 1: Compute effective permittivity for skin + fat
% Hb1 = (3*VF_SF_fat - 1).*ed + (3*VF_SF_skin - 1).*es;
% E_SF = (Hb1 + sqrt(Hb1.^2 + 8.*es.*ed)) ./ 4;
% 
% % Step 2: Compute effective permittivity for MFM + muscle
% Hb2 = (3*VF_MMFM_mfm - 1).*emfm + (3*VF_MMFM_muscle - 1).*em;
% E_MMFM = (Hb2 + sqrt(Hb2.^2 + 8.*emfm.*em)) ./ 4;
% 
% % Step 3: Compute final effective permittivity for (skin-fat) + (muscle-MFM)
% Hb3 = (3*VF_total_MMFM - 1).*E_MMFM + (3*VF_total_SF - 1).*E_SF;
% E_eff = (Hb3 + sqrt(Hb3.^2 + 8.*E_SF.*E_MMFM)) ./ 4;
% 
% % Compute real permittivity and conductivity
% E_eff_real = real(E_eff);
% E_eff_cond = -1 .* imag(E_eff) .* (ang_freq * eo);
% E_eff_imag = -1 .* imag(E_eff);
% 
% % Save results
% E_perm_out = [freq, E_eff_real];
% E_cond_out = [freq, E_eff_cond];
% E_iperm_out = [freq, E_eff_imag];
% 
% writematrix(E_perm_out, 'Effective_Permittivity.csv');
% writematrix(E_cond_out, 'Effective_Conductivity.csv');
% writematrix(E_iperm_out, 'Effective_IPermittivity.csv');