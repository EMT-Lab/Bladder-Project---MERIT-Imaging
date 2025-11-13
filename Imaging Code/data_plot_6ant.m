%% Load
S_Full = sparameters('data_6ant_final\Full_6ant_1GHz_5GHz.s6p');
S_Empty = sparameters('data_6ant_final\Empty_6ant_1GHz_5GHz.s6p');

Freq = S_Full.Frequencies;
threshold_dB = -50; %%%%%%%%%%%%%%%%%%%%%%%%
matched_threshold = -10;
%% Diffs

CatF = [];
CatE = [];
Max_Diffs = [];
Max_Diffs_Lin = [];
Max_Diffs_Lin_filt = [];
Max_Diffs_Lin_fullfilt = [];
Max_Diffs_Lin_topfilt = [];

for k = 1:6
    for l = 1:6
        if l >= k % >= for all, > for just transmission, == for just reflection 
            SF_dB = 20*log10(abs(squeeze(S_Full.Parameters(k,l,:))));
            SF_dB_Refl1 = 20*log10(abs(squeeze(S_Full.Parameters(k,k,:))));
            SF_dB_Refl2 = 20*log10(abs(squeeze(S_Full.Parameters(l,l,:))));
            SE_dB = 20*log10(abs(squeeze(S_Empty.Parameters(k,l,:))));
            SE_dB_Refl1 = 20*log10(abs(squeeze(S_Empty.Parameters(k,k,:))));
            SE_dB_Refl2 = 20*log10(abs(squeeze(S_Empty.Parameters(l,l,:))));
            low_mag_mask_F = SF_dB < threshold_dB;
            low_mag_mask_E = SE_dB < threshold_dB;
            SF_dB_filt = SF_dB;
            SE_dB_filt = SE_dB;
            SF_dB_filt(low_mag_mask_F) = threshold_dB;
            SE_dB_filt(low_mag_mask_E) = threshold_dB;
            Diff = abs(SF_dB_filt - SE_dB_filt);
            for f = 1:length(Freq)
                if ((SE_dB_Refl1(f) > matched_threshold) || (SF_dB_Refl1(f) > matched_threshold) || (SE_dB_Refl2(f) > matched_threshold) || (SF_dB_Refl2(f) > matched_threshold))
                    Diff(f) = 0;
                end
            end
            [MD,Idx] = max(Diff);
            DFreq = Freq(Idx);
            SL = SE_dB(Idx);
            CH = [num2str(k),num2str(l)];
            CH = str2num(CH);
            Max_Diffs = [Max_Diffs;[MD,SL,DFreq,CH]];


            SF = (abs(squeeze(S_Full.Parameters(k,l,:))));
            SE = (abs(squeeze(S_Empty.Parameters(k,l,:))));
            SF_filt = SF;
            SE_filt = SE;
            SF_topfilt = SF;
            SE_topfilt = SE;
            SF_filt(low_mag_mask_F) = 10^(threshold_dB/20);
            SE_filt(low_mag_mask_E) = 10^(threshold_dB/20);
            SF_fullfilt = SF_filt;
            SE_fullfilt = SE_filt;
            % Diff = abs(SF - SE); % lin mag diff
            Diff = abs(SF - SE).*100./SE; %percent lin diff
            Diff_filt = abs(SF_filt - SE_filt).*100./SE_filt;
            for f = 1:length(Freq)
                if ((SE_dB_Refl1(f) > matched_threshold) || (SF_dB_Refl1(f) > matched_threshold) || (SE_dB_Refl2(f) > matched_threshold) || (SF_dB_Refl2(f) > matched_threshold))
                    SF_fullfilt(f) = 0;
                    SE_fullfilt(f) = 0;
                    SF_topfilt(f) = 0;
                    SE_topfilt(f) = 0;
                end
            end
            Diff_fullfilt = abs(SF_fullfilt - SE_fullfilt).*100./SE_fullfilt;
            Diff_topfilt = abs(SF_topfilt - SE_topfilt).*100./SE_topfilt;

            [MD,Idx] = max(Diff);
            [MD_filt,Idx_filt] = max(Diff_filt);
            [MD_fullfilt,Idx_fullfilt] = max(Diff_fullfilt);
            [MD_topfilt,Idx_topfilt] = max(Diff_topfilt);
            DFreq = Freq(Idx);
            DFreq_filt = Freq(Idx_filt);
            DFreq_fullfilt = Freq(Idx_fullfilt);
            DFreq_topfilt = Freq(Idx_topfilt);
            SL = SE_dB(Idx);
            SL_filt = SE_dB(Idx_filt);
            SL_fullfilt = SE_dB(Idx_fullfilt);
            SL_topfilt = SE_dB(Idx_topfilt);
            Max_Diffs_Lin = [Max_Diffs_Lin;[MD,SL,DFreq,CH]];
            Max_Diffs_Lin_filt = [Max_Diffs_Lin_filt;[MD_filt,SL_filt,DFreq_filt,CH]];
            Max_Diffs_Lin_fullfilt = [Max_Diffs_Lin_fullfilt;[MD_fullfilt,SL_fullfilt,DFreq_fullfilt,CH]];
            Max_Diffs_Lin_topfilt = [Max_Diffs_Lin_topfilt;[MD_topfilt,SL_topfilt,DFreq_topfilt,CH]];

            % if k ~= l
            % 
            %         SF_comp = squeeze(S_Full.Parameters(k,l,:));
            %         SE_comp = squeeze(S_Empty.Parameters(k,l,:));
            %         Diff_Comp = 
        end
    end
end

%% dB Diff
figure;
bar(Max_Diffs(:,1));
xlabel('Channel');
ylabel('Maximum Difference [dB]');
title('Max Difference per Channel with Frequency and Signal Level');
ylim([0, max(Max_Diffs(:,1)) * 1.2]);

% Set x-axis ticks to actual channel numbers
xticks(1:length(Max_Diffs));
xticklabels(string(Max_Diffs(:,4)));  % Use actual channel numbers as labels

% Add text labels above each bar
for i = 1:length(Max_Diffs)
    label = sprintf('%.1f dB\n%.2f GHz', Max_Diffs(i,2), Max_Diffs(i,3)/1E9);
    text(i, Max_Diffs(i,1) + 0.01, label, ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
         'FontSize', 8);
end
%% Remove 16 and 34
Max_Diffs_Lin_fullfilt([6,13],:) = [];
Max_Diffs_Lin_filt([6,13],:) = [];
Max_Diffs_Lin_topfilt([6,13],:) = [];

%% Lin Diff
figure;
% bar(Max_Diffs_Lin(:,1));
% hold on;
% bar([Max_Diffs_Lin_topfilt(:,1),Max_Diffs_Lin_filt(:,1),Max_Diffs_Lin_fullfilt(:,1)]);
% hold on;
% bar(Max_Diffs_Lin_filt(:,1));
bar(Max_Diffs_Lin_fullfilt(:,1));
hold on;
% for i = 1:length(Max_Diffs_Lin_topfilt)
%     if Max_Diffs_Lin_topfilt(i,1) > 20
%         Max_Diffs_Lin_topfilt(i,1) = 20;
%     end
%     if Max_Diffs_Lin_filt(i,1) > 20
%         Max_Diffs_Lin_filt(i,1) = 20;
%     end
% end
scatter(1:length(Max_Diffs_Lin_filt),Max_Diffs_Lin_filt(:,1),600,'filled','square');
scatter(1:length(Max_Diffs_Lin_topfilt),Max_Diffs_Lin_topfilt(:,1),480,'filled','diamond');
xlabel('Channel');
ylabel('Maximum Difference [Linear Percent]');
% title('Max Difference per Channel with Frequency and Signal Level');
% ylim([0, max([max(Max_Diffs_Lin_topfilt(:,1)),max(Max_Diffs_Lin_filt(:,1)),max(Max_Diffs_Lin_fullfilt(:,1))]) * 1.1]);
axis tight;
ylim([0, 20]);
% Set x-axis ticks to actual channel numbers
xticks(1:length(Max_Diffs_Lin_fullfilt));
xticklabels(string(Max_Diffs_Lin_fullfilt(:,4)));  % Use actual channel numbers as labels
fontsize(25,'points');
% Add text labels above each bar
for i = 1:length(Max_Diffs_Lin_fullfilt)
    % label = sprintf('%.1f dB\n%.2f GHz', Max_Diffs_Lin(i,2), Max_Diffs_Lin(i,3)/1E9);
    % text(i, min([Max_Diffs_Lin(i,1),max([max(Max_Diffs_Lin_filt(:,1)),max(Max_Diffs_Lin_fullfilt(:,1))])*1.05]) + 0.0001, label, ...
    %      'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
    %      'FontSize', 8);
    % label_filt = sprintf('%.1f dB\n%.2f GHz', Max_Diffs_Lin_filt(i,2), Max_Diffs_Lin_filt(i,3)/1E9);
    % text(i, Max_Diffs_Lin_filt(i,1) + 0.0001, label_filt, ...
    %      'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
    %      'FontSize', 9);
    label_fullfilt = sprintf('%.1f \n%.2f ', Max_Diffs_Lin_fullfilt(i,2), Max_Diffs_Lin_fullfilt(i,3)/1E9);
    text(i, 21, label_fullfilt, ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
         'FontSize', 20);
    label_signal_freq = sprintf('Signal(dB)\nFrequency(GHz)');
    text(-0.7, 21, label_signal_freq, ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
         'FontSize', 18);
    % label_topfilt = sprintf('%.1f dB\n%.2f GHz', Max_Diffs_Lin_topfilt(i,2), Max_Diffs_Lin_topfilt(i,3)/1E9);
    % text(i-0.5, min([Max_Diffs_Lin_topfilt(i,1),max([max(Max_Diffs_Lin_filt(:,1)),max(Max_Diffs_Lin_fullfilt(:,1))]*1.12)]) + 0.0001, label_topfilt, ...
    %      'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
    %      'FontSize', 9);
end
lgd = legend('Full Filtered','Floor Filtered','Match Filtered','FontSize',16);
lgd.IconColumnWidth = 7;
% hold off;
% alpha(0.1);
%% MIC
Sik_F = [];
Sjk_F = [];
Sik_E = [];
Sjk_E = [];
MIC_D = ones(6);

for i = 1:6
    for j = 1:6
        if j>i
            for k = 1:6
                if ((k ~= i) && (k ~= j))
                    Sik_Ft = squeeze(S_Full.Parameters(i,k,:));
                    Sjk_Ft = squeeze(S_Full.Parameters(j,k,:));
                    Sik_Et = squeeze(S_Empty.Parameters(i,k,:));
                    Sjk_Et = squeeze(S_Empty.Parameters(j,k,:));
                    Sik_Ft_dB = 20*log10(abs(Sik_Ft));
                    Sjk_Ft_dB = 20*log10(abs(Sjk_Ft));
                    Sik_Et_dB = 20*log10(abs(Sik_Et));
                    Sjk_Et_dB = 20*log10(abs(Sjk_Et));
                    low_mag_mask_Sik_Ft = Sik_Ft_dB < threshold_dB;
                    low_mag_mask_Sjk_Ft = Sjk_Ft_dB < threshold_dB;
                    low_mag_mask_Sik_Et = Sik_Et_dB < threshold_dB;
                    low_mag_mask_Sjk_Et = Sjk_Et_dB < threshold_dB;
                    Sik_Ft(low_mag_mask_Sik_Ft) = 0;
                    Sjk_Ft(low_mag_mask_Sjk_Ft) = 0;
                    Sik_Et(low_mag_mask_Sik_Et) = 0;
                    Sjk_Et(low_mag_mask_Sjk_Et) = 0;

                    SF_dB_Refl1 = 20*log10(abs(squeeze(S_Full.Parameters(i,i,:))));
                    SF_dB_Refl2 = 20*log10(abs(squeeze(S_Full.Parameters(j,j,:))));
                    SF_dB_Refl3 = 20*log10(abs(squeeze(S_Full.Parameters(k,k,:))));
                    SE_dB_Refl1 = 20*log10(abs(squeeze(S_Empty.Parameters(i,i,:))));
                    SE_dB_Refl2 = 20*log10(abs(squeeze(S_Empty.Parameters(j,j,:))));
                    SE_dB_Refl3 = 20*log10(abs(squeeze(S_Empty.Parameters(k,k,:))));
                    for f = 1:length(Freq)
                        if ((SE_dB_Refl1(f) > matched_threshold) || (SF_dB_Refl1(f) > matched_threshold) || (SE_dB_Refl3(f) > matched_threshold) || (SF_dB_Refl3(f) > matched_threshold))
                            Sik_Ft(f) = 0;
                            Sik_Et(f) = 0;
                        end
                        if ((SE_dB_Refl2(f) > matched_threshold) || (SF_dB_Refl2(f) > matched_threshold) || (SE_dB_Refl3(f) > matched_threshold) || (SF_dB_Refl3(f) > matched_threshold))
                            Sjk_Ft(f) = 0;
                            Sjk_Et(f) = 0;
                        end
                    end
                    Sik_F = [Sik_F;Sik_Ft];
                    Sjk_F = [Sjk_F;Sjk_Ft];
                    Sik_E = [Sik_E;Sik_Et];
                    Sjk_E = [Sjk_E;Sjk_Et];
                end
            end
            diffVal1 = Sik_F - Sik_E;
            diffVal2 = Sjk_F - Sjk_E;
            
            numeratorD = abs(sum(diffVal1.*conj(diffVal2)));
            denominatorD = sqrt(sum(abs(diffVal1).^2).*sum(abs(diffVal2).^2));
            MIC_D(i,j) = numeratorD/denominatorD;
            MIC_D(j,i) = numeratorD/denominatorD;
            Sik_F = [];
            Sjk_F = [];
            Sik_E = [];
            Sjk_E = [];
        end
    end
end
% for i = 0.1:0.1:0.9
    figure;
    colormap(flipud(sky));
    imagesc(1:6, 1:6, MIC_D(:,:));
    c = colorbar;
    clim([0 0.5]);
    c.Ticks = [0 0.1 0.2 0.3 0.4 0.5];
    title('6 Antenna MIC Matrix - Difference');
    xlabel('TX');
    ylabel('RX');
    fontsize(20,'points');
    axis square;
% end