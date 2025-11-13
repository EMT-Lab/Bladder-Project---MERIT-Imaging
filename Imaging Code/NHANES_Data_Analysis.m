data = xptread('BMX_I.XPT');

Waist_data = table2array(data(:,19));
Waist_mean = mean(Waist_data,'omitnan');

SAD_data = table2array(data(:,25));
SAD_mean = mean(SAD_data,'omitnan');
