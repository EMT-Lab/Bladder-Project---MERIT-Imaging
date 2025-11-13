
%%%%%%%% EXAMPLE FROM GITHUB %%%%%%%%%%%

% Clone a repository from GitHub with its URL
%repoURL = 'https://github.com/EMFMed/MERIT.git'; % Replace with the actual repository URL
%system(['git clone ' repoURL]);

%% Load sample data_high_loc3_no_refl (antenna locations, frequencies and signals)
% frequencies = dlmread('data_high_loc3/frequencies.csv');
% antenna_locations = dlmread('data_high_loc3/antenna_locations.csv');
% channel_names = dlmread('data_high_loc3/channel_names.csv');
% 
% load('data_high_loc3/scans.mat')
% 
% signals = scan2 - scan1;
Folder = 'Array';
frequencies = readmatrix([Folder '/frequencies.csv']);
antenna_locations = readmatrix([Folder '/antenna_locations.csv']);
channel_names = readmatrix([Folder '/channels.csv']);
signals = readmatrix([Folder '/scan.csv']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

normalize_intensity = true;
% parameters for elliptical domain
b = 0.1705;  % Semi-major axis of the ellipse (x direction)
a = 0.1064;  % Semi-minor axis of the ellipse (y direction)
resolution = 2e-3;  % Resolution
height = 0.15;  % Height of the elliptical cylinder


cutY = 0.140;
% create a grid of points within the elliptical domain
% [X, Y, Z] = meshgrid(-a:resolution:a, -b:resolution:b, 0:resolution:height);
% create grid of points for cut model
[X, Y, Z] = meshgrid(0:resolution:a, -cutY:resolution:cutY, 0:resolution:height);
% filter points to be within the elliptical radius
mask = (X/a).^2 + (Y/b).^2 <= 1;
X = X(mask);
Y = Y(mask);
Z = Z(mask);

% combine points into single array
points = [X(:), Y(:), Z(:)];

%% Calculate delays for synthetic focusing
delays = merit.beamform.get_delays(channel_names, antenna_locations, ...
  'relative_permittivity', 19.72);%19.72-1i*4.55

%% Perform imaging
img = abs(merit.beamform(signals, frequencies, points, delays, ...
        merit.beamformers.DAS)); %DMAS more advanced

%% Plot
% B_wall_top_p = xlsread("Full_Bladder_Wall_Points_topdown_80mm.xlsx").*0.001;
% B_wall_top_x = B_wall_top_p(:,1);
% B_wall_top_x = interp(B_wall_top_x,3);
% B_wall_top_y = B_wall_top_p(:,2);
% B_wall_top_y = interp(B_wall_top_y,3);
% z value for the slice
overall_max = 0;

% Calculate the maximum intensity among all slides
for z_slice_value = 0:0.01:0.15
    slice_axes_x = unique(points(:, 1));
    slice_axes_y = unique(points(:, 2));
    slice_axes_z = unique(points(:, 3));
    axes_ = {slice_axes_x, slice_axes_y, slice_axes_z};
    % slice
    im_slice = merit.visualize.get_slice(img, points, {slice_axes_x, slice_axes_y, slice_axes_z}, 'z', z_slice_value);
    max_slice = max(img);
    if max_slice > overall_max
        overall_max = max_slice; % Update overall maximum intensity
    end

end

for z_slice_value = 0:0.01:0.15
    % z_slice_value = 0.01;

    % % indices of the points at the desired z value
    % slice_indices = abs(Z - z_slice_value) < resolution/2;
    % slice_points = points(slice_indices, :);

    % corresponding axes for the slice
    slice_axes_x = unique(points(:, 1));
    slice_axes_y = unique(points(:, 2));
    slice_axes_z = unique(points(:, 3));
    axes_ = {slice_axes_x, slice_axes_y, slice_axes_z};
    % slice
    im_slice = merit.visualize.get_slice(img, points, {slice_axes_x, slice_axes_y, slice_axes_z}, 'z', z_slice_value);
    if normalize_intensity
        im_slice = im_slice / overall_max; % Normalize the intensity
    end
    % display
    figure;
    imagesc(slice_axes_x, slice_axes_y, im_slice');
    % hold on;
    % plot(B_wall_top_x, B_wall_top_y, 'k-', 'LineWidth', 2);
    % hold off;
    xlabel('X-axis');
    ylabel('Y-axis');
    title(['Slice at Z = ', num2str(z_slice_value)]);
    axis image
    % set the limits to match the elliptical domain dimensions with added padding
    padding = 0.01; % Adjust the padding as needed
    xlim([0 - padding, a + padding]);
    ylim([-cutY - padding, cutY + padding]);

    % set the aspect ratio to match the ellipse
    pbaspect([a b 1]);

    % antenna positions equally spaced around the ellipse
    % num_antennas = 8; 
    % theta = linspace(0, 2*pi, num_antennas + 1);
    % theta(end) = []; % remove last element to avoid overlap
    % x_locations = a * cos(theta);
    % y_locations = b * sin(theta);
    % z_locations = zeros(size(x_locations)); % z=0 plane

    % get current axes handle
    ax = gca;

    % adjust axis limits and ticks
    x_limits = get(ax, 'XLim');
    y_limits = get(ax, 'YLim');
    ticks_x = linspace(x_limits(1), x_limits(2), 5);
    ticks_y = linspace(y_limits(1), y_limits(2), 5);

    % set the ticks for both axes
    set(ax, 'XTick', ticks_x, 'YTick', ticks_y);
    colorbar;
    if normalize_intensity
        clim([0 1]);
    end
    %clim([(max(img)/150) (max(img)/2.2)]);
    % hold on;
    
    % Save the current figure as a PNG file
    saveas(gcf, fullfile(Folder, ['slice_' num2str(round(z_slice_value*100)) '.png']));

    % Initialize GIF filename
    gifFilename = 'slice_images.gif';

    % Capture the plot as an image
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    
    % Write to the GIF file
    if z_slice_value == 0
        imwrite(imind, cm, gifFilename, 'gif', 'LoopCount', inf, 'DelayTime', 1.0);
    else
        imwrite(imind, cm, gifFilename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.25);
    end
    
    %close(gcf); % Close the figure after saving
  
end


% plot antenna locations
% scatter(x_locations, y_locations, 'r', 'filled');
% 
% % labels
% for i = 1:length(x_locations)
%     % adjust label positioning
%     text(x_locations(i), y_locations(i), sprintf('A%d', i), ...
%         'Color', 'red', 'FontSize', 12, 'FontWeight', 'bold', ...
%         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
%         'Margin', 2);
% end
maxValues = zeros(size(slice_axes_z));
xMax = zeros(size(slice_axes_z));
yMax = zeros(size(slice_axes_z));

points_img = [points,img];

for i = 1:length(slice_axes_z)
    z = slice_axes_z(i);

    points_img_atZ = points_img(abs(points_img(:,3) - z) < 1e-6,:);

    [maxVal, idx] = max(points_img_atZ(:,4));

    maxValues(i) = maxVal;
    xMax(i) = points_img_atZ(idx, 1);
    yMax(i) = points_img_atZ(idx, 2);
end

resultTable = table(slice_axes_z, maxValues,xMax,yMax);


%% Side View Plot of Max Intensity
%overlay ellipse
center = [0, 0];    % ellipse center (x0, y0)
a = 0.1064;                % semi-major axis
b = 0.1705;                % semi-minor axis
theta = linspace(0, 2*pi, 100);  % angle range
z_layer = 0.09;        % the Z height at which to plot the ellipse
% Parametric equations
x_ellipse = center(1) + a*cos(theta);
x_ellipse2 = center(1) + (a-0.02)*cos(theta);
y_ellipse = center(2) + b*sin(theta);
y_ellipse2 = center(2) + (b-0.02)*sin(theta);
z_ellipse = z_layer * ones(size(theta));  % keep Z constant

% B_wall_p = xlsread("Full_Bladder_Wall_Points.xlsx").*0.001;
% B_wall_x = B_wall_p(:,1);
% B_wall_x = interp(B_wall_x,3);
% B_wall_z = B_wall_p(:,2);
% B_wall_z = interp(B_wall_z,3);
% 
% Urine_p = xlsread("Full_Urine_Points.xlsx").*0.001;
% Urine_x = Urine_p(:,1);
% Urine_x = interp(Urine_x,3);
% Urine_z = Urine_p(:,2);
% Urine_z = interp(Urine_z,3);

figure();
% scatter(xMax,slice_axes_z,50,maxValues,'filled');
scatter(slice_axes_z,xMax,100,maxValues,'filled');
hold on;
% scatter(antenna_locations(:,1),antenna_locations(:,3),'k','square','filled');
% plot(xMax,slice_axes_z,'b--','LineWidth',1.5);
% plot(B_wall_x, B_wall_z, 'r-', 'LineWidth', 2);
% plot(Urine_x, Urine_z, 'y-', 'LineWidth', 2);
% xlabel('X');
% ylabel('Z');
scatter(antenna_locations(:,3), antenna_locations(:,1),300,'k','square','filled');
plot(slice_axes_z, xMax,'b--','LineWidth',1.5);
%plot(B_wall_z, B_wall_x, 'r-', 'LineWidth', 2);
%plot(Urine_z, Urine_x, 'y-', 'LineWidth', 2);
xlabel('Z');
ylabel('X');
% title('Max Intensity at Each Z Layer');
axis('image');
grid on;
fontsize(20,'points');
colorbar;
legend('Max Intensities','Antenna Locations','','Bladder Wall','Urine');
hold off;