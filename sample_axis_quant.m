%Evan Underhill
%Quantification code for paper "Control of gastruloid patterning 
%and morphogenesis by the Erk and Akt signaling pathways" Development 2023

%This code takes as input (1) a .tif stack with n images of stained gastruloids, 
%(2) a folder with n .txt files, each containing the axis coordinates of an
%individual gastruloid, (3) a folder with n .txt files, each containing the
%edges coordinates of an individual gastruloid. 200 equally spaced minor
%axis curves are generated perpendicular to the specified major axis, the
%mean pixel intensity is calculated for each minor axis curve, and this mean
%pixel intensity is plotted as a function of major axis length. Must
%specify pixel_micron_ratio for particular dataset. 

clc 
clear all
hold off

%% Quantify major axis profiles for individual gastruloids
%determine total number of files
axis_file_dir = '/axis_directory';
axis_files = dir(fullfile(axis_file_dir, '*.txt'));
nfiles = length(axis_files)

%specify total_num_orth for number of major axis curves
total_num_orth = 200;

%MODIFY PIXEL/MICRON RATIO FOR YOUR IMAGE SET
pixel_micron_ratio = 922/1293.41;

%initialize matrix for aggregating intensity profiles
WT_aggregate_001 = zeros(nfiles, total_num_orth);

for a = 1:nfiles

    %indexing works for file format '..._001.txt' 
    current_num = axis_files(a).name(end-6:end-4)
    
    axis = readmatrix(strcat('/axis_directory/filename_',num2str(current_num), '.txt'));
    edges = readmatrix(strcat('/edges_directory/filename_',num2str(current_num), '.txt'));
    staining_tif = imread('/tif_directory/filename.tif',str2num(current_num));
    
    xaxis = axis(:,1);
    yaxis = axis(:,2);
    xedges = edges(:,1);
    yedges = edges(:,2);
    
    %close the edges polygon
    xedges(end + 1) = xedges(1);
    yedges(end + 1) = yedges(1);
    
    %convert microns to pixels
    xaxis = xaxis .* pixel_micron_ratio;
    yaxis = yaxis .* pixel_micron_ratio;
    xedges = xedges .* pixel_micron_ratio;
    yedges = yedges .* pixel_micron_ratio;
    
    %define x values for orthogonal lines, bounded by the min and max xedge
    x = min(xedges):.1:max(xedges);
    
    %determine length of each major axis segment and total length of axis
    lengths = zeros(1, length(xaxis)-1);
    
    for i = 1:(length(xaxis)-1)
        lengths(i) = sqrt((xaxis(i)-xaxis(i+1))^2 + (yaxis(i)-yaxis(i+1))^2);
    end
    
    sum_length = sum(lengths);
    
    %Determine number of orthogonal lines for each segment (proportionally allocated)
    num_orth = zeros(1, length(xaxis)-1);
    
    for i = 1:(length(xaxis)-1)
        num_orth(i) = round(total_num_orth * lengths(i) / sum_length);
    end
    
    %Account for rounding to nearest integer
    if sum(num_orth) > total_num_orth
        num_orth(1) = num_orth(1) - 1;
    end
    
    if sum(num_orth) < total_num_orth
        num_orth(1) = num_orth(1) + 1;
    end
    
    %initialize array for mean intensity of each minor axis curve
    intensity_mean = zeros(1, total_num_orth);
    
    count = 1;
    orth_count = 0;
    
    %calculate intensity along minor axis curves
    for i = 1:(length(xaxis)-1)  
    
        x0 = xaxis(i);
        y0 = yaxis(i);

           for j = 0:(num_orth(i)-1)
           
             x_temp = x0 + j * (1/num_orth(i))*(xaxis(i+1) - xaxis(i));
             y_temp = y0 + j * (1/num_orth(i))*(yaxis(i+1) - yaxis(i));
             m = -(xaxis(i+1) - xaxis(i))/(yaxis(i+1) - yaxis(i));
             
             %if orthogonal lines are vertical, set m to arbitrary high value
             if m == -Inf
                 m = 100;
             end
             
             y = m*(x - x_temp) + y_temp; 
             in = inpolygon(x,y,xedges,yedges);
             
             %finds where 'in' switches from 0 to 1 or 1 to 0
             f = find(diff([false,in==1,false])~=0);
             
             %find the index where x crosses the axis line
             [minValue,closestIndex] = min(abs(x-x_temp));
             
             %set any ones outside of the desired segment equal to zero
             for p = 1:(length(f)-1)
                 if (closestIndex > f(p)) && (closestIndex < f(p+1))
                     in(1:f(p)) = 0;
                     in(f(p+1):end) = 0;
                 end
             end
                   
             xinside = x(in);
             yinside = y(in);
             
             %initialize array for intensity along minor axis curves
             intensity = zeros(1, length(x(in)));
                      
             for k = 1:length(x(in))
                 %y is the row, x is the column
                 intensity(k) = staining_tif(ceil(yinside(k)),ceil(xinside(k)));
             end
             
             %assign mean intensity of minor axis curve to aggregate array
             intensity_mean(count) = mean(intensity);
    
             orth_count = orth_count + 1;
             count = count + 1;  
    
          end
    end
    
    normalized_x = linspace(0, 1, total_num_orth);
    
    %plot major axis profiles for individual gastruloids
    figure
    plot(normalized_x, intensity_mean, 'k', 'LineWidth',2)
    xlabel('major axis (normalized)')
    ylabel('intensity (a.u.)')
    set(gca,'FontSize',15)
    
    %combine major axis profiles into matrix
    WT_aggregate_001(a, :) = intensity_mean;

end


%% Aggregate and Plot Data

WT_exp_001 = mean(WT_aggregate_001);
WT_std_001 = std(WT_aggregate_001);

figure
plot(normalized_x, WT_exp_001, 'k', 'LineWidth',4)
hold on
plot(normalized_x, WT_exp_001 - WT_std_001, 'k', 'LineWidth',2)
plot(normalized_x, WT_exp_001 + WT_std_001, 'k', 'LineWidth',2)
xlabel('major axis (normalized)')
ylabel('intensity (a.u.)')
set(gca,'FontSize',15)