%% plotting difference plots between two HX datasets
% Written by Levani Zandarashvili, Sep/12/2018
% Code has been tested on MATLAB R2018a. It does not require any additional
% libraries.
%
% This code allows to plot the difference in HX rates for each peptide.
% Input format should be HDExaminer's .csv file.
% the program will take two files as inputs and plot a diff-plot based on
% Deut% column.
%
% The code only works if the .csv files have matching peptides only. In
% other words in they were processed using the same peptide pools. It is
% possible to modify this code to work with non-matching peptide pool
% files, but it will need an additional modification.

%% Initiate the program
% clear all previous data, close all previos figures, clear MATLAB's main
% window.
clear all
close all
clc

%% Enter the protein information
% Following prompt will allow user to enter the following information
% regarding the protein:
% Protein name
% Amino acid cut off (to limit the length of the protein)
% N-terminal amino acid offset (e.g. to correct the amino acid # due to N-terminal tag)

prompt = {'Protein Name', 'Amino Acit Cutoff', 'Offset'};
dlg_title = 'Protein information';
defaultans = {'Protein', '500', '0'};
dims = [1 50];
protein_info = (inputdlg(prompt,dlg_title,dims,defaultans));

protein = char(protein_info(1));
AAcutoff = str2double(protein_info{2});
offset = str2double(protein_info{3});

%% prompt user to select the color function
% current color function (difference_plot_color_scheme.m) takes the difference 
% in HX as an input (real float or integer) and outputs an array of three
% numbers corresponding to discrete set of colors. The function can be
% changed to have a smoother gradient of color, but the current format
% makes it easier to modify
[fname_color, pname_color] = uigetfile('*.m', 'Select Color File');
colorFunc = fname_color(1:end-2);

fig_dir = [protein '_difference plot' date];
mkdir(fig_dir)
%% Enter information regarding the samples
% User will be prompted to enter the information regarding the samples:
% Name of the first sample
% Name of the second sample
% Time of the exchange (while the assumption is that the files have same exchange times, user can enter any value)

% This information will only be used in the title of the final figure.


prompt = {'Sample 1 name', 'Sample 2 name', 'timepoint'};
dlg_title = 'Enter sample names and the timepoint';
dims = [1 75]; % Prompt windows dimensions
defaultans = {'set 1', 'set 2', '0'};
input_info = (inputdlg(prompt,dlg_title,dims,defaultans));

set_1_name = char(input_info(1));
set_2_name = char(input_info(2));
timepoint = char(input_info(3));

% Following is the name of the final figure. Will be used later in the
% program.
fig_name1 = [fig_dir, '\' set_1_name '-' set_2_name ' ' timepoint 's'];

%% Following section will allow user to select the data files.
% Files should be outputs of HDExaminer program.

%Prompt user for the first file
[fname_1, pname_1] = uigetfile('*.csv', 'Select the First File'); 
[fname_2, pname_2] = uigetfile('*.csv', 'Select the Second File'); 

%Create fully-formed filename as a string
filename_1 = fullfile(pname_1, fname_1);
filename_2 = fullfile(pname_2, fname_2);

fileID = fopen(filename_1);
dataset_1 = textscan(fileID,'%f %f %s %f %f %s %f %f %f %f %f %f %f %f %f %f %s','headerlines', 1,'delimiter', ',');
fclose(fileID);

fileID = fopen(filename_2);
dataset_2 = textscan(fileID,'%f %f %s %f %f %s %f %f %f %f %f %f %f %f %f %f %s','headerlines', 1,'delimiter', ',');
fclose(fileID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maxSeq = length(dataset_1{1});          % number of peptides in the pool
tdeltanormPercent = zeros(maxSeq,2);    % col 1: pept_num; col 2: deut_diff
tdeltanormPercent(:,1) = (1:maxSeq)';   % keep track of peptide ix (col 1 in tnormPercent

%% calculate difference in deuteration
% Following section will calculate the difference in deuteration between
% peptides which satisfy the following two conditions: their scores are
% above 0 and they don't have low confidence levels in either of the two
% files.

% For the peptides which don't satisfy either of the two conditions in any
% of files, difference in deuteration will be defines as nan.

% apply conditions
for n = 1:maxSeq
    if dataset_1{13}(n) ~= 0 && dataset_2{13}(n) ~= 0 && ... % score is not 0
            ~strcmp(char(dataset_1{17}(n)),'Low') && ~strcmp(char(dataset_2{17}(n)),'Low') && ...  % confidence is not Low
            max(dataset_1{1}(n),dataset_1{2}(n)) < AAcutoff
        tdeltanormPercent(n,2) = dataset_1{16}(n) - dataset_2{16}(n);   % calculate difference in deuteration
    else
        tdeltanormPercent(n,2) = nan; % set difference in deuteration as nan for the remaining peptides
    end
end

%% Exclude additional peptides
% In addition to filtering "bad" peptides, user can specify more peptides
% to remove from the analysis. Peptides numbers should be specified by
% their number in the peptide pool or the specified .csv files.

excld = zeros(maxSeq,1);

% list of peptides to exclude. For example, to exclude peptide 3, 7 and 89
% user will have to change the code below to excld([3 7 89]) = 1;
excld([]) = 1;

tdeltanormPercent(excld==1,:) = []; % FIRST remove excluded by INDEX
tdeltanormPercent(any(isnan(tdeltanormPercent),2),:) = []; % remove peptides with difference in deuteration equal to nan

%% Create the figure to show difference in deuteration levels for each peptide

fig = figure;
fig.Units = 'inches';
fig.Position = [2 2 8 6];

ax1 = axes('units', 'inches', 'position', [ 0.5, 2.5,  7,  3]);
ax2 = axes('units', 'inches', 'position', [ 0.5, 0.5,  5.5,  1.5]);
ax3 = axes('units', 'inches', 'position', [ 6.5, 0.5,  1,  1.5]);

subplot(ax1);
yticks([])
xlim([-offset AAcutoff-offset])
box on
hold on

for k = 1:length(tdeltanormPercent) % each sequence
    l = tdeltanormPercent(k,1); % position within peptide pool (u)
    x1 = dataset_1{1}(l) - offset; % AAstart
    x2 = dataset_1{2}(l) - offset; % AAend
    
    r = -(2*k + 2); % rectangle position
    val = tdeltanormPercent(k,2);
    if ~isnan(val)
        eval(['color = ' colorFunc '(val);']);
        rectangle('Position',[x1,r,x2-x1+1,1],'FaceColor',color,'LineStyle','none')
    end
end

txt_string = [fname_1(1:end-4) ' vs. ' fname_2(1:end-4)];
title(txt_string,'Interpreter', 'none')
hold off

% Plot the histogram
subplot(ax2);
box on
histogram(tdeltanormPercent(:,2),20)
title('Distribution of difference rates','Interpreter', 'none')
set(gcf,'Units','inches')
%set(gcf,'Position', [.5 .5 6 4])

% Plot the color legend for difference plot
subplot(ax3);
title('Legend','Interpreter', 'none')
ylim([-100 100])
xticks([])
yticks(-100:20:100)
colors_ranges = -100:1:100;

for k = 1:length(colors_ranges) % each sequence
    
    y1 = colors_ranges(k); % y position of left lower corner
    val = colors_ranges(k);
    if ~isnan(val)
        eval(['color = ' colorFunc '(val);']);
        rectangle('Position',[0,y1,1,1],'FaceColor',color,'LineStyle','none')
    end
end

% Save figure in pdf and fig formats
saveas(gcf , fig_name1, 'fig')
saveas(gcf , fig_name1, 'pdf')





















