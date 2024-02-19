%
% NAME:
%               Particle_measurement_AICL_v2
%
% PURPOSE:
%               For batch analysis of particles in individual images.
%
%               Individual images pass a rolling-ball & bpass filter for removing background and noises. Regions of interest (ROIs) are then
%               identified after thresholding.
%               The number\area\intensity\signal-to-background ratio
%               (SBR)\length of individual particles in a single image are calculated.
%
%               Input:
%               Single .tif images (no videos) in a folder and a subfolder
%
%               Output:
%               Analysed images and the ones with numbering in .png file in a ROI subfolder; a Count and Result .txt files.
%               Area is calculated in pixel length squared within a ROI.
%               Intensity is the sum of values of individual pixels.
%               SBR is the corrected intensity, revealing as the sum of the difference values between individual pixels and average value of
%               background, divided by the backround.
%               Length is a rough approach, calculated by thinning individual ROIs to an one-dimensional skeleton.
%               The more the value close to the diffraction limit, the more discrete and distorted the value is.
%
%
%               Require the script of Length_calculation_v2.m, bpass.m, and Matlab R2015b
%
%
%               Written by Jason C Sang, University of Cambridge,
%               June 2016
%
%               Updated on 2017\04\08
%                       Optimised the calculating speed by ten-times by rewriting looping structures and the length calculation method.
%                       Approximately less than 10 sec for analysing one image.
%                       User-friendly interface. Now all parameters can be fine-tuned in the Parameter setting section.
%
%               Updated on 2017\04\08
%                       Bug fixed.
%                       Applied Gaussian fit for ROI finding before bw thresholding.
%
%               Updated on 2017\04\08
%                       Now pixel size can be altered in the setting. The previous setting is approximated to a pixel size of 225 nm.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
clear all

% Append the SPT codes
% addpath('C:\Users\mjmor\Documents\MATLAB\')
% addpath('C:\Users\mjmor\Documents\MATLAB\Aggregate size\main_tracking source code\')


%MainPath = 'C:\Users\mjmor\Documents\Imperial data\Cameron\Test for Hailey\Data'; % Select the folder for averaged images
MainPath = '/rds/general/user/hg1222/projects/imaging/live/Reconstructed Images Brain soak/DLB';
FolderList = dir(MainPath); %% Require a subfolder (eg. averaged) under the initial folder
FolderList = FolderList(3:end);
FolderList;
fprintf('Start %s\n', datestr(now,'HH:MM'))

%%
% Change folder
for folderN = 1:length(FolderList)
    
    pathname = [MainPath '/' FolderList(folderN).name '/'];
    pathname1 = [MainPath '/' FolderList(folderN).name];
    NameList = dir(pathname);
    NameList = NameList(3:end);
    dname = {NameList(:).name};
    
    mydir  = MainPath;
    idcs   = strfind(mydir,filesep);
    newdir = mydir(1:idcs(end)-1);
    
    %% Parameter setting
    
    
    localised_tracking = 0;                  % Set 1 for tracking restricted area, specified by the centre (X,Y) and the area size; 0 for entire images.
    
    X = 198;                                  % Coordinate of X
    Y = 405;                                   % Coordinate of Y
    stepper = 40;                          % Size of ROI in pixels
    
    mag = 2;                                        % Magnification of images during calculation was 10
    
    pixel_size = 12;                            % pixel size for length calculation (nm)
    
    RB_size = 10;                                 % Rolling-ball size in pixels, was 5.
    
    lnoise = 1;                                        % was 1 Level of noise size in pixels; usually 3. May be set to 0 or false, in which case only the highpass "background subtraction" operation is performed.
    
    lobject = 100;                                     % Level of object size in pixels somewhat larger than a typical object; usually 30-50, Syn=50, PrP=30 750
    
    th = 0.00000001;                                   % Threshold for b\w conversion; usually 0.005-0.015. Highly dependent to the original image intensity. Was=0.00001
    
    save_status = 1;                              % Set 1 for saving figures\data; 0 for not saving.
    
    gaus = 0.01;                              % Gaus filter setting, was 3 (2)
    
    %% Parameter file
    Parameter(1)=mag;
    Parameter(2)=pixel_size;
    Parameter(3)=RB_size;
    Parameter(4)=lnoise;
    Parameter(5)=lobject;
    Parameter(6)=th;
    Parameter(7)=gaus;
    
   
    
    
    %%
    if localised_tracking==1 % Stepper parameters for localized tracking
        
        X = yc;
        Y = xc;
        
        x1 = xc-stepper;
        x2 = xc+stepper;
        y1 = yc-stepper;
        y2 = yc+stepper;
    end
    
    for imageN = 1:size(dname,2)


    figuredir = [MainPath '_Results/Figures/' FolderList(folderN).name '/ROI/ROI_' num2str(imageN) '/']; % Destination for the figure saving; don't have to change the magenta part anymore (MM 17\11\20)
    
    largefiguredir = [MainPath '_Results/LargeFigures/' FolderList(folderN).name '/ROI/ROI_' num2str(imageN) '/']; % Destination for the figures that have a length of greater than 100
    
    datadir = [MainPath '_Results/Data/' FolderList(folderN).name '/']; % Destination for the data saving; don't have to change the magenta part anymore (MM 17\11\20)
    
        %% Image input and processing
        
        display(['Processing ', num2str(imageN) ,'th image in the "',  FolderList(folderN).name, '" folder:   ', dname{imageN}])
        fprintf('Processing %s\n', datestr(now,'HH:MM'))
        % Check existence of folder
        if save_status==1
            if exist(figuredir, 'dir')==0
                warning('off','MATLAB:MKDIR:DirectoryExists')
                mkdir(figuredir);
                mkdir(datadir);
            end
        end
        
        save ([datadir 'Parameter.txt'], 'Parameter', '-ascii');
        
        fname = [pathname dname{imageN}];
        %fname = [pathname1];
        info = imfinfo(fname, 'tif');
        images = numel(info);
        
        % Input image
        clear imDataOri
        
        for i = 1:images
            % Bin Up
            imTemp = imread(fname, i);
            imTemp(imTemp == 0) = mean2(imTemp);
            
            if localised_tracking==1
                imDataOri(:,:) = imTemp(x1:x2,y1:y2); % if using stepper
            elseif localised_tracking==0
                imDataOri(:,:) = imTemp(:,:); % if NOT using stepper
            end
            clear imTemp
        end
        clear i
        
        list = [];
        boundList = [];
        iList = [];
        sbrList = [];
        roiList = [];
        %%
        for i = 1:images
            
            se = strel('disk', RB_size);
            imDataB = imtophat(imDataOri, se); % rolling-ball filtered image
            
            imDataG(:,:) = imDataOri - imDataB; % background per unit
            
            imDataBG = imresize(imDataG, mag,'nearest'); %10x background
            
            imDataF = imresize(imDataB, mag,'nearest'); % 10x rolling-ball filtered image
            
            imData = imresize(imDataOri, mag,'nearest'); % 10x original image
            
            
            imDataGauss = imgaussfilt(imData,gaus); % Gaussian filter for images
            
            imDataF(:,:)=bpass(imgaussfilt(imDataF,gaus),lnoise,lobject); % Filter for image thresholding; usually 30-50, Syn=50, PrP=30
            
            clear imDataB imDataG
            
            
            
            if max(max(imDataF))>0
                %% Identify and track spots
                
                
                
                BW = im2bw(imDataF, th); % Transform to BW for boundaries; set the threshold; usually 0.005-0.015, Syn=0.019, PrP=0.0055
                BW = imfill(BW,'holes');
                
                clear out intensity sbr image blank background
                
                
                
                image = imDataGauss;
                blank = imDataBG;
                
                clear imT
                
                % Measure features
                [imT roiN] = bwlabel(BW);
                boundaries = bwboundaries(BW, 'noholes');
                roi=[];
                roi=roiN;
                labelOri=imT.*double(imData);
                testStats=regionprops(labelOri,'all');
                
                F_area=([testStats.Area])';
                F_solidity=([testStats.Solidity])';
                F_eccentricity=([testStats.Eccentricity])';
                
                imT_bw = imbinarize(imT);
                imT_bw_skel = bwskel(imT_bw);
                imT_bw_skel_label = double(imT_bw_skel).*imT;
                skelStats = regionprops(imT_bw_skel_label,'all');
                
                F_length = ([testStats.MajorAxisLength])';
                
                Count = size(testStats,1);
                
                Result=zeros(Count,7);
                Result(:,1)=imageN;
                Result(:,2)=(1:Count)';
                Result(:,3)=F_area./mag;
                Result(:,4)=F_solidity;
                Result(:,5)=F_eccentricity;
                Result(:,6)=F_length;
                             
                if save_status==1
                    save ([datadir 'Count.txt'], 'Count', '-ascii', '-append');
                    save ([datadir 'Result' num2str(imageN) '.txt'], 'Result', '-ascii', '-append');
                end
                
                display(['Found ', num2str(roiN), ' particles.'])
                fprintf('Found %s\n', datestr(now,'HH:MM'))
                
                for agg=1:size(testStats,1)
                    

                    
                    xcorner=floor(testStats(agg).BoundingBox(1));
                    ycorner=floor(testStats(agg).BoundingBox(2));
                    width=floor(testStats(agg).BoundingBox(3));
                    height=floor(testStats(agg).BoundingBox(4));
                    
                    space=40;              %RB_size-1;
                    % TEST
                    if testStats(agg).Area<35
                        space=6;
                    end
                    
                    
                    close all
                    
                    aggIm=imData(ycorner-space:ycorner+height-1+space,xcorner-space:xcorner+width-1+space);
                    
                    hull=testStats(agg).Image;
                    hullSpace=zeros(size(aggIm));
                    hullSpace(space+2:space+height+1,space+2:space+width+1)=hull;
                    aggIm=imclearborder(aggIm,4);
                    
                    if max(aggIm(:))==0
                        aggIm=imData(ycorner-space:ycorner+height-1+space,xcorner-space:xcorner+width-1+space);
                    end
                    
                    aggIm=imresize(aggIm,1/mag);
                    aggIm=uint16(aggIm);
                    
                    % Check if length is greater than 100 and set directory accordingly 
                    if F_length(agg) > 100
                        savedir = largefiguredir;
                    else
                        savedir = figuredir;
                    end
                    
                    % Ensure the directory exists and create it if it doesn't 
                    if ~exist(savedir, 'dir')
                        mkdir(savedir);
                    end 
                    
                    % Save the image 
                    imwrite(aggIm, [savedir 'FOV_' num2str(i) '_aggreg_' num2str(agg) '.tif']);
                    
                    
                    
                end

                
                
            end
            clear roi out intensity sbr
            
            
            
            
            
            
            
            
            
            
            
            
            
        end
        
        
        
    end
    fprintf('End of original analysis %s\n', datestr(now,'HH:MM'))
    
    
 
    fprintf('Finished spool %s\n', datestr(now,'HH:MM'))
    close all
end


%clear boundListN roiListN iListN sbrListN lengthListN NList frame imDataBG imDataG imDataF imDataGauss BW roi boundaries image blank out intensity sbr results


% Concat_results_v3
% Heat_map

display('Execution completed!')
toc


