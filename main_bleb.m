%main_bleb has been written by Thomas S van Zanten last update 15th April 2021
%% Section 0: Initialisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
clc
addpath('path where all the general functions are located')
addpath('path where all the specific functions are located')
pathname=uigetdir; pathname=[pathname '/']; cd(pathname)
%% Section 1: Variable and image loading and conversion %%%%%%%%%%%%%%%%%%%
ws=80;%width of the ROI. Is either 60 or 80 during bleb selection process
minR=5;%minimum radius
maxR=30;%maximum radius: depends largely on the ROI size
blur=3;%gaussian blur during for circle determination 
av=6;%size rolling average along the membrane for re-allignment 
%--> should be an EVEN number
width=10;%width around the bleb membrane to allow deformations to be fitted 
%--> should be an EVEN number
avM=6;%thickness of the membrane + PFS (to include in the membrane summation) 
%--> should be an EVEN number and smaller than width
Bleb.area=[];%initializing the Bleb array to prevent errors in ROIselect.m

temp_im1=tiffread2b([pathname 'Ldn_blue.tif'],1,40); 
for i=1:length(temp_im1), im1(i,:,:)=temp_im1(i).data; end
temp_im2=tiffread2b([pathname 'Ldn_red.tif'],1,40);
for i=1:length(temp_im2), im2(i,:,:)=temp_im2(i).data; end
%temp_im3=tiffread2b([pathname 'Mov640_xlink.tif'],1,40);
%for i=1:length(temp_im3), im3(i,:,:)=temp_im3(i).data; end

cam.mode=0; cam.serialNo='0000'; cam.tInt='200ms';%camera settings for more info
cam.ROI=[0 0 size(im1,2) size(im1,3)];%see PhotoConvertIm.m

[im1, gain_ch1, offset_ch1] = PhotoConvertIm (im1, cam);
[im2, gain_ch2, offset_ch2] = PhotoConvertIm (im2, cam);
%[im3] = PhotoConvertIm (im3, cam, gain_ch2, offset_ch2);
%assuming ch3 is taken with camera of ch2

im1BG=tiffread2([pathname 'BG/BG_blue.tif']);im1BG=double(im1BG.data);
im2BG=tiffread2([pathname 'BG/BG_red.tif']);im2BG=double(im2BG.data);
%im3BG=tiffread2([pathname 'BG/BG_640.tif');im3BG=double(im3BG.data);
%im3BG = PhotoConvertIm (im3BG, cam, gain_ch2, offset_ch2);
im1BG = PhotoConvertIm (im1BG, cam, gain_ch1, offset_ch1); 
im2BG = PhotoConvertIm (im2BG, cam, gain_ch2, offset_ch2);
BG(1,:,:) = im1BG; BG(2,:,:) = im2BG;
im1GF=tiffread2([pathname 'BG/GF_blue.tif']);im1GF=double(im1GF.data);
im2GF=tiffread2([pathname 'BG/GF_red.tif']);im2GF=double(im2GF.data);
im1GF = PhotoConvertIm (im1GF, cam, gain_ch1, offset_ch1);
im2GF = PhotoConvertIm (im2GF, cam, gain_ch2, offset_ch2);
GF(1,:,:) = im1GF; GF(2,:,:) = im2GF;

	im1=ImageCorrection(im1,im1BG);
    [im2 GFactor]=ImageCorrection(im2,BG,'Ldn',GF);
    %im3=ImageCorrection(im3,im3BG);

im1(find(im1<0))=0;im1(find(im1==Inf))=0;im1(isnan(im1))=0;
im2(find(im2<0))=0;im2(find(im2==Inf))=0;im2(isnan(im2))=0;
%im3(find(im3<0))=0;im3(find(im3==Inf))=0;im3(isnan(im3))=0;
clear temp_im1 temp_im2 temp_im3 im1BG im2BG BG im3BG im1GF im2GF GF
cd(pathname), save('bleb_fullimages.mat')
%% Section 2: Indicating the Blebs to be analysed per image %%%%%%%%%%%%%%%
for j=1:size(im1,1)
    
   	ch1(:,:)=im1(j,:,:);
    ch2(:,:)=im2(j,:,:);
    %ch3(:,:)=im3(j,:,:);
 	
    No=size(Bleb,2);

    figure(9)
    set(gcf, 'Position',  [400 50 950 900])
    imagesc(ch1+ch2);
    %imagesc(ch1+ch2+ch3);
        title(['Image ' num2str(j) ' of ' num2str(size(im1,1)) ' images: PLEASE DRAW BACKGROUND'])
 
roi = imfreehand%%%Draw a single region that will become the image-associated background

BW=createMask(roi);
BG(1)=sum(sum(BW.*ch1))/sum(sum(BW));
BG(2)=sum(sum(BW.*ch2))/sum(sum(BW));
%BG(3)=sum(sum(BW.*ch3))/sum(sum(BW));
delete(roi);

    fin=0;
    while fin==0
       [Bleb, fin] = ROIselect(ch1+ch2, ws, 'object', Bleb);
    end
        
    for i=No+1:size(Bleb,2)
        Bleb(i).BG=BG;
        Bleb(i).frame=j;
        Bleb(i).im1=ch1(Bleb(i).area(3):Bleb(i).area(4),Bleb(i).area(1):Bleb(i).area(2));
        Bleb(i).im2=ch2(Bleb(i).area(3):Bleb(i).area(4),Bleb(i).area(1):Bleb(i).area(2));
        %Bleb(i).im3=ch3(Bleb(i).area(3):Bleb(i).area(4),Bleb(i).area(1):Bleb(i).area(2));
    end
        
    clear No BG BW ch1 ch2 ch3 roi fin
    close Figure 9
    
end
Bleb(1)=[]; cd(pathname), save('bleb_ROIs.mat', 'Bleb')
%% Section 3: Centering, aligning and opening up the Bleb for analysis %%%%
for i=1:size(Bleb,2)
%% Alligning different channels in Bleb ROI 
    ch1(:,:)=Bleb(i).im1; ch1(find(ch1<0))=0;
    ch2(:,:)=Bleb(i).im2; ch2(find(ch2<0))=0;
    %ch3(:,:)=Bleb(i).im3; ch3(find(ch3<0))=0;
ws=size(ch2,1); 

        [ch2] = DualCh_align (ch1, ch2, 'translation', 1);
        %[ch3] = DualCh_align (ch1, ch3, 'translation');
        im(:,:)=ch1+ch2;%+ch3;
%for freestanding GPMVs imfindcircles will work, but for cell-attached blebs manual intervention can be faster
%% FREESTANDING / AUTOMATED Bleb find
imblur = imgaussfilt(im,blur);
[centers_out, radii_out, metric_out] = imfindcircles(imblur,[minR maxR],...
        'ObjectPolarity','bright','Sensitivity',0.85,'EdgeThreshold', 0.1);
%'Sensitivity' and/or 'EdgeThreshold' both need a value between 0:1 and proper 
%identification is dependent on raw data
tempmin=mean((ws/2-centers_out).^2,2);%in the possible case that multiple 
%circles are found the one closest to the slected center position is chosen
centers_out=centers_out(find(tempmin==min(tempmin)),:);
radii_out=radii_out(find(tempmin==min(tempmin)),:);
metric_out=metric_out(find(tempmin==min(tempmin)),:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(11)
    set(gcf, 'Position',  [730, 55, 950, 900])
    imagesc(im)
    title(['Bleb ' num2str(i) ' of ' num2str(size(Bleb,2))])
    
    hold on
    h=viscircles(centers_out, radii_out,'EdgeColor','r');    
%% Possibility to change the calculated values
%       A{1,1} = num2str(1);
%                 A{1,1} = num2str(centers_out(1));
%                 A{2,1} = num2str(centers_out(2));
%                 A{3,1} = num2str(radii_out);
% 
%                 A{1,2} = 'x-pos';
%                 A{2,2} = 'y-pos';
%                 A{3,2} = 'radius';
% prompt = A(1:3,2); def = A(1:3,1); AddOpts.Resize = 'on'; 
% AddOpts.WindowStyle = 'normal'; AddOpts.Interpreter = 'tex';
% Parameters(1:3) = inputdlg(prompt,'Please change if required:', 1, def, AddOpts); 
% centers_out=[str2num(Parameters{1}), str2num(Parameters{2})]; radii_out=str2num(Parameters{3});
%% CELL ATTACHED / MANUAL Bleb find
% r=imrect;
% Parameters=r.getPosition; delete(r)
% centers_out=[Parameters(1), Parameters(2)]; radii_out=Parameters(3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close Figure 11
    figure(1)
    set(gcf, 'Position',  [730, 55, 950, 900])
    subplot(2,3,1)
    imagesc(im)
    title(['Bleb ' num2str(i) ' of ' num2str(size(Bleb,2))])
    
    hold on
    h=viscircles(centers_out, radii_out,'EdgeColor','b');
%% Open up Bleb circle to a Bleb rectangular image
    [Bleb(i).ini, Bleb(i).ind] = OpenUpBleb(im, av, width, centers_out, radii_out);
    [Bleb(i).ch1, ind] = OpenUpBleb(ch1, av, width, centers_out, radii_out, Bleb(i).ind);
    [Bleb(i).ch2, ind] = OpenUpBleb(ch2, av, width, centers_out, radii_out, Bleb(i).ind);
    %[Bleb(i).ch3, ind] = OpenUpBleb(ch3, av, width, centers_out, radii_out, Bleb(i).ind);

Bleb(i).centers=centers_out;
Bleb(i).radii=radii_out;
Bleb(i).AdjustedRadii=mean(round(radii_out+width/2)-Bleb(i).ind);
Bleb(i).circularity=1-(std(round(radii_out+width/2)-Bleb(i).ind)/Bleb(i).AdjustedRadii);
Bleb(i).metric=metric_out; 

clear r A def Parameters prompt h AddOpts imblur im ch1 ch2 ch3...
    ind centers_out radii_out metric_out tempmin
close Figure 1
end
Bleb = rmfield(Bleb, 'ini');
%Bleb([11, 18])=[];%This will delete the indicated blebs from the file
cd(pathname), save('blebs_identified.mat','Bleb')
%% Section 4: Analyzing all the Blebs and extracting relevant data %%%%%%%%
j=1;
for i=1:length(Bleb)  
    if length(Bleb(i).ch1)>0
%trace extraction        
    Bleb(i).ch1Trace=sum(Bleb(i).ch1((width/2)+1-(avM/2):...
        (width/2)+1+(avM/2),:),1)-(avM+1)*Bleb(i).BG(1);
    Bleb(i).ch2Trace=sum(Bleb(i).ch2((width/2)+1-(avM/2):...
        (width/2)+1+(avM/2),:),1)-(avM+1)*Bleb(i).BG(2);
    %Bleb(i).ch3Trace=sum(Bleb(i).ch3((width/2)+1-(avM/2):...
        %(width/2)+1+(avM/2),:),1)-(avM+1)*Bleb(i).BG(3);
    Bleb(i).totTrace=Bleb(i).ch1Trace+Bleb(i).ch2Trace;

%An value for anisotropy calculation
    %Bleb(i).anTrace=(Bleb(i).ch1Trace-Bleb(i).ch2Trace)./Bleb(i).totTrace;
%GP value for Laurdan calculation 
    Bleb(i).gpTrace=(Bleb(i).ch1Trace-Bleb(i).ch2Trace)./Bleb(i).totTrace;   
    
%trace correlation        
    %[Bleb(i).IntR, Bleb(i).IntP] = corrcoef(Bleb(i).totTrace,Bleb(i).ch3Trace);
    %[Bleb(i).anR, Bleb(i).anP] = corrcoef(Bleb(i).anTrace, Bleb(i).ch3Trace);
    %[Bleb(i).gpR, Bleb(i).gpP] = corrcoef(Bleb(i).gpTrace, Bleb(i).ch3Trace);
    
%whole bleb values     
    Bleb(i).meanInt=mean(Bleb(i).totTrace); Bleb(i).stdInt=std(Bleb(i).totTrace);
    Bleb(i).totInt=sum(Bleb(i).totTrace);
    %Bleb(i).AN=mean(Bleb(i).anTrace); Bleb(i).ANstd=std(Bleb(i).anTrace);
    Bleb(i).GP=mean(Bleb(i).gpTrace); Bleb(i).GPstd=std(Bleb(i).gpTrace);

%rotational dependence of the measured values  
    CircTOT(:,j) = interp1([1:length(Bleb(i).totTrace)],Bleb(i).totTrace,...
        linspace(1, length(Bleb(i).totTrace), 360),'spline');
    %CircCH3(:,j) = interp1([1:length(Bleb(i).ch3Trace)],Bleb(i).ch3Trace,...
        %linspace(1, length(Bleb(i).ch3Trace), 360),'spline');
    %CircAN(:,j) = interp1([1:length(Bleb(i).anTrace)],Bleb(i).anTrace,...
        %linspace(1, length(Bleb(i).anTrace), 360),'spline');
    CircGP(:,j) = interp1([1:length(Bleb(i).gpTrace)],Bleb(i).gpTrace,...
        linspace(1, length(Bleb(i).gpTrace), 360),'spline');
    j=j+1; 
    
    clear CC CCg
    end    
end
CircTOT=CircTOT./nanmean(CircTOT,1); %CircCH3=CircCH3./nanmean(CircCH3,1);
%CircAN=CircAN./nanmean(CircAN,1); 
CircGP=CircGP./nanmean(CircGP,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%DATA ON BACKGROUND%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for i=1:length(bleb)
%   BG(i,:)=Bleb(i).BG;
%end
%BGratio=unique(BG(:,1)./BG(:,2));
%BGtotal=unique(BG(:,1)+BG(:,2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
cd(pathname), save('blebs_analysed_final.mat')
