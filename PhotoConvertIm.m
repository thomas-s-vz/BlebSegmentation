%% CONVERT AN IMAGE FROM ADU TO PHOTONS
%INPUT  (1)image: the image or image stack (frame,:,:) to be converted from ADU to photo-electrons
%       (2)cam variable:
%           cam.mode:       EMCCD=0, Prime95B 12b Sensitivity mode=1, Prime95B 12b Balanced mode=2,
%                           Prime95B 12b Full well mode=3, Prime95B 16b HDR mode=4
%           cam.serialNo:   Serial number of the camera, if not know use '0000'
%           cam.ROI:        ROI of the camera that is used, default should be [0 0 size(im,1) size(im,2)]
%       (3)&(4) single value gain and offset for EMCCD images. If not given
%       the values will be automatically calculated using PCFO
%
%OUTPUT An image or image stack (frame,:,:) where the pixel values are now number of photons
%       & average gain & average offset values that were used
%
%Written by Thomas S. van Zanten, last modified 15th April 2021
function [image_photons, gain, offset] = PhotoConvertIm (im, cam, gain, offset)

addpath('/Users/thomas/Documents/MATLAB/added_codes/tiffread2/')
cameraPath='/Users/thomas/Documents/MATLAB/added_codes/sCMOS_Cal/';

mode=cam.mode; serialNo=cam.serialNo; ROI=cam.ROI;

if length(size(im))==3
    frames=size(im,1);
else
    frames=1;
end
%% EMCCD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<4 && mode==0
    addpath('/Users/thomas/Documents/MATLAB/added_codes/SingleShotGain/')
    addpath(genpath('/Applications/dip/'))%this initialises the DIPlib software within MatLab
    dip_initialise
        for i=1:frames
            if frames==1
                temp_im(:,:)=double(im);
            elseif frames>1
                temp_im(:,:)=double(im(i,:,:));
            end
            [gain(i), offset(i)] = pcfo(temp_im, 0.9, 0, 1, 0, [3 3]);
        end
gain=nanmean(gain); 
offset=nanmean(offset);
end

    wb = waitbar(0,'Converting the images to photo-e');
if mode==0
            for i=1:frames
                waitbar(i/frames,wb);
                if frames==1
                    temp_im=double(im);
                    image_photons=(double(temp_im)-offset)/gain;
                elseif frames>1
                    temp_im(:,:)=im(i,:,:);
                    image_photons(i,:,:)=(double(temp_im)-offset)/gain;
                end
            end
%% sCMOS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif mode==1
    td=tiffread2([cameraPath serialNo '_Sensitivity_dark.tif']);
    offset=td.data(ROI(2)+1:ROI(2)+ROI(4),ROI(1)+1:ROI(1)+ROI(3)); offset=double(offset);
    tg=tiffread2([cameraPath serialNo '_Sensitivity_gain.tif']);
    gain=tg.data(ROI(2)+1:ROI(2)+ROI(4),ROI(1)+1:ROI(1)+ROI(3)); gain=double(gain)/10000;
            for i=1:frames
                waitbar(i/frames,wb);
                if frames==1
                    image_photons=(double(im)-offset)./gain;
                elseif frames>1
                    temp_im(:,:)=double(im(i,:,:));
                    image_photons(i,:,:)=(temp_im-offset)./gain;
                end
            end    
gain=nanmean(nanmean(gain));
offset=nanmean(nanmean(offset));
elseif mode==2
    td=tiffread2([cameraPath serialNo '_Balanced_dark.tif']);
    offset=td.data(ROI(2)+1:ROI(2)+ROI(4),ROI(1)+1:ROI(1)+ROI(3)); offset=double(offset);
    tg=tiffread2([cameraPath serialNo '_Balanced_gain.tif']);
    gain=tg.data(ROI(2)+1:ROI(2)+ROI(4),ROI(1)+1:ROI(1)+ROI(3)); gain=double(gain)/10000;
            for i=1:frames
                waitbar(i/frames,wb);
                if frames==1
                    image_photons=(double(im)-offset)./gain;
                elseif frames>1
                    temp_im(:,:)=double(im(i,:,:));
                    image_photons(i,:,:)=(temp_im-offset)./gain;
                end
            end 
gain=nanmean(nanmean(gain));
offset=nanmean(nanmean(offset));
elseif mode==3
    td=tiffread2([cameraPath serialNo '_Fullwell_dark.tif']);
    offset=td.data(ROI(2)+1:ROI(2)+ROI(4),ROI(1)+1:ROI(1)+ROI(3)); offset=double(offset);
    tg=tiffread2([cameraPath serialNo '_Fullwell_gain.tif']);
    gain=tg.data(ROI(2)+1:ROI(2)+ROI(4),ROI(1)+1:ROI(1)+ROI(3)); gain=double(gain)/10000;
             for i=1:frames
                waitbar(i/frames,wb);
                if frames==1
                    image_photons=(double(im)-offset)./gain;
                elseif frames>1
                    temp_im(:,:)=double(im(i,:,:));
                    image_photons(i,:,:)=(temp_im-offset)./gain;
                end
            end 
gain=nanmean(nanmean(gain));
offset=nanmean(nanmean(offset));
elseif mode==4
    td=tiffread2([cameraPath serialNo '_HDR_dark.tif']);        
    offset=td.data(ROI(2)+1:ROI(2)+ROI(4),ROI(1)+1:ROI(1)+ROI(3)); offset=double(offset);
    tg=tiffread2([cameraPath serialNo '_HDR_gain.tif']);     
    gain=tg.data(ROI(2)+1:ROI(2)+ROI(4),ROI(1)+1:ROI(1)+ROI(3)); gain=double(gain)/10000; 
            for i=1:frames
                waitbar(i/frames,wb);
                if frames==1
                    image_photons=(double(im)-offset)./gain;
                elseif frames>1
                    temp_im(:,:)=double(im(i,:,:));
                    image_photons(i,:,:)=(temp_im-offset)./gain;
                end
            end 
gain=nanmean(nanmean(gain));
offset=nanmean(nanmean(offset));
end
close(wb)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%image_photons=round(image_photons);
image_photons(find(image_photons<0))=0;
end