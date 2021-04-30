%DualCh_allign reads two images and generates a transform matrix
%for aligning the two images. Originally from Pontus Nordenfelt 07/2013
%Adapted by Thomas S. van Zanten 04/2021

%% allign two channels with the possibility to check and re-adjust
%INPUT: 1)ch1: image of channel 1
%       2)ch2: image of channel 2 which will be alligned on top of channel 1
%       3)method: give the method for the two channel allignment:
%           'translation'	(x,y) translation.
%           'rigid'         Rigid transformation consisting of translation 
%                           and rotation.
%           'similarity'	Nonreflective similarity transformation consisting 
%                           of translation, rotation, and scale.
%           'affine'        Affine transformation consisting of translation, 
%                           rotation, scale, and shear.
%       4) 0 for no visualisation and adjustment, 1 for interactive adjustment
%       5)tFORM: if provided will directly be used for alignment
%OUTPUT: the corrected channel 2 & tFORM used
function [ch2_aligned, tFORM] = DualCh_align (ch1, ch2, method, vis, tFORM)

OK=0;%setting up the condition for adjustments

%set up image registration parameters
[optimizer, metric] = imregconfig('multimodal');
optimizer.MaximumIterations = 1000;
optimizer.InitialRadius = 0.0002;

%%%create transform if required
if nargin<5
    tFORM = imregtform(ch2, ch1, method, optimizer, metric);
end
%%%align ch2 image onto ch1 for visualization
ch2_aligned = imwarp(ch2,tFORM,'OutputView',imref2d(size(ch1)));

    if vis==1
        while OK==0
%%%CONVERT TO RGB FOR VISUALISATION OF THE POTENTIAL ALLIGNMENT SUCCESS%%%%
before(:,:,1)=imadjust(uint16(ch2));
before(:,:,2)=imadjust(uint16(ch1));
before(:,:,3)=zeros(size(ch1));
after(:,:,1)=imadjust(uint16(ch2_aligned));
after(:,:,2)=imadjust(uint16(ch1));
after(:,:,3)=zeros(size(ch1));

%%%SHOW DIFFERENCE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(12), set(gcf, 'Position',  [400, 500, 950, 450])
subplot(1,2,1),imshow(before);
subplot(1,2,2),imshow(after);

            button = questdlg('Is the transform successful?','Yes');
            if strcmp(button, 'Yes')
                OK=1;%done
                close Figure 12
            else
                optimizer.InitialRadius = optimizer.InitialRadius + 0.0001;%increase the transform range
                clear before after
                close Figure 12
            end    
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end