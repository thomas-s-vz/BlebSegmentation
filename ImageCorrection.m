%% Background correct and/or Gfactor correct an image
%INPUT  (1)im:      the image or image stack (frame,:,:) to be corrected
%       (2)imBG:    background image, should be a single image if only
%                   background correction needs to be performed. Needs to contain the
%                   2 channels if Gfactor correction is also performed (ch,:,:).
%       (3)type:    'Ldn' for laurdan Gfactor correction, 'An' for
%                   anisotropy Gfactor correction
%       (4)imGF:    the two channel Gfactor images (ch,:,:)
%       (5)tFORM:   Image transform to align all ch2's
%
%OUTPUT An image or image stack (frame,:,:) which is background and Gfactor
%       corrected & Gfactor image if used
%
%Written by Thomas S. van Zanten, last modified 15th April 2021
function [im_corr GFactor] = ImageCorrection (im, imBG, type, imGF, tFORM)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(size(im))==3
    frames=size(im,1);
else
    frames=1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wb = waitbar(0,'Image Correction');
if nargin==2
    imBG=medfilt2(imBG,[9 9]);
        for i=1:frames
            waitbar(i/frames,wb);
            if frames==1
                im_corr=double(im)-imBG;
            elseif frames>1
                temp_im(:,:)=double(im(i,:,:));
                im_corr(i,:,:)=temp_im-imBG;
            end   
        end
else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if nargin==4
                    tFORM=affine2d([1 0 0; 0 1 0; 0 0 1]);
                end
    BGch1(:,:)=imBG(1,:,:);BGch1(:,:)=medfilt2(BGch1,[9 9]);
    BGch2(:,:)=imBG(2,:,:); BGch2 = imwarp(BGch2,tFORM,'OutputView',imref2d(size(BGch1))); BGch2(:,:)=medfilt2(BGch2,[9 9]);

        if type=='Ldn'
            imGFch1(:,:)=imGF(1,:,:); imGFch1=imGFch1-BGch1; imGFch1=medfilt2(imGFch1,[9 9]);
            
            imGFch2(:,:)=imGF(2,:,:); 
            imGFch2 = imwarp(imGFch2,tFORM,'OutputView',imref2d(size(imGFch1)));
            imGFch2=imGFch2-BGch2; imGFch2=medfilt2(imGFch2,[9 9]);
            
            gp=(imGFch1-imGFch2)./(imGFch1+imGFch2);
            GFactor=(0.208+0.208*gp-gp-1)./(gp+0.208*gp-0.208-1);
            GFactor(find(abs(GFactor)==Inf))=0; GFactor(isnan(GFactor))=0;
        elseif type=='An'
            imGFch1(:,:)=imGF(1,:,:); imGFch1=imGFch1-BGch1; imGFch1=medfilt2(imGFch1,[9 9]);
            
            imGFch2(:,:)=imGF(2,:,:);
            imGFch2 = imwarp(imGFch2,tFORM,'OutputView',imref2d(size(imGFch1)));
            imGFch2=imGFch2-BGch2; imGFch2=medfilt2(imGFch2,[9 9]);
            
            GFactor=imGFch1./imGFch2;
            GFactor(find(abs(GFactor)==Inf))=0; GFactor(isnan(GFactor))=0;
        end
        
	for i=1:frames
        waitbar(i/frames,wb);
     	if frames==1
        	im_corr=(double(im)-BGch2).*GFactor;
      	elseif frames>1
          	temp_im(:,:)=double(im(i,:,:));
           	im_corr(i,:,:)=(temp_im-BGch2).*GFactor;
        end  
 	end
        
end
close(wb)
im_corr(find(im_corr<0))=0;
end