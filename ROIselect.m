%% Allow multiple ROIs to be selected in an image
%INPUT:     1)im: the image onto which to select multiple ROIs
%           2)ws: window size of the ROI
%           3)identity: 'object' will indicate a white '+'-mark at the
%           object of interest and 'square' will indicate a red squared ROI
%           3)ROI: give the pre-existing ROIs(=[] if it is the start)
%OUTPUT:    1)ROI: the updated ROIs
%           2)fin: using this in an outside while loop allows the user to
%only exit multiple ROI selection when fin==1
function [ROI, fin] = ROIselect(im, ws, identity, ROI)
ws=ws/2;
%%%%%%%END PEAK SELECTION ONLY WHEN OK BUTTON IS PUSHED%%%%%%%%%%%%%%%%%%
button = questdlg('Do you want to select another single ROI?','selecting','Yes','No','Yes');
        if strcmp(button,'No') 
            fin=1;
        elseif strcmp(button,'Yes')
            fin=0; i=length(ROI)+1;
            figure(9)
                title(['Currently you have ' num2str(i-1) ' ROIs : please select another']);
            [x,y] = ginput(1);
            x_c = round(x);
            y_c = round(y);
%%%%%%%CORRECTION FOR POSSIBLE EDGE EFFECTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
                if x_c-ws<=1
                    xl1=1; xshift=abs(x_c-ws)+1;
                    xl2=x_c+ws+xshift;
                elseif x_c+ws>=size(im,1)
                    xl2=size(im,1); xshift=size(im,1)-abs(x_c+ws);  
                    xl1=x_c-ws+xshift;
                else
                    xl1=x_c-ws;
                    xl2=x_c+ws;
                end

                if y_c-ws<=1
                    yl1=1; yshift=abs(y_c-ws)+1;
                    yl2=y_c+ws+yshift;
                elseif y_c+ws>=size(im,1)
                    yl2=size(im,1); yshift=size(im,1)-abs(y_c+ws);  
                    yl1=y_c-ws+yshift;
                else
                    yl1=y_c-ws;
                    yl2=y_c+ws;
                end    
%%%%%%%EVERY ROI CENTER IS MARKED IN THE ORIGINAL LARGE IMAGE%%%%%%%%%%%%%%
            figure(9);
            hold on
            if identity=='object'
                plot(x,y,'w+')
            elseif identity=='square'
                rectangle('Position', [xl1 yl1 2*ws 2*ws], 'Edgecolor', 'r')
            end
            ROI(i).area=[xl1 xl2 yl1 yl2];
        end     
end  