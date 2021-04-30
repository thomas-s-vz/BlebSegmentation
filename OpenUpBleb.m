%% open up the Bleb circle to a Bleb line
%Part of the code below was adapted from: https://in.mathworks.com/matlabcentral/
%answers/299530-how-to-convert-2d-image-of-circle-shape-to-2d-straight-line


%INPUT: 1)im: the image containing a circle
%       2)av: the rolling average for smoothening the circle fitting
%       3)width: width above and below the circle once opened up
%       4)centers_out: the x,y position of the circle in the image
%       5)radii_out: the radius of the circle
%4 and 5 can be determined through the Circle Hough Transform used in imfindcircles
%       6)ind: the index of the pixels on the circle periphery. If this
%input is empty the ind will be determined and figures will be made. When 
%ind is given the function will use them for giving the imRcorr output
%
%OUTPUT:1)imRcorr: the rectangular image of the circle periphery having a
%length of the perimeter and a width set by the user
%       2)ind: the index of all the pixels on the circle periphery
function [imRcorr, ind, R] = OpenUpBlebCP(im, av, width, centers_out, radii_out, ind)

fig=0;%with ind as an input argument the function will return the data at ind

if nargin==5%with ind ABSENT as an input argument the function will perform the fitting
    ind=[];
    fig=1;
end

L=floor(2*pi*radii_out);%expected length of the morphed graph
cx=centers_out(1);%coordinates of circle centres
cy=centers_out(2);

    if fig==1
        figure(1)
        subplot(2,3,1)
        hold on,
        plot(cx,cy,'yo');%indicate centre circle in figure (yellow circle)
        plot([cx 1],[cy 1],'r-.');%indicate the radial line over which the data is integrated (red line-dot)
    end

t=[0:360/L:360-360/L];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                  
R=radii_out+width/2;%take extra pixels for the outer measurement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=R*sind(t)+cx; y=R*cosd(t)+cy;% build outer perimeter
hL1=plot(x,y,'m');%automatically determined circle perimeter + the extra pixels (magenta)
hold all;
plot(hL1.XData,hL1.YData,'ro');%indicate the points along the line from which the profile is taken (red circles)

x_ref=hL1.XData;y_ref=hL1.YData;
% Sx=zeros(ceil(R),1);Sy=zeros(ceil(R),1);

Sx={};Sy={};
for k=1:1:numel(hL1.XData)
    Lx=floor(linspace(x_ref(k),cx,ceil(R)));
    Ly=floor(linspace(y_ref(k),cy,ceil(R)));
%       plot(Lx,Ly,'go');    % check
%       plot([cx x(k)],[cy y(k)],'r');
%       L1=unique([Lx;Ly]','rows');
    Sx=[Sx Lx'];Sy=[Sy Ly'];
end

sx=cell2mat(Sx);sy=cell2mat(Sy);
[s1 s2]=size(sx);

imRect=zeros(s1,s2);

for n=1:1:s2
    for k=1:1:s1
        imRect(k,n)=im(sy(k,n),sx(k,n));%%%be aware that x and y are inverted here 
    end
end

    if fig==1
        figure(1)
        subplot(2,3,2)
        imagesc(imRect)
    end
%% allign the Bleb membrane plane and substract local background
imRcorr=zeros(width+1,size(imRect,2));
tempIM=[imRect(:,end-(av/2)+1:end) imRect imRect(:,1:av/2)];%since its a circle the data can be extended

for k=1:size(imRect,2)
    peakF=mean(tempIM(:,k:k+av),2);%smoothen the data to make a clearer maximum
    
    if fig==1
        figure(1), subplot(2,3,3), plot(peakF,'b-')
        hold on, plot(mean(tempIM(:,k:k+av),2),'b-')
        %pause
        ind(k)=min(find(peakF==max(peakF(1+width/2:end-max([width/2 round((length(peakF)-10)/3)])))));%peak determination should remain within the assigned width
        figure(1), subplot(2,3,2), hold on, plot(k,ind(k),'wo')%indicate determined peak position
    end
    
    imRcorr(:,k)=imRect(ind(k)-width/2:ind(k)+width/2,k);%allow some space around the membrane    
end

    if fig==1
        figure(1)
        subplot(2,3,4:6)
        imagesc(imRcorr)
        %pause
    end
    
end