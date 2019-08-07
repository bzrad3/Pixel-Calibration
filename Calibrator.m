
clear;home;
data=dlmread('/Volumes/BRad/Microscope Images/TIRFM data/2019/20190124_calibration/Prime95B/Profile Plots/100X.csv',',',2);
yy=data(:,2); %Read pixel value of thresholded image
xx=data(:,1); %Read pixel number
pixel=[];
hh=yy==255;%Create vector from pixel Value that is either 1 or 0. Removes weird values.
mm=diff(hh) >= 1; %Look for front edge of each marker
for ii=1:length(mm)
    if mm(ii) == 1
        pixel=[pixel xx(ii)];
    end
end
pixel=transpose(pixel);

distance=transpose([1:length(pixel)].*10); %each division is 10 micrometers
line=polyfit(pixel,distance,1);
fitline=polyval(line,pixel);

yresid = distance-fitline;
SSresid = sum(yresid.^2);
SStotal = (length(distance)-1 * var(distance));
rsq = 1- SSresid/SStotal;

n=length(pixel);
SSxx=sum(pixel.^2)-(n*mean(pixel)^2);
SSyy=sum(distance.^2)-(n*mean(distance)^2);
SSxy=sum(pixel.*distance)-(n.*mean(pixel).*mean(distance));

s=sqrt((SSyy-(SSxy^2/SSxx))/(n-1));
errslope=s/sqrt(SSxy);

slope=line(1);


'r-squared' 
rsq
'Conversion (pixel/micrometer)'
slope 
errslope
plot(pixel,distance,'o',pixel,fitline,'r')
%Title Graph
title('Plot Values 100X Prime95B')
% Create ylabel
ylabel('Distance (nm)');

% Create xlabel
xlabel('Pixels');

saveas(gcf,'/Volumes/BRad/Microscope Images/TIRFM data/2019/20190124_calibration/Prime95B/Profile Plots/Plot Values 100X Prime95B.png')

