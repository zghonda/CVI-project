clear all
img = imread('imgsrc/Moedas2.jpg'); 
gray_img = rgb2gray(img);
% figure
% imshow(gray_img);

L= imsegkmeans(img,2);
bw1 = L==1;


%remove very small objects
bw1 = bwareaopen(bw1,20);
bw = ~bw1;

% flip background if needed
[B,L,N,A] = bwboundaries(bw);
if N==1
    bw = ~bw;
end

D = -bwdist(~bw);
Ld = watershed(D);
bw2 = bw;
bw2(Ld ==0 ) =0;
mask = imextendedmin(D,4);
D2 = imimposemin(D,mask);
Ld2 = watershed(D2);
BW = bw;
BW(Ld2 == 0) = 0;


%remove very small objects
BW = bwareaopen(BW,4);

%Remove border components
BW = imclearborder(BW);

%Detect boundaries
[B,L,N,A] = bwboundaries(BW);

% %plot boundaries
% figure
% imshow(BW)
% hold on
% for k=1:length(B)
%    boundary = B{k};
%    plot(boundary(:,2), boundary(:,1), 'r','LineWidth',1);
% end
% hold off

%statistics on regions
stats = regionprops('table',BW,'Area','Centroid','Perimeter','BoundingBox');
BoundingBox = stats.BoundingBox;
Centroid = stats.Centroid;
Area = stats.Area;
Perimeter = stats.Perimeter;
Circularity = min(4*pi* Area./(Perimeter.^2),1);

%Table containing stats on each object (sorted by Area) 
rowNames = strings(N,1);
for i=1:N
    rowNames(i) = sprintf('Object%d',i);
end

tblsorted = sortrows(table(Area,Perimeter,Circularity,Centroid,BoundingBox));
tblsorted.Row = rowNames;

tblsorted


%compute relative distances between images(centroids)
relDist = zeros(N);
c=Centroid;
for i=1 : size(Centroid)
    for j = 1 : size(Centroid)
        relDist(i,j) = sqrt((c(2*j)-c(2*i))^2+(c(2*j-1)-c(2*i-1))^2);
    end
end



% Derivative of Object Boundary

[Gmag, Gdir] = imgradient(BW,'prewitt');
figure
imshowpair(Gmag, Gdir, 'montage');
title('Gradient Magnitude, Gmag (left), and Gradient Direction, Gdir (right), using Prewitt method')


%Plot Centroids, Bounding boxes, Areas and Perimeters
figure 
imshow(img)
hold on 
plot(Centroid(:,1),Centroid(:,2),'b*')

 for k = 1 : size(tblsorted)
  thisBB = tblsorted.BoundingBox(k,:);
  rectangle('Position', [thisBB(1),thisBB(2),thisBB(3),thisBB(4)],...
  'EdgeColor','g','LineWidth',2 );
    str = sprintf('%s\nArea: %0.2f\nPerimeter: %0.2f',tblsorted.Row{k},...
        tblsorted.Area(k),tblsorted.Perimeter(k));
  text(thisBB(1)+3,thisBB(2)+thisBB(4)-35,str,'FontSize',6,'Color','g');
  
 end
 
 
tblrelDist = array2table(relDist,'RowNames',rowNames,'VariableNames',rowNames);

tblrelDist
 
filename = 'ImageObjectStatistics.xlsx';
writetable(tblsorted(:,1:N-1),filename,'Sheet',1,'Range','B2','WriteRowNames',true);
writetable(tblrelDist,filename,'Sheet',2,'Range','B2');

 
 %Object selection and image order relative to Area 
 [c,r,~] = impixel;
 hold off
 mask = bwselect(BW,c,r,8);
 mask3 = cat(3, mask, mask, mask);
 J  = img;
 J(~mask3) = 0;
 I=J;

 stats = regionprops(mask,'Area');
 tblsorted{:,1} = abs(tblsorted{:,1}-stats.Area);
 
 toDelete = tblsorted.Area == 0;
 tblsorted(toDelete,:) = [];
 
 while(size(tblsorted)>0)
  [~,idx] = min(tblsorted.Area);
  centroid = tblsorted.Centroid(idx,:);
  tblsorted(idx,:) = [];
  mask = bwselect(BW,centroid(1),centroid(2),8);
  mask3 = cat(3, mask, mask, mask);
 J  = img;
 J(~mask3) = 0;
 I = [I J];   
    
 end
 
 figure
 imshow(I)
