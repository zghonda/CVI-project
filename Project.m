% We used 2 differents methods to identify the number of object and theirs
% properties. The first one is using the k-means algorithm and the second 
% one is using the edge detection algorithm. We choose to use the k-means
% method because when computing the results with the edge detection algorithm, 
% we have some difficulties to identify 2 touching objects and objects 
% with a hole in it. Therefore, the results with the first algorithm are more reliable.


%We exported the data relative to the objects in an excel file inside the
%project directory, please see this file and take a look at both first and
%second sheet for every image.

clear all
imagename = 'Moedas4.jpg';
img = imread(strcat('imgsrc/',imagename)); 
gray_img = rgb2gray(img);

%image segmentation using kmeans algorithm : separation of the background
%and the forground 
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


% 
% %detect edges 
% [~, threshold] = edge(gray_img,'canny');
% fudgeFactor = 6.5;
% BWedge = edge(gray_img,'canny',fudgeFactor *threshold);
% 
% %remove very small objects
% BWclean = bwareaopen(BWedge,45);
% 
% %perform dilate operation
% se90 = strel('line',2,0);
% se0 = strel('line',2,90);
% BWdilate = imdilate(BWclean,[se90 se0]);
% 
% %fill the holes
% BWfilled = imfill(BWdilate,'holes');


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
stats = regionprops('table',BW,'Area','Centroid','MajorAxisLength',...
 'MinorAxisLength','Perimeter','Image','BoundingBox');
BoundingBox = stats.BoundingBox;
Centroid = stats.Centroid;
Area = stats.Area;
Perimeter = stats.Perimeter;
Image = stats.Image;
Circularity = min(4*pi* Area./(Perimeter.^2),1);

%Table containing stats on each object (sorted by Area) 
rowNames = strings(N,1);
for i=1:N
    rowNames(i) = sprintf('Object%d',i);
end

tblsorted = sortrows(table(Area,Perimeter,Circularity,Centroid,BoundingBox));
TABLE = tblsorted;
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


%plot Centroids, Bounding boxes, Areas and Perimeters
figure
imshow(img)
hold on 
plot(Centroid(:,1),Centroid(:,2),'b*')

 for k = 1 : size(tblsorted)
  thisBB = tblsorted.BoundingBox(k,:);
  rectangle('Position', [thisBB(1),thisBB(2),thisBB(3),thisBB(4)],...
  'EdgeColor','g','LineWidth',1 );
    str = sprintf('%s\nArea: %0.2f\nPerimeter: %0.2f',tblsorted.Row{k},...
        tblsorted.Area(k),tblsorted.Perimeter(k));
  text(thisBB(1)+3,thisBB(2)+thisBB(4)-35,str,'FontSize',8,'Color','g');
  
 end
 
 
 tblrelDist = array2table(relDist,'RowNames',rowNames,'VariableNames',rowNames);

tblrelDist
 
filename = strcat('ImageObjectStatistics_',imagename,'.xlsx');
writetable(tblsorted,filename,'Sheet',1,'Range','B2','WriteRowNames',true);
fclose('all');
writetable(tblrelDist,filename,'Sheet',2,'Range','B2');
fclose('all');

 
 %Object selection and image order relative to Area 
 [c,r,~] = impixel;
 hold off
 mask = bwselect(BW,c,r,8);
 mask3 = cat(3, mask, mask, mask);
 J  = img;
 J(~mask3) = 0;
 %J = regionfill(gray_img,~mask);
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
 %J = regionfill(gray_img,~mask);
 I = [I J];   
    
 end
 
 
 figure
 imshow(I)


% Compute Sum of money
sum = 0;
MARGIN = 300;
for n = 1:size(TABLE)
    circ = TABLE{n, 3};
    if circ > 0.99;
        area = TABLE{n, 1};
        if and(area > 10364 - MARGIN, area < 10364 + MARGIN)
            sum = sum + 0.01;
        end
        if and(area > 14209 - MARGIN, area < 14209 + MARGIN)
            sum = sum + 0.02;
        end
        if and(area > 17656 - MARGIN, area < 17656 + MARGIN)
            sum = sum + 0.05;
        end
        if and(area > 15335 - MARGIN, area < 15335 + MARGIN)
            sum = sum + 0.10;
        end
        if and(area > 19513 - MARGIN, area < 19513 + MARGIN)
            sum = sum + 0.20;
        end
        if and(area > 23173 - MARGIN, area < 23173 + MARGIN)
            sum = sum + 0.50;
        end
        if and(area > 21363 - MARGIN, area < 21363 + MARGIN)
            sum = sum + 1;
        end
    end    
end

sum
