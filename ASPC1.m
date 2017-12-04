function [cellfinal,results] = ASPC1()
%%

%Going through each image seperately for image processing

i = 2;

imname = strcat('./ASPC1/Picture',num2str(i),'.png');
reader1 = bfGetReader(imname);
iplane = reader1.getIndex(1-1,1-1,1-1)+1;
img = bfGetPlane(reader1,iplane);

img2 = im2double(img);
img_dilate = imdilate(img2,strel('disk',3));
img_norm = img2./img_dilate;

mask = img_norm < 0.95;

edgem = edge(mask,'sobel');
edgem = imdilate(edgem,strel('disk',2));
mask = imdilate(edgem,strel('disk',5));


mask = bwmorph(mask,'clean');
mask = imdilate(mask,strel('disk',1));
mask = imclose(mask,strel('disk',1));
mask = imdilate(mask,strel('disk',1));
mask = imerode(mask,strel('disk',1));

measurements1 = regionprops( edgem,img_norm,'EulerNumber', 'ConvexArea', 'BoundingBox','MeanIntensity','Area','Perimeter','Eccentricity','EquivDiameter','WeightedCentroid','Orientation','ConvexArea','Extent','Centroid');
measurements1 = struct2table(measurements1);
measurements1 = measurements1(measurements1.Area > 10,:);


figure;
imshowpair(mask,img);
hold on;
scatter(measurements1.WeightedCentroid(:,1),measurements1.WeightedCentroid(:,2));


measurements = regionprops( mask,img_norm,'EulerNumber', 'ConvexArea', 'BoundingBox','MeanIntensity','Area','Perimeter','Eccentricity','EquivDiameter','WeightedCentroid','Orientation','ConvexArea','Extent','Centroid');
measurements = struct2table(measurements);

[o,~] = size(measurements);
[o2,~] = size(measurements1);
cents = zeros(o,1);
mindist = zeros(o2,1);

for i = 1:o
counter = 0;
box = measurements.BoundingBox(i,:,:,:);    
for l = 1:o2
    dist = sqrt((measurements1.WeightedCentroid(l,1)-measurements1.WeightedCentroid(:,1)).^2+(measurements1.WeightedCentroid(l,2)-measurements1.WeightedCentroid(:,2)).^2);
    dist = dist(dist > 0,:);
    mindistance = (dist);
    [x1,y1] = size(img);
    mindist(l) = min(mindistance)/(x1*y1);
    temp = ((box(1) < measurements1.WeightedCentroid(l,1))  & ...
        (measurements1.WeightedCentroid(l,1) <  (box(1)+ box(4)) &...
        (box(2) < measurements1.WeightedCentroid(l,2)) &...
        (measurements1.WeightedCentroid(l,2)) < box(2)+ box(3)));
    counter = counter + temp;
end
    cents(i) = counter;
    
end

measurements.nuclei = cents;
measurements.meancents = measurements.nuclei./measurements.Area;

%cells = {'MeanConnections','MinCentroidDist','MeanCentroidCluster','PercentsCentsInClusters'};
con = bwconncomp(mask);
cells(1,1:5) = {con.Connectivity/o,mean(mindist),max(mean(measurements.meancents)),sum(measurements.nuclei)/o2,'mesenchymal'};
resultmes = measurements;



%%

i = 4;
imname = strcat('./ASPC1/Picture',num2str(i),'.png');
reader1 = bfGetReader(imname);
iplane = reader1.getIndex(1-1,1-1,1-1)+1;
img = bfGetPlane(reader1,iplane);

img2 = im2double(img);
img_dilate = imdilate(img2,strel('disk',3));
img_norm = img2./img_dilate;

mask = img_norm < 0.95;

edgem = edge(mask,'sobel');
edgem = imdilate(edgem,strel('disk',1));
mask = imdilate(edgem,strel('disk',5));


mask = bwmorph(mask,'clean');
mask = imdilate(mask,strel('disk',1));
mask = imclose(mask,strel('disk',1));
mask = imdilate(mask,strel('disk',1));
mask = imerode(mask,strel('disk',1));


measurements1 = regionprops( edgem,img_norm,'EulerNumber', 'ConvexArea', 'BoundingBox','MeanIntensity','Area','Perimeter','Eccentricity','EquivDiameter','WeightedCentroid','Orientation','ConvexArea','Extent','Centroid');
measurements1 = struct2table(measurements1);

measurements1 = measurements1(measurements1.Area > 10,:);

figure;
imshowpair(mask,img);
hold on;
scatter(measurements1.WeightedCentroid(:,1),measurements1.WeightedCentroid(:,2));


measurements = regionprops( mask,img_norm,'EulerNumber', 'ConvexArea', 'BoundingBox','MeanIntensity','Area','Perimeter','Eccentricity','EquivDiameter','WeightedCentroid','Orientation','ConvexArea','Extent','Centroid');
measurements = struct2table(measurements);

[o,~] = size(measurements);
[o2,~] = size(measurements1);
cents = zeros(o,1);
mindist = zeros(o2,1);

for i = 1:o
counter = 0;
box = measurements.BoundingBox(i,:,:,:);    
for l = 1:o2
    dist = sqrt((measurements1.WeightedCentroid(l,1)-measurements1.WeightedCentroid(:,1)).^2+(measurements1.WeightedCentroid(l,2)-measurements1.WeightedCentroid(:,2)).^2);
    dist = dist(dist > 0,:);
    mindistance = (dist);
    [x1,y1] = size(img);
    mindist(l) = min(mindistance)/(x1*y1);
    temp = ((box(1) < measurements1.WeightedCentroid(l,1))  & ...
        (measurements1.WeightedCentroid(l,1) <  (box(1)+ box(4)) &...
        (box(2) < measurements1.WeightedCentroid(l,2)) &...
        (measurements1.WeightedCentroid(l,2)) < box(2)+ box(3)));
    counter = counter + temp;
end
    cents(i) = counter;
    
end

measurements.nuclei = cents;
measurements.meancents = measurements.nuclei./measurements.Area;

%cells = {'MeanConnections','MinCentroidDist','MeanCentroidCluster','PercentsCentsInClusters'};
con = bwconncomp(mask);
cells(2,1:5) = {con.Connectivity/o,mean(mindist),max(mean(measurements.meancents)),sum(measurements.nuclei)/o2,'mesenchymal'};
resultmes = vertcat(resultmes,measurements);


%%
i = 5;
imname = strcat('./ASPC1/Picture',num2str(i),'.png');
reader1 = bfGetReader(imname);
iplane = reader1.getIndex(1-1,1-1,1-1)+1;
img = bfGetPlane(reader1,iplane);

img2 = im2double(img);
img_dilate = imdilate(img2,strel('disk',3));
img_norm = img2./img_dilate;

mask = img_norm < 0.95;

edgem = edge(mask,'sobel');
edgem = imdilate(edgem,strel('disk',1));
edgem = imerode(edgem,strel('disk',1));

mask = imdilate(edgem,strel('disk',5));


mask = bwmorph(mask,'clean');
mask = imdilate(mask,strel('disk',1));
mask = imclose(mask,strel('disk',1));
mask = imdilate(mask,strel('disk',1));
mask = imerode(mask,strel('disk',1));


%mask = imfill(mask,'holes');

measurements1 = regionprops( edgem,img_norm,'EulerNumber', 'ConvexArea', 'BoundingBox','MeanIntensity','Area','Perimeter','Eccentricity','EquivDiameter','WeightedCentroid','Orientation','ConvexArea','Extent','Centroid');
measurements1 = struct2table(measurements1);

measurements1 = measurements1(measurements1.Area > 10,:);
figure;
imshowpair(mask,img);
hold on;
scatter(measurements1.WeightedCentroid(:,1),measurements1.WeightedCentroid(:,2));

measurements = regionprops( mask,img_norm,'EulerNumber', 'ConvexArea', 'BoundingBox','MeanIntensity','Area','Perimeter','Eccentricity','EquivDiameter','WeightedCentroid','Orientation','ConvexArea','Extent','Centroid');
measurements = struct2table(measurements);

[o,~] = size(measurements);
[o2,~] = size(measurements1);
cents = zeros(o,1);
mindist = zeros(o2,1);

for i = 1:o
counter = 0;
box = measurements.BoundingBox(i,:,:,:);    
for l = 1:o2
    dist = sqrt((measurements1.WeightedCentroid(l,1)-measurements1.WeightedCentroid(:,1)).^2+(measurements1.WeightedCentroid(l,2)-measurements1.WeightedCentroid(:,2)).^2);
    dist = dist(dist > 0,:);
    mindistance = (dist);
    [x1,y1] = size(img);
    mindist(l) = min(mindistance)/(x1*y1);
    temp = ((box(1) < measurements1.WeightedCentroid(l,1))  & ...
        (measurements1.WeightedCentroid(l,1) <  (box(1)+ box(4)) &...
        (box(2) < measurements1.WeightedCentroid(l,2)) &...
        (measurements1.WeightedCentroid(l,2)) < box(2)+ box(3)));
    counter = counter + temp;
end
    cents(i) = counter;
    
end

measurements.nuclei = cents;
measurements.meancents = measurements.nuclei./measurements.Area;

%cells = {'MeanConnections','MinCentroidDist','MeanCentroidCluster','PercentsCentsInClusters'};
con = bwconncomp(mask);
cells(3,1:5) = {con.Connectivity/o,mean(mindist),max(mean(measurements.meancents)),sum(measurements.nuclei)/o2,'mesenchymal'};
resultmes = vertcat(resultmes,measurements);


%%

i = 6;
imname = strcat('./ASPC1/Picture',num2str(i),'.png');
reader1 = bfGetReader(imname);
iplane = reader1.getIndex(1-1,1-1,1-1)+1;
img = bfGetPlane(reader1,iplane);

img2 = im2double(img);
img_dilate = imdilate(img2,strel('disk',3));
img_norm = img2./img_dilate;

mask = img_norm < 0.90;

edgem = edge(mask,'canny');
edgem = imdilate(edgem,strel('disk',1));


mask = imdilate(edgem,strel('disk',5));


mask = bwmorph(mask,'clean');
mask = imdilate(mask,strel('disk',1));
mask = imclose(mask,strel('disk',1));
mask = imdilate(mask,strel('disk',1));
mask = imerode(mask,strel('disk',1));


%mask = imfill(mask,'holes');

measurements1 = regionprops( edgem,img_norm,'EulerNumber', 'ConvexArea', 'BoundingBox','MeanIntensity','Area','Perimeter','Eccentricity','EquivDiameter','WeightedCentroid','Orientation','ConvexArea','Extent','Centroid');
measurements1 = struct2table(measurements1);

measurements1 = measurements1(measurements1.Area > 10,:);


figure;
imshowpair(mask,img);
hold on;
scatter(measurements1.WeightedCentroid(:,1),measurements1.WeightedCentroid(:,2));
measurements = regionprops( mask,img_norm,'EulerNumber', 'ConvexArea', 'BoundingBox','MeanIntensity','Area','Perimeter','Eccentricity','EquivDiameter','WeightedCentroid','Orientation','ConvexArea','Extent','Centroid');
measurements = struct2table(measurements);

[o,~] = size(measurements);
[o2,~] = size(measurements1);
cents = zeros(o,1);
mindist = zeros(o2,1);

for i = 1:o
counter = 0;
box = measurements.BoundingBox(i,:,:,:);    
for l = 1:o2
    dist = sqrt((measurements1.WeightedCentroid(l,1)-measurements1.WeightedCentroid(:,1)).^2+(measurements1.WeightedCentroid(l,2)-measurements1.WeightedCentroid(:,2)).^2);
    dist = dist(dist > 0,:);
    mindistance = (dist);
    [x1,y1] = size(img);
    mindist(l) = min(mindistance)/(x1*y1);
    temp = ((box(1) < measurements1.WeightedCentroid(l,1))  & ...
        (measurements1.WeightedCentroid(l,1) <  (box(1)+ box(4)) &...
        (box(2) < measurements1.WeightedCentroid(l,2)) &...
        (measurements1.WeightedCentroid(l,2)) < box(2)+ box(3)));
    counter = counter + temp;
end
    cents(i) = counter;
    
end

measurements.nuclei = cents;
measurements.meancents = measurements.nuclei./measurements.Area;

%cells = {'MeanConnections','MinCentroidDist','MeanCentroidCluster','PercentsCentsInClusters'};
con = bwconncomp(mask);
cells(4,1:5) = {con.Connectivity/o,mean(mindist),max(mean(measurements.meancents)),sum(measurements.nuclei)/o2,'mesenchymal'};
resultmes = vertcat(resultmes,measurements);


%%

i = 7;
imname = strcat('./ASPC1/Picture',num2str(i),'.png');
reader1 = bfGetReader(imname);
iplane = reader1.getIndex(1-1,1-1,1-1)+1;
img = bfGetPlane(reader1,iplane);

img2 = im2double(img);
img_dilate = imdilate(img2,strel('disk',3));
img_norm = img2./img_dilate;

mask = img_norm < 0.92;

edgem = edge(mask,'canny');
edgem = imdilate(edgem,strel('disk',1));
edgem = imerode(edgem,strel('disk',1));
%edgem = imclose(edgem,strel('disk',1));

mask = imdilate(edgem,strel('disk',5));


mask = bwmorph(mask,'clean');
mask = imdilate(mask,strel('disk',1));
mask = imclose(mask,strel('disk',1));
mask = imdilate(mask,strel('disk',1));
mask = imerode(mask,strel('disk',1));


%mask = imfill(mask,'holes');

measurements1 = regionprops( edgem,img_norm,'EulerNumber', 'ConvexArea', 'BoundingBox','MeanIntensity','Area','Perimeter','Eccentricity','EquivDiameter','WeightedCentroid','Orientation','ConvexArea','Extent','Centroid');
measurements1 = struct2table(measurements1);

measurements1 = measurements1(measurements1.Area > 10,:);

figure;
imshowpair(mask,img);
hold on;
scatter(measurements1.WeightedCentroid(:,1),measurements1.WeightedCentroid(:,2));
measurements = regionprops( mask,img_norm,'EulerNumber', 'ConvexArea', 'BoundingBox','MeanIntensity','Area','Perimeter','Eccentricity','EquivDiameter','WeightedCentroid','Orientation','ConvexArea','Extent','Centroid');
measurements = struct2table(measurements);


[o,~] = size(measurements);
[o2,~] = size(measurements1);
cents = zeros(o,1);
mindist = zeros(o2,1);

for i = 1:o
counter = 0;
distance = inf;
box = measurements.BoundingBox(i,:,:,:);    
for l = 1:o2
    dist = sqrt((measurements1.WeightedCentroid(l,1)-measurements1.WeightedCentroid(:,1)).^2+(measurements1.WeightedCentroid(l,2)-measurements1.WeightedCentroid(:,2)).^2);
    dist = dist(dist > 0,:);
    mindistance = (dist);
    [x1,y1] = size(img);
    mindist(l) = min(mindistance)/(x1*y1);
    temp = ((box(1) < measurements1.WeightedCentroid(l,1))  & ...
        (measurements1.WeightedCentroid(l,1) <  (box(1)+ box(4)) &...
        (box(2) < measurements1.WeightedCentroid(l,2)) &...
        (measurements1.WeightedCentroid(l,2)) < box(2)+ box(3)));
    counter = counter + temp;
end
    cents(i) = counter;
    
end

measurements.nuclei = cents;
measurements.meancents = measurements.nuclei./measurements.Area;

%cells = {'MeanConnections','MinCentroidDist','MeanCentroidCluster','PercentsCentsInClusters'};
con = bwconncomp(mask);
cells(5,1:5) = {con.Connectivity/o,mean(mindist),max(mean(measurements.meancents)),sum(measurements.nuclei)/o2,'mesenchymal'};
resultmes = vertcat(resultmes,measurements);


%Combining the total data together 

[x,~] = size(resultmes);
[fakecolumn{1:x}] = deal('Mesenchymal');
fakecolumn = cell2table(fakecolumn');
resultmes.class = fakecolumn.Var1;
resultmes.ratio = resultmes.Area./resultmes.Perimeter;
resultmes = resultmes(resultmes.Area > 10,:);
results = resultmes;
clear fakecolumn; 

cellfinal = cells;

end