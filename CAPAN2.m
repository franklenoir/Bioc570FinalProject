function [cellfinal,results] = CAPAN2()

%%

i = 1;
imname = strcat('./capan2/Picture',num2str(i),'.png');
reader1 = bfGetReader(imname);
iplane = reader1.getIndex(1-1,1-1,1-1)+1;
img = bfGetPlane(reader1,iplane);

img2 = im2double(img);
img_dilate = imdilate(img2,strel('disk',3));
img_norm = img2./img_dilate;

mask = img_norm < 0.95;
mask = imopen(mask,strel('disk',2));

edgem = edge(mask,'sobel');

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
resultmes = measurements;
%cells = {'MeanConnections','MinCentroidDist','MeanCentroidCluster','PercentsCentsInClusters'};
con = bwconncomp(mask);
cells(1:5) = {con.Connectivity/o,mean(mindist),max(mean(measurements.meancents)),sum(measurements.nuclei)/o2,'Epithelial'};



%%

i = 5;
imname = strcat('./capan2/Picture',num2str(i),'.png');
reader1 = bfGetReader(imname);
iplane = reader1.getIndex(1-1,1-1,1-1)+1;
img = bfGetPlane(reader1,iplane);

img2 = im2double(img);
img_dilate = imdilate(img2,strel('disk',3));
img_norm = img2./img_dilate;

mask = img_norm < 0.95;
mask = imopen(mask,strel('disk',2));

edgem = edge(mask,'sobel');

mask = bwmorph(mask,'clean');
mask = imopen(mask,strel('arbitrary',3));
mask = imclose(mask,strel('disk',5));
mask = edge(mask,'sobel',[]);
mask = imdilate(mask,strel('disk',2));
mask = imfill(mask,'holes');
mask = imerode(mask,strel('disk',2));
mask = bwmorph(mask,'clean');

mask = imdilate(mask,strel('disk',15));
mask = imfill(mask,'holes');
mask = imerode(mask,strel('disk',10));

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
cells(2,1:5) = {con.Connectivity/o,mean(mindist),max(mean(measurements.meancents)),sum(measurements.nuclei)/o2,'Epithelial'};


resultmes = vertcat(resultmes,measurements);

%%

i = 6;
imname = strcat('./capan2/Picture',num2str(i),'.png');
reader1 = bfGetReader(imname);
iplane = reader1.getIndex(1-1,1-1,1-1)+1;
img = bfGetPlane(reader1,iplane);

img2 = im2double(img);
img_dilate = imdilate(img2,strel('disk',3));
img_norm = img2./img_dilate;

mask = img_norm < 0.98;

mask = imopen(mask,strel('disk',2));

edgem = edge(mask,'sobel');


mask = bwmorph(mask,'clean');
mask = imopen(mask,strel('arbitrary',3));
mask = imclose(mask,strel('disk',5));
mask = edge(mask,'sobel',[]);
mask = imdilate(mask,strel('disk',2));
mask = imfill(mask,'holes');
mask = imerode(mask,strel('disk',2));
mask = bwmorph(mask,'clean');

mask = imdilate(mask,strel('disk',15));
mask = imfill(mask,'holes');
mask = imerode(mask,strel('disk',10));

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
cells(3,1:5) = {con.Connectivity/o,mean(mindist),max(mean(measurements.meancents)),sum(measurements.nuclei)/o2,'Epithelial'};


resultmes = vertcat(resultmes,measurements);


%%

i = 7;
imname = strcat('./capan2/Picture',num2str(i),'.png');
reader1 = bfGetReader(imname);
iplane = reader1.getIndex(1-1,1-1,1-1)+1;
img = bfGetPlane(reader1,iplane);

img2 = im2double(img);
img_dilate = imdilate(img2,strel('disk',2));
img_norm = img2./img_dilate;
mask = img_norm < 0.991;

mask = imopen(mask,strel('disk',2));

edgem = edge(mask,'sobel');

mask = bwmorph(mask,'clean');
mask = imopen(mask,strel('disk',3));
mask = imfill(mask,26,'holes');
mask = edge(mask,'canny',[]);
mask = imdilate(mask,strel('disk',15));
mask = imfill(mask,'holes');
mask = imerode(mask,strel('disk',8));

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
cells(4,1:5) = {con.Connectivity/o,mean(mindist),max(mean(measurements.meancents)),sum(measurements.nuclei)/o2,'Epithelial'};


resultmes = vertcat(resultmes,measurements);

%%

i = 8;
imname = strcat('./capan2/Picture',num2str(i),'.png');
reader1 = bfGetReader(imname);
iplane = reader1.getIndex(1-1,1-1,1-1)+1;
img = bfGetPlane(reader1,iplane);

img2 = im2double(img);
img_dilate = imdilate(img2,strel('disk',2));
img_norm = img2./img_dilate;
mask = img_norm < 0.991;

mask = imopen(mask,strel('disk',2));

edgem = edge(mask,'canny');
%edgem = imclose(mask,strel('disk',1));


mask = bwmorph(mask,'clean');
mask = imopen(mask,strel('disk',3));
mask = imfill(mask,26,'holes');
mask = edge(mask,'canny',[]);
mask = imdilate(mask,strel('disk',10));
mask = imfill(mask,'holes');
mask = imerode(mask,strel('disk',6));

measurements1 = regionprops( edgem,img_norm,'EulerNumber', 'ConvexArea', 'BoundingBox','MeanIntensity','Area','Perimeter','Eccentricity','EquivDiameter','WeightedCentroid','Orientation','ConvexArea','Extent','Centroid');
measurements1 = struct2table(measurements1);
measurements1 = measurements1(measurements1.Area > 25,:);

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
cells(5,1:5) = {con.Connectivity/o,mean(mindist),max(mean(measurements.meancents)),sum(measurements.nuclei)/o2,'Epithelial'};


resultmes = vertcat(resultmes,measurements);

%%
[x,~] = size(resultmes);
[fakecolumn{1:x}] = deal('Epithelial');
fakecolumn = cell2table(fakecolumn');
resultmes.class = fakecolumn.Var1;
resultmes.ratio = resultmes.Area./resultmes.Perimeter;
resultmes = resultmes(resultmes.Area > 10,:);
results = resultmes;
clear fakecolumn; 
cellfinal = cells;

end