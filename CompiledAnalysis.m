% Training Data


%Running Training Date Together
[cells1,aspc1] = ASPC1();
[cells2,capan2] = CAPAN2();
[cells3,panc1] = Panc1();

%save('training.mat','aspc1','capan2','panc1','cells1','cells2','cells3');
%Saved data 

load('training.mat');

%combine data for analysis. Cells is individual picture analysis, resultmes
%is each regionprops result. 
resultmes = vertcat(aspc1,capan2,panc1);
resultcells = vertcat(cells1,cells2,cells3);
resultcells = cell2table(resultcells);

%Generate normalized metrics (trying to avoid image bias)
resultmes.NormArea = resultmes.ConvexArea./resultmes.Area;
resultmes.meanareapernucs = resultmes.Area./resultmes.nuclei;

%Deleting unnecessary columns
resultmes.BoundingBox = [];
resultmes.Centroid=[];
resultmes.WeightedCentroid=[];


%save('Model.mat',trainedModel);
%Resulting Trained Model 

%Test Data
[cells4,patc53] = Patc53();
[cells5,patc69] = Patc69();
[cells6,patc124] = Patc124();

%save('Test.mat')
load('Test.mat');

%Generating the same columns as above
patc53.NormArea = patc53.ConvexArea./patc53.Area;
patc53.meanareapernucs = patc53.Area./patc53.nuclei;
patc69.NormArea = patc69.ConvexArea./patc69.Area;
patc69.meanareapernucs = patc69.Area./patc69.nuclei;
patc124.NormArea = patc124.ConvexArea./patc124.Area;
patc124.meanareapernucs = patc124.Area./patc124.nuclei;

load('Model.mat'); % first pass at model

%Running Trained Model on test with results (m vs e)
yfit53 = trainedModel.predictFcn(patc53); %280 - e, 143 - m
yfit69 = trainedModel.predictFcn(patc69); %610 - e, 205 - m
yfit124 = trainedModel.predictFcn(patc124); %227 - e, 59 - m

% Code Used to count trainedModel
%[uniqueXX, ~, J]=unique(yfit53);
%occ = histc(J, 1:numel(uniqueXX))


% Ran Classifcation learner again on training data using only columns 
% selected by NCA analysis.
save('trainedModel2.mat','trainedModel2');
load('trainedModel2.mat');

%Selected columns 9,11,12,14 for test data
 patc53_2 = [patc53.MeanIntensity,patc53.meancents,patc53.ratio,patc53.meanareapernucs];
 patc69_2 = [patc69.MeanIntensity,patc69.meancents,patc69.ratio,patc69.meanareapernucs];
 patc124_2 = [patc124.MeanIntensity,patc124.meancents,patc124.ratio,patc124.meanareapernucs];

%improved
yfit53_2 = trainedModel2.predictFcn(patc53_2); %253 - e, 170 - m
yfit69_2 = trainedModel2.predictFcn(patc69_2); %734 - e, 81 - m
yfit124_2 = trainedModel2.predictFcn(patc124_2); %254 - e, 32 - m

%[uniqueXX, ~, J]=unique(yfit124_2);
%occ = histc(J, 1:numel(uniqueXX));

%Analysis on Image breakdown (per image rather than per individual cell)
%save('trainedcells.mat','trained_cells');
load('trainedcells.mat');

%training data had 75% epithelial selection

cells4 = cell2table(cells4);
cells4.Properties.VariableNames = {'resultcells1'  'resultcells2'  'resultcells3'  'resultcells4' 'class'};
cells5 = cell2table(cells5);
cells5.Properties.VariableNames = {'resultcells1'  'resultcells2'  'resultcells3'  'resultcells4' 'class'};
cells6 = cell2table(cells6);
cells6.Properties.VariableNames = {'resultcells1'  'resultcells2'  'resultcells3'  'resultcells4' 'class'};


yfit53cells = trained_cells.predictFcn(cells4); %100% mes
yfit69cells = trained_cells.predictFcn(cells5); %100% mes
yfit124cells = trained_cells.predictFcn(cells6); % 100% mes

%All results mes 

%[uniqueXX, ~, J]=unique(yfit53cells);
%occ = histc(J, 1:numel(uniqueXX));

