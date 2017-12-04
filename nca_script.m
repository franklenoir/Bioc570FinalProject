%%%NCA 
clear all;


%%creating the model
load('/training.mat');

resultmes = vertcat(aspc1,panc1, capan2);
resultmes.BoundingBox = [];
resultmes.Centroid=[];
resultmes.WeightedCentroid=[];
resultcells = vertcat(cells1,cells2,cells3);
resultcells = cell2table(resultcells);
resultmes.NormArea = resultmes.ConvexArea./resultmes.Area;
resultmes.meanareapernucs = resultmes.Area./resultmes.nuclei;
temp = resultmes;
temp1 = temp(temp.meanareapernucs  < inf,:);

temp2 = temp1.class;

temp1.class = [];
temp1 = table2array(temp1);

for ii = 1:length(temp2)
    if temp2(ii)=="Epithelial"
        temp3(ii)=1;
    else
        temp3(ii)=-1;
    end
end

temp3=temp3';

test = fscnca(temp1,temp3);

figure();
plot(test.FeatureWeights,'ro');
grid on;
xlabel('Feature index');
ylabel('Feature weight');

%test_r=fsrnca(temp1, temp3, 'FitMethod', 'exact');

%figure();
%plot(test_r.FeatureWeights,'ro');
%grid on;
%xlabel('Feature index');
%ylabel('Feature weight');



%% pulling in test data

%after running base nca, get the L value without any fitting first

load('/Test.mat');

%creating a combined matrix of all patc lines
patc=vertcat(patc53, patc69, patc124);
patc_x=patc;
patc_x.NormArea = patc_x.ConvexArea./patc_x.Area;
patc_x.meanareapernucs = patc_x.Area./patc_x.nuclei;
patc_x=patc_x(patc_x.meanareapernucs  < inf,:);
patc_class=patc_x.class;
patc_x.BoundingBox=[];
patc_x.Centroid=[];
patc_x.WeightedCentroid=[];
patc_x.class=[];
patc_x=table2array(patc_x);


% epithelial coded as 1, mesenchymal coded as -1 
for ii = 1:length(patc_class)
    if patc_class(ii)=="Epithelial"
        patc_y(ii)=1;
    else
        patc_y(ii)=-1;
    end
end
patc_y=patc_y';

L=loss(test, patc_x, patc_y);


%% correcting for high L values

%checking to see what the impact of fitting actually is

nofit=fscnca(temp1,temp3, 'FitMethod', 'none');

L_nofit = loss(nofit, patc_x, patc_y);

L_zerolambda=fscnca(temp1,temp3,'FitMethod','exact');
L_lambda = loss(L_zerolambda,patc_x,patc_y);


%identifying best lambda in a loop

cvp = cvpartition(temp3,'kfold',5);
numvalidsets = cvp.NumTestSets;

n = length(temp3);
lambdavals = linspace(0,20,20)/n;
lossvals = zeros(length(lambdavals),numvalidsets);

for i = 1:length(lambdavals)
    for k = 1:numvalidsets
        X = temp1(cvp.training(k),:);
        y = temp3(cvp.training(k),:);
        Xtest = temp1(cvp.test(k),:);
        ytest = temp3(cvp.test(k),:);

        nca = fscnca(X,y,'FitMethod','exact', ...
        'Solver','sgd', 'Lambda',lambdavals(i), ...
             'IterationLimit',30,'GradientTolerance',1e-4, ...
             'Standardize',true);

        lossvals(i,k) = loss(nca,Xtest,ytest,'LossFunction','classiferror');
    end
end

meanloss = mean(lossvals,2);

figure()
plot(lambdavals,meanloss,'ro-')
xlabel('Lambda')
ylabel('Loss (MSE)')
grid on

[~,idx] = min(meanloss) % Find the index
bestlambda = lambdavals(idx) % Find the best lambda value
bestloss = meanloss(idx)

nca_best = fscnca(temp1,temp3,'FitMethod','exact','Solver','sgd',...
    'Lambda',bestlambda,'Standardize',true,'Verbose',1);

figure();
plot(nca_best.FeatureWeights,'ro');
grid on;
xlabel('Feature index');
ylabel('Feature weight');

tol    = 0.02;
selidx = find(nca_best.FeatureWeights > tol*max(1,max(nca_best.FeatureWeights)))

temp1=temp1(:,selidx);

%svmMdl = fitcsvm(temp1,temp3);

%L = loss(svmMdl,Xtest(:,selidx),ytest);

%% t test on features selected

temp_epi=resultmes(resultmes.class=="Epithelial",[9 10 13 15]);
temp_mes=resultmes(resultmes.class=="Mesenchymal",[9 10 13 15]);
temp_epi = temp_epi(temp_epi.meanareapernucs  < inf,:);
temp_mes = temp_mes(temp_mes.meanareapernucs  < inf,:);


col_var=vartest2(table2array(temp_epi), table2array(temp_mes));

%all but column 3 have equal variance

for i = 1:4
    col_var(i)=vartest2(table2array(temp_epi(:,i)), table2array(temp_mes(:,i)));
    if col_var(i) == 1
        [h(i), p(i)] = ttest2(table2array(temp_epi(:,i)), table2array(temp_mes(:,i)));
    else
        [h(i), p(i)] = ttest2(table2array(temp_epi(:,i)), table2array(temp_mes(:,i)), 'Vartype','unequal');
    end
end
        
        
