%% Load Experiments

[loadList, loadPath ]= uigetfile('Z:\ioldenburg\outputdata','MultiSelect','on');
%%
numExps = numel(loadList);

clear All
for ind = 1:numExps
    pTime =tic;
    fprintf(['Loading Experiment ' num2str(ind) '...']);
    All(ind) = load(fullfile(loadPath,loadList{ind}),'out');
    fprintf([' Took ' num2str(toc(pTime)) 's.\n'])
end

%% Known Error Manual Catching
All(1).out.exp.stimParamsBackup = All(1).out.exp.stimParams;
All(1).out.exp.holoTargetsBackup = All(1).out.exp.holoTargets;
%%
All(1).out.exp.stimParams.Seq([3 4])=[];
All(1).out.exp.stimParams.numPulse([3 4])=[];
All(1).out.exp.stimParams.roi([3 4])=[];
All(1).out.exp.stimParams.Hz([2 3])=[];
All(1).out.exp.stimParams.numCells([2 3])=[];
All(1).out.exp.holoTargets([2 3])=[];



%% Clean Data and stuff i dunno come up with a better title someday

FRDefault=6;
recWinSec = [1.25 2.5];

 for ind =1:numExps
     pTime =tic;
     fprintf(['Processing Experiment ' num2str(ind) '...']);
     
     All(ind).out.anal.numCells = size(All(ind).out.exp.zdfData,1);
     numCells(ind) = size(All(ind).out.exp.zdfData,1);
     
     if ~isfield(All(ind).out.info,'FR')
         All(ind).out.info.FR=FRDefault;
     end
     
     sz = size(All(ind).out.exp.zdfData);
     
     winToUse = min(round(recWinSec*All(ind).out.info.FR),[inf sz(2)]) ;
     rdata = squeeze(mean(All(ind).out.exp.zdfData(:,winToUse,:),2));
     bwinToUse = max(round([0 recWinSec(1)]*All(ind).out.info.FR),[1 1]);
     bdata = squeeze(mean(All(ind).out.exp.zdfData(:,bwinToUse,:),2));
     
     All(ind).out.exp.rdData=rdata;
     All(ind).out.exp.bdata=bdata;

     sz2 = size(All(ind).out.vis.zdfData);
     winToUse = min(round(recWinSec*All(ind).out.info.FR),[inf sz2(2)]) ;
     bwinToUse = max(round([0 recWinSec(1)]*All(ind).out.info.FR),[1 1]);

     rdata = squeeze(mean(All(ind).out.vis.zdfData(:,winToUse,:),2));
     bdata = squeeze(mean(All(ind).out.vis.zdfData(:,bwinToUse,:),2));
     
     All(ind).out.vis.rdata=rdata;
     All(ind).out.vis.bdata=bdata;
     
     temp = unique([All(ind).out.exp.holoTargets{:}]);
     temp(isnan(temp))=[];
     All(ind).out.anal.targets = temp;
     numUniqueTargets(ind) =numel(temp);
     
     %ensure has a visID
     if ~isfield(All(ind).out.exp,'visID')
         All(ind).out.exp.visID = ones(size(All(ind).out.exp.stimID));
         disp(['Added visID to Exp ' num2str(ind)]);
     end
     
     if numel(All(ind).out.vis.visID) ~= numel(All(ind).out.vis.lowMotionTrials)
         All(ind).out.vis.lowMotionTrials(end+1:numel(All(ind).out.vis.visID))= 0 ;
     end
     
     fprintf([' Took ' num2str(toc(pTime)) 's.\n'])

 end
%% Determine the OSI from the Vis section of each cell.

for ind=1:numExps
    pTime =tic;
    fprintf(['Processing Experiment ' num2str(ind) '...']);
    
    uVisID = unique(All(ind).out.vis.visID);
    uVisID(uVisID==0)=[];
    
    oriCurve=[];
    for i=1:numel(uVisID)
        v= uVisID(i);
 
        trialsToUse = All(ind).out.vis.visID==v & All(ind).out.vis.lowMotionTrials;
        
        oriCurve(i,:)=mean(All(ind).out.vis.rdata(:,trialsToUse),2);
    end
    
    All(ind).out.anal.oriCurve = oriCurve;
    
    [maxOriVal maxOriIndex]= max(oriCurve);
    All(ind).out.anal.prefOri = maxOriIndex;
    
    % (Rpref ? Rorth)/(Rpref + Rorth)
    % visID 1 is always catch (I Hope... Will is more confident...ish)
    
    prefOri = maxOriIndex;
    orthoOri = prefOri-2;
    orthoOri(orthoOri<2)=orthoOri(orthoOri<2)+8;
    
    orthoOri2 = orthoOri+4;
    orthoOri2(orthoOri2>9) = orthoOri2(orthoOri2>9)-8;
    
    orthoOri = cat(1,orthoOri, orthoOri2);
    
    
    
    %     orthoOri(prefOri==1)=NaN;
    
    oriCurveBL = oriCurve - min(oriCurve);%oriCurve(1,:);
    
    OSI=[];
    for i=1:numel(prefOri)
        OSI(i) = (oriCurveBL(prefOri(i),i) - mean(oriCurveBL(orthoOri(:,i)',i)) ) / (oriCurveBL(prefOri(i),i)+ mean(oriCurveBL(orthoOri(:,i)',i)) );
        %     OSI = (oriCurveBL(prefOri) - oriCurveBL(orthoOri) ) ./ ( oriCurveBL(prefOri)+oriCurveBL(orthoOri) )
        OSI(prefOri==1)=nan;
    end
    
    All(ind).out.anal.OSI=OSI;
    
    
    pVisR=[];pVisT=[];
    for i=1:All(ind).out.anal.numCells
        trialsToUse = All(ind).out.vis.visID~=0 & All(ind).out.vis.lowMotionTrials;
        pVisR(i) = anova1(All(ind).out.vis.rdata(i,trialsToUse),All(ind).out.vis.visID(trialsToUse),'off');
        
        trialsToUse = All(ind).out.vis.visID~=0 & All(ind).out.vis.visID~=1 & All(ind).out.vis.lowMotionTrials;
        pVisT(i) = anova1(All(ind).out.vis.rdata(i,trialsToUse),All(ind).out.vis.visID(trialsToUse),'off');
    end
    
    All(ind).out.anal.pVisR = pVisR;
    All(ind).out.anal.pVisT = pVisT;
    
    alpha = 0.05;
    
    All(ind).out.anal.visPercent = sum(pVisR<alpha) / numel(pVisR);
    visPercent(ind) =  All(ind).out.anal.visPercent;

    meanOSI=[];ensembleOSI=[];ensembleOriCurve =[];ensemblePref=[];
    
    for i=1:numel(All(ind).out.exp.holoTargets)
        ht = All(ind).out.exp.holoTargets{i};
        ht(isnan(ht))=[];
        meanOSI(i)=nanmean(OSI(ht));
        
        ensOriCurve = mean(oriCurve(:,ht),2);
        ensembleOriCurve(i,:)= ensOriCurve;
        [maxOriVal maxOriIndex]= max(ensOriCurve);
        prefOri = maxOriIndex;
        
        orthoOri = prefOri-2;
        orthoOri(orthoOri<2)=orthoOri(orthoOri<2)+8;
        orthoOri2 = orthoOri+4;
        orthoOri2(orthoOri2>9) = orthoOri2(orthoOri2>9)-8;
        orthoOri = cat(1,orthoOri, orthoOri2);
        oriCurveBL = ensOriCurve - min(ensOriCurve);
        
        ensembleOSI(i) = (oriCurveBL(prefOri)- mean(oriCurveBL(orthoOri'))) / (oriCurveBL(prefOri)+ mean(oriCurveBL(orthoOri')));
        ensemblePref(i) = prefOri;
    end
    
    deg = [nan 0:45:315];
    All(ind).out.anal.ensembleOSI=ensembleOSI;
    All(ind).out.anal.meanOSI=meanOSI;
    All(ind).out.anal.ensembleOSI=ensembleOSI;
    All(ind).out.anal.ensembleOriCurve=ensembleOriCurve;
    All(ind).out.anal.ensemblePrefOri=ensemblePref;
    All(ind).out.anal.ensemblePrefDeg=deg(ensemblePref);
    
    fprintf([' Took ' num2str(toc(pTime)) 's.\n'])
end

%% Pretty plots of OSI and tunings

clear allOSI ensOSI meanOSI ensNum roiNum
% OSI across all cells, all experiments
for i = 1:numel(All)
    allOSI{i} = All(i).out.anal.OSI(:);
    ensOSI{i} = All(i).out.anal.ensembleOSI(:);
    meanOSI{i} = All(i).out.anal.meanOSI(:);

    ensNum{i} = cellfun(@(x) sum(~isnan(x)),All(i).out.exp.holoTargets)'; %number of discovered Cells in ensemble
    roiNum{i} = cellfun(@(x) sum(~isnan(x)),All(i).out.exp.rois)'; %number of shot targets in ensemble
end

% unroll
allOSI = cell2mat(allOSI(:));
ensOSI = cell2mat(ensOSI(:));
meanOSI = cell2mat(meanOSI(:));
ensNum = cell2mat(ensNum(:));
roiNum = cell2mat(roiNum(:));

% plot
figure(1)
clf(1)
colors = {rgb('royalblue'), rgb('firebrick'), rgb('coral')};

hold on
h(1) = histogram(allOSI, 25);
h(2) = histogram(ensOSI, 25);
h(3) = histogram(meanOSI, 25);

for i=1:numel(h)
    h(i).Normalization = 'pdf';
    h(i).BinWidth = 0.05;
    h(i).FaceColor = colors{i};
    h(i).FaceAlpha = 0.44;
    kde(i) = fitdist(h(i).Data, 'kernel');
    p = plot(h(i).BinEdges, pdf(kde(i),h(i).BinEdges));
    p.LineWidth=2;
    p.Color= colors{i};
end

hold off

set(gcf(),'Name','OSI Across All Expts')
title('OSI Across All Expts')
xlabel('Orientation Selectivity Index')
ylabel('PDF')
legend('All Cells', 'Ensemble', 'Ensemble Mean')
 


%% now look at tuned ensembles
% get plots for the 2 different methods with higher bin count
f2 = figure(2);
clf(f2)
hold on
h(1) = histogram(ensOSI, 50);
%h(2) = histogram(meanOSI, 50);

% try using ensembleOSI > 0.3 for "tuned", this is arbitrary
OSIthreshold = 0.3;
for i = 1:numel(All)
    istuned = All(i).out.anal.ensembleOSI > OSIthreshold;
    tunedEnsembles = All(i).out.exp.holoTargets(istuned);
    untunedEnsembles = All(i).out.exp.holoTargets(~istuned);
    tunedEnsembleIdx = find(All(i).out.anal.ensembleOSI >= OSIthreshold);
    untunedEnsembleIdx = find(All(i).out.anal.ensembleOSI < OSIthreshold);
    All(i).out.anal.tunedEnsembles = tunedEnsembles;
    All(i).out.anal.untunedEnsembles = untunedEnsembles;
    All(i).out.anal.tunedEnsembleIdx = tunedEnsembleIdx;
    All(i).out.anal.untunedEnsembleIdx = untunedEnsembleIdx;
end

% plot the OSIs of tuned vs untuned ensembles
clear allOSIT ensOSIT meanOSIT
for i = 1:numel(All)
    allOSIT{i} = All(i).out.anal.OSI(:);
    ensOSIT{i} = All(i).out.anal.ensembleOSI(:);
    meanOSIT{i} = All(i).out.anal.meanOSI(:);
end

% unroll
allOSIT = cell2mat(allOSIT(:));
ensOSIT = cell2mat(ensOSIT(:));
meanOSIT = cell2mat(meanOSIT(:));

% f3 = subplots(1,2,1)

% and thier preferred oris

    
    
%% Get the number of spikes in each stimulus

clear numSpikesEachStim numCellsEachEns
for ind = 1:numExps
    temp = All(ind).out.exp.stimParams.numPulse;
    numSpikes=[];
    c=0;
    for i=1:numel(temp); %overly complicated way of aligning 0s to be safe if we have 0s that aren't in the begining
        if temp(i)==0
            numSpikes(i)=0;
        else
            c=c+1;
            numSpikes(i) = temp(i)*All(ind).out.exp.stimParams.numCells(c);
        end
    end
    
    
    All(ind).out.anal.numSpikesAddedPerCond = numSpikes;
    numSpikesEachStim{ind} = numSpikes;
    numCellsEachEns{ind} = All(ind).out.exp.stimParams.numCells;
end
numSpikesEachStim=cell2mat(numSpikesEachStim(:)');
numSpikesEachEns = numSpikesEachStim;
numSpikesEachEns(numSpikesEachStim==0)=[];

numCellsEachEns=cell2mat(numCellsEachEns(:)');
    
%% Make all dataPlots into matrixes of mean responses


clear popResponse pVisR pVisT
ensIndNumber=[];
for ind=1:numExps
    pTime =tic;
    fprintf(['Processing Experiment ' num2str(ind) '...']);

    trialsToUse = All(ind).out.exp.lowMotionTrials;

    clear respMat baseMat %Order stims,vis,cells
    for i=1:numel(unique(All(ind).out.exp.stimID))
        us = unique(All(ind).out.exp.stimID);
        s = us(i);

        for k= 1 : numel(unique(All(ind).out.exp.visID))
            vs = unique(All(ind).out.exp.visID);
            v = vs(k);

            respMat(i,k,:) = mean(All(ind).out.exp.rdData(:,...
                trialsToUse & All(ind).out.exp.stimID ==s &...
                All(ind).out.exp.visID ==v), 2) ;
            baseMat(i,k,:) = mean(All(ind).out.exp.bdata(:,...
                trialsToUse & All(ind).out.exp.stimID ==s &...
                All(ind).out.exp.visID ==v), 2) ;
        end
    end

    All(ind).out.anal.respMat = respMat;
    All(ind).out.anal.baseMat = baseMat;


    %%offtargetRisk
    stimCoM = All(ind).out.exp.stimCoM;
    numCells = size(All(ind).out.exp.zdfData,1);
    allCoM = All(ind).out.exp.allCoM;
    stimDepth = All(ind).out.exp.stimDepth;
    allDepth = All(ind).out.exp.allDepth;
    muPerPx = 800/512;

    allLoc = [allCoM*muPerPx (allDepth-1)*30];
    stimLoc = [stimCoM*muPerPx (stimDepth-1)*30];

    roisTargets = All(ind).out.exp.rois;
    holoTargets = All(ind).out.exp.holoTargets;

    thisPlaneTolerance = 15;10; %in pixels
    onePlaneTolerance = 25;20;

    radialDistToStim=zeros([size(stimCoM,1) numCells]);
    axialDistToStim = zeros([size(stimCoM,1) numCells]);
    StimDistance = zeros([size(stimCoM,1) numCells]);
    for i=1:size(stimCoM,1);
        for k=1:numCells;
            D = sqrt(sum((stimCoM(i,:)-allCoM(k,:)).^2));
            radialDistToStim(i,k)=D;
            z = stimDepth(i)-allDepth(k);
            axialDistToStim(i,k) = z;
            StimDistance(i,k) = sqrt(sum((stimLoc(i,:)-allLoc(k,:)).^2));

        end
    end

    offTargetRisk = zeros([numel(roisTargets) numCells]);
    for i=1:numel(roisTargets)
        Tg = roisTargets{i};
        try
        TgCells = holoTargets{i};
        catch;end;
        
        if numel(Tg) == 1
            temp = radialDistToStim(Tg,:)<thisPlaneTolerance & axialDistToStim(Tg,:) ==0;
            temp2 = radialDistToStim(Tg,:)<onePlaneTolerance & abs(axialDistToStim(Tg,:)) ==1;
        else
            temp = any(radialDistToStim(Tg,:)<thisPlaneTolerance & axialDistToStim(Tg,:) ==0);
            temp2 = any(radialDistToStim(Tg,:)<onePlaneTolerance & abs(axialDistToStim(Tg,:)) ==1);
        end
        offTargetRisk(i,:) = temp | temp2;
    end
    All(ind).out.anal.offTargetRisk = offTargetRisk;


    %%ROIinArtifact
    try
        yoffset = -All(ind).out.info.offsets(2);
    catch
        yoffset = 0 ;
    end

    ArtifactSizeLeft = 100;
    ArtifactSizeRight = 100;
    ROIinArtifact = allCoM(:,2)<ArtifactSizeLeft-yoffset | allCoM(:,2)>511-(ArtifactSizeRight+yoffset);
    All(ind).out.anal.ROIinArtifact = ROIinArtifact;
    pVisR = All(ind).out.anal.pVisR;
    pVisT = All(ind).out.anal.pVisT;

    %%Get Pop Responses
    %         v=1; %best bet for no vis stim.
    clear popResp popRespDist
    for v = 1:numel(vs)
        for i= 1:numel(All(ind).out.exp.stimParams.Seq)
            %             try
            holo =All(ind).out.exp.stimParams.Seq(i) ;% roi{i}{1};
            %             catch
            %                 holo =All(ind).out.exp.stimParams.roi{i};
            %             end

            if i==1;
                cellsToUse = ~ROIinArtifact' & pVisR<0.05;
            else
                cellsToUse = ~ROIinArtifact' & pVisR<0.05 & ~offTargetRisk(holo,:);
            end
            popResp(i,v) = mean(squeeze(respMat(i,v,cellsToUse) - baseMat(i,v,cellsToUse)));
            
            if i~=1
                Tg=All(ind).out.exp.rois{holo};
                dists = StimDistance(Tg,:);
                minDist = mean(dists);
                
                distBins = [0:50:500];
                for d = 1:numel(distBins)-1
                    cellsToUse = ~ROIinArtifact' &...
                        pVisR<0.05 &...
                        ~offTargetRisk(holo,:) &...
                        minDist > distBins(d) &...
                        minDist <= distBins(d+1) ;
                    popRespDist(i,v,d) = mean(squeeze(respMat(i,v,cellsToUse) - baseMat(i,v,cellsToUse)));
                end
            end
        
        
        end
    end
    
    VisCondToUse = 1; %1 is no vis
    if VisCondToUse > size(popResp,2) 
        popResponse{ind} = single(nan(size(popResp(:,1))));
        popResponseDist{ind} = single(nan(size(squeeze(popRespDist(:,1,:)))));
    else
        popResponse{ind} = popResp(:,VisCondToUse);
        popResponseDist{ind} = squeeze(popRespDist(:,VisCondToUse,:));
    end
    popResponseAll{ind} = popResp;
    
    ensIndNumber = [ensIndNumber ones(size(popResp(:,1)'))*ind];
    
    fprintf([' Took ' num2str(toc(pTime)) 's.\n'])
end

popResponse = cell2mat(popResponse(:));
popResponseEns=popResponse;
popResponseEns(numSpikesEachStim==0)=[];

ensIndNumber(numSpikesEachStim==0)=[];

noStimPopResp = popResponse(numSpikesEachStim==0);

    
%% Plot
f3 = figure(3);
clf(3)

ensemblesToUse = numSpikesEachEns > 75 & numSpikesEachEns <125;% & ensIndNumber==15; %& numCellsEachEns>10 ;
%scatter(meanOSI(ensemblesToUse),popResponseEns(ensemblesToUse),[],numCellsEachEns(ensemblesToUse),'filled')
scatter(ensOSI(ensemblesToUse),popResponseEns(ensemblesToUse),[],numCellsEachEns(ensemblesToUse),'filled')

xlabel('Ensemble OSI')
ylabel('Population Mean Response')
title('OSIs by Ensemble Size')
set(gcf(),'Name','OSIs by Ensemble Size')
cb = colorbar('Ticks', unique(numCellsEachEns(ensemblesToUse)));
cb.Label.String = 'Number of Cells in Ensemble';
r = refline(0);
r.LineStyle =':';

%% group conditions
% not proud of this but it did prevent me from re-writing a bunch of code

numCellsEachEnsBackup = numCellsEachEns;

numCellsEachEns(numCellsEachEns <=5) = 5;
numCellsEachEns(numCellsEachEns > 10) = 20;

%% fit ensembles of different sizes

f5 = figure(5);
clf(f5)
numEns = numel(unique(numCellsEachEns(ensemblesToUse)));

uniqueEns = unique(numCellsEachEns(ensemblesToUse));

for i=1:numEns
    ens2plot = find(numCellsEachEns(ensemblesToUse)==uniqueEns(i));
    p = polyfit(ensOSI(ens2plot),popResponseEns(ens2plot),1);
    f = polyval(p, ensOSI(ens2plot));
    fits(i,:) = fit;
    subplot(1,numEns,i)
    plt = plot(ensOSI(ens2plot),popResponseEns(ens2plot), '.', 'MarkerSize',12);
    hold on
    fline = plot(ensOSI(ens2plot), f, 'LineWidth', 1);
    xlabel('OSI')
    ylabel('Pop Response')
    title(['Ensembles of size ' num2str(uniqueEns(i))])
end

linkaxes

%% more simple, take the means, population response by ensemble size
clear avg err ns ens2plt
f6 = figure(6);
clf(f6)
numEns = numel(unique(numCellsEachEns(ensemblesToUse)));
uniqueEns = unique(numCellsEachEns(ensemblesToUse));

x = 1:numEns;
clear data names
for i=1:numEns
    ens2plot = find(numCellsEachEns==uniqueEns(i) & ensemblesToUse);
    data{i} = popResponseEns(ens2plot);
    names{i} = string(uniqueEns(i));
    avg(i) = mean(popResponseEns(ens2plot));
    err(i) = sem(popResponseEns(ens2plot));
    ns(i) = numel(popResponseEns(ens2plot));
end

data{end+1} = noStimPopResp;
names{end+1} = 'No Stim';

cmap=colormap(viridis(numEns));
cmap(end+1,:)=rgb('grey');
p = plotSpread(data, 'xNames', names, 'showMM', 4, 'distributionColors',cmap);
% bar(x, avg)
% hold on
% er = errorbar(x, avg, err);
% er.Color = [0 0 0];
% er.LineStyle = 'none';
% hold off
ylabel('Population Response (vis responsive)')
% xticklabels(uniqueEns)
% xticks = 1:6;
title('Mean population response to holo')
xlabel('Ensemble Size')
set(gcf(),'Name','Mean population response to holo')
% ns

ax=p{3};
set(findall(gcf(),'type','line'),'markerSize',16)
p{2}(1).Color = rgb('darkgrey');
p{2}(2).Color = rgb('darkgrey');
p{2}(1).LineWidth = 1;
p{2}(2).LineWidth = 1;

 r = refline(0);
    r.LineStyle=':';
    r.Color = rgb('grey');

pValEnselbeSize = anovan(popResponseEns(ensemblesToUse),numCellsEachEns(ensemblesToUse)','display','off')

ranksum(noStimPopResp,popResponseEns(ensemblesToUse & numCellsEachEns==5))
ranksum(noStimPopResp,popResponseEns(ensemblesToUse & numCellsEachEns==10))
ranksum(noStimPopResp,popResponseEns(ensemblesToUse & numCellsEachEns==20))


%% look at just the 10s data for each mouse


% allens2plt = popResponseEns(numCellsEachEns(ensemblesToUse))';

f7 = figure(7);
clf(f7)
k=0;
ens_ids = ensIndNumber(ensemblesToUse);
ens_sizes = numCellsEachEns(ensemblesToUse);
popResponseClip = popResponseEns(ensemblesToUse); %indexing error need to subselect first 

clear sp
for s=unique(ens_sizes)
    clear ens2plt expid exp2plt names
    k=k+1;
    hold on
    sp(k) = subplot(1,numel(unique(numCellsEachEns(ensemblesToUse))),k);
    
    expid = ens_ids(ens_sizes==s);
    ens2plt = popResponseClip(ens_sizes==s)'; %indexing error need to subselect first 

    c=0;
    for i=unique(expid)
        
        c = c+1;
        exp2plt{c} = ens2plt(expid==i);
        names{c}=strrep(All(i).out.info.mouse, '_', '.');
    end

    cmap=colormap(viridis(numel(exp2plt)));
    p=plotSpread(exp2plt,'xNames',names,'showMM',4,'distributionColors',cmap);
    ax=p{3};
    set(findall(gcf(),'type','line'),'markerSize',16)
    p{2}(1).Color = rgb('darkgrey');
    p{2}(2).Color = rgb('darkgrey');
    p{2}(1).LineWidth = 1.5;
    p{2}(2).LineWidth = 1.5;
    uistack(p{2},'bottom')
    xtickangle(45)
    title(['Ensembles of ' num2str(s)])
    
    r = refline(0);
    r.LineStyle=':';
    r.Color = rgb('grey');
    
end

linkaxes(sp(:), 'y')
ax = findobj(sp(1), 'type', 'axes');
set([ax.YLabel], 'string', 'Population Response')
set(gcf(),'Name','Population response to holo by expt and size')
sgtitle('Population response to holo by expt and size')

%% what do the mean/normalized tuning look like for low, middle, and high OSI values?
clear oriShifted lowOSIidx midOSIidx highOSIidx lowOSIcurve midOSIcurve highOSIcurve alignedOris
% set bounds for OSI
low = 0.2;
high = 0.8;

for i = 1:numel(All)
    
    lowOSIidx{i} = find(All(i).out.anal.ensembleOSI <= low);
    midOSIidx{i} = find(All(i).out.anal.ensembleOSI > low & All(i).out.anal.ensembleOSI < high);
    highOSIidx{i} = find(All(i).out.anal.ensembleOSI >= high);
    
    ensembleOriCurve = All(i).out.anal.ensembleOriCurve;
    ensemblePref = All(i).out.anal.ensemblePrefOri;
    
    % peak align to the 3rd position
    % should this be normalized somehow?
    for j = 1:size(ensembleOriCurve,1)
        oriShifted(j,:) = circshift(ensembleOriCurve(j,:),-ensemblePref(j)+3);
    end
    
    alignedOris{i} = oriShifted;
    
    lowOSIcurve{i} = oriShifted(lowOSIidx{i},:);
    midOSIcurve{i} = oriShifted(midOSIidx{i},:);
    highOSIcurve{i} = oriShifted(highOSIidx{i},:);
    
    All(i).out.anal.lowOSIcurve = oriShifted(lowOSIidx{i},:);
    All(i).out.anal.midOSIcurve = oriShifted(midOSIidx{i},:);
    All(i).out.anal.highOSIcurve = oriShifted(highOSIidx{i},:);

end

% unroll and computer errors
lowOSIcurveAll = cell2mat(lowOSIcurve(:));
err1 = std(lowOSIcurveAll)/sqrt(size(lowOSIcurveAll,1));
midOSIcurveAll = cell2mat(midOSIcurve(:));
err2 = std(midOSIcurveAll)/sqrt(size(midOSIcurveAll,1));
highOSIcurveAll = cell2mat(highOSIcurve(:));
err3 = std(highOSIcurveAll)/sqrt(size(highOSIcurveAll,1));

% plot them
f4 = figure(4);
clf(4)
hold on
errorbar(nanmean(lowOSIcurveAll,1), err1, 'linewidth',2);
errorbar(nanmean(midOSIcurveAll,1), err2, 'linewidth', 2);
errorbar(nanmean(highOSIcurveAll,1), err3, 'linewidth', 2);
hold off

title('Mean OSI Curves')
ylabel('Mean Response')
xlabel('Ori (preferred centered at 3)')
legend('Low OSIs', 'Mid OSIs', 'High OSIs')
set(gcf(),'Name','Mean OSI Curves')



%% Plot Pop Response by Distance
popDist = cell2mat(popResponseDist');


