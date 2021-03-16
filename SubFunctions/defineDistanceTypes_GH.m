function [All, outVars] = defineDistanceTypes_GH(All, outVars)
numExps = numel(All);
ensMeaD=[];ensGeoD=[];ensMaxD=[];ensMinD=[]; ensCenterD=[]; ensCenterV=[];

for ind = 1:numExps
    pTime =tic;
    fprintf(['Processing Experiment ' num2str(ind) '...']);
    
    stimCoM = All(ind).out.exp.stimCoM;
    numCells = size(All(ind).out.exp.zdfData,1);
    allCoM = All(ind).out.exp.allCoM;
    stimDepth = All(ind).out.exp.stimDepth;
    allDepth = All(ind).out.exp.allDepth;
    muPerPx = 800/512;
    
    allLoc = [allCoM*muPerPx (allDepth-1)*30];
    stimLoc = [stimCoM*muPerPx (stimDepth-1)*30];

    StimDistance = zeros([size(stimCoM,1) numCells]);
    for i=1:size(stimCoM,1);
        for k=1:numCells;
            StimDistance(i,k) = sqrt(sum((stimLoc(i,:)-allLoc(k,:)).^2));
        end
    end
    All(ind).out.anal.StimDistance = StimDistance;
    
    allDistance = zeros([numCells numCells]);
    for i=1:numCells;
        for k=1:numCells;
            allDistance(i,k) = sqrt(sum((allLoc(i,:)-allLoc(k,:)).^2));
        end
    end
    All(ind).out.anal.allDistance = allDistance;
    
    roiDistance = zeros([size(stimCoM,1) size(stimCoM,1)]);
    for i=1:size(stimCoM,1);
        for k=1:size(stimCoM,1);
            roiDistance(i,k) = sqrt(sum((stimLoc(i,:)-stimLoc(k,:)).^2));
        end
    end
    All(ind).out.anal.roiDistance = roiDistance;
    
    fprintf([' Took ' num2str(toc(pTime)) 's.\n'])
% end

%% Determine the Ensemble CoCorrelation

% for ind = 1:numExps
    numStims = numel(All(ind).out.exp.stimParams.Seq);
    
    clear ensembleMeanDist ensembleGeoMeanDist ensembleMaxDist ensembleMinDist ensembleCenterDist ensembleCenterVar
    %     for i =1:numel(All(ind).out.exp.holoTargets)
    %         ht = All(ind).out.exp.holoTargets{i};
    %         ht(isnan(ht))=[];
    c=0;
    for i= 1:numStims
        holo = All(ind).out.exp.stimParams.roi{i}; % Better Identifying ensemble
        if holo>0
            c=c+1;
            ht = All(ind).out.exp.holoTargets{holo};
            ht(isnan(ht))=[];
            
            rt = All(ind).out.exp.rois{holo};
            
            if size(stimLoc(All(ind).out.exp.rois{holo},:),1)==1
                center_of_mass = NaN;
            else
                center_of_mass = mean(stimLoc(All(ind).out.exp.rois{holo},:));
            end
            
%             holo
%             All(ind).out.exp.rois{holo}
%             size(allLoc)
%             allLoc(All(ind).out.exp.rois{holo},:)

            ensembleCenterVar(c) = sqrt(sum(diag(cov((center_of_mass - stimLoc(All(ind).out.exp.rois{holo},:))))));            
            ensembleCenterDist(c) = mean(sqrt(sum((center_of_mass - stimLoc(All(ind).out.exp.rois{holo},:)).^2,2)));

            distToUse = All(ind).out.anal.roiDistance;
            distMat = distToUse(rt,rt);
            distMat(logical(eye(numel(rt))))=nan;
            
            ensembleMeanDist(c)     = nanmean(distMat(:));
            tempDist = distMat(~isnan(distMat));
            try
                ensembleGeoMeanDist(c)  = geo_mean(tempDist(:));
            catch
                ensembleGeoMeanDist(c)  = geomean(tempDist(:));
            end
            ensembleMaxDist(c)      = max(distMat(:));
            ensembleMinDist(c)      = min(distMat(:));
            
        end
    end
 
    All(ind).out.anal.ensembleMeanDist      = ensembleMeanDist;
    All(ind).out.anal.ensembleGeoMeanDist   = ensembleGeoMeanDist;
    All(ind).out.anal.ensembleMaxDist       = ensembleMaxDist;
    All(ind).out.anal.ensembleMinDist       = ensembleMinDist;
    All(ind).out.anal.ensembleCenterDist    = ensembleCenterDist;
    All(ind).out.anal.ensembleCenterVar    = ensembleCenterVar;
    
    ensMeaD = cat(2,ensMeaD,ensembleMeanDist);
    ensGeoD = cat(2,ensGeoD,ensembleGeoMeanDist);
    ensMaxD = cat(2,ensMaxD,ensembleMaxDist);
    ensMinD = cat(2,ensMinD,ensembleMinDist);
    ensCenterD = cat(2,ensCenterD,ensembleCenterDist);
    ensCenterV = cat(2,ensCenterV,ensembleCenterVar);
end

outVars.ensMeaD=ensMeaD;
outVars.ensGeoD=ensGeoD;
outVars.ensMaxD=ensMaxD;
outVars.ensMinD=ensMinD;
outVars.ensCenterD = ensCenterD;
outVars.ensCenterV = ensCenterV;


