%% Create histogram of distances between alleles in a nucleus 
% Load your constructs - listed as they are in DataStatus.xlsx
ConstructList= {'KrDist','KrProx','KrBothSep','KrDistDuplicN','KrProxDuplic','KrBoth','KrDist32C','KrProx32C','KrBothSep32C','KrBoth32C','Kr2xProx32C','KrDistDuplicN','KrDist17C','Kr2xDist32C','KrBoth17C','Kr2xDist17C','KrProx17C','Kr2xDistdBcd','Kr2xProxdHb','Kr2xProx17C','KrSEdBcd','KrInvSE','Kr2xDistLessMS2','KrSEdHb','KrEndogDist','KrEndogDist32C','KrEndogDist17C','Kr4_Chr2_Chr3','Kr2xDistChrom3','Kr3_Chr2_Chr3'} %{'KrDist','KrProx','KrBothSep', 'KrDistDuplicN', 'KrProxDuplic', 'KrBoth'};% %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');

%Colors for constructs
DistalColor=[1 64 172]./255;
DistalEmptyColor=[8 210 238] ./ 255;
Distal32CColor=[118 180 238] ./ 255;
DoubDistColor=[73 184 253] ./ 255;
ProxColor=[238 123 23]./255;
ProxEmptyColor=[251 250 50] ./255;
Proximal32CColor=[251 150 10] ./ 255;
DoubProxColor=[215 183 58] ./ 255;
DoubProxEmptyColor=[251 220 50] ./ 255;
BothSepColor=[94 250 81] ./ 255;
BothSep32CColor=[150 250 81] ./255;
BothColor=[52 119 71]./255;
BothEmptyColor=[12 250 100] ./ 255;
Both32CColor=[120 195 82] ./ 255;
DoubProx32CColor=[200 150 100] ./ 255;
PDistalColor = [129 161 214]./255;
grey = [0.5 0.5 0.5];
Grey = [0.5 0.5 0.5];

Colors(1).Color = DistalColor;
Colors(2).Color = ProxColor;
Colors(3).Color = BothSepColor;
Colors(4).Color = DoubDistColor;
Colors(5).Color = DoubProxColor;
Colors(6).Color = BothColor;

fontsize=15;
fontname='Arial';
x_width=3; y_width=2.25;
x_widthsplit=1.5; y_widthsplit=1.125;
xSize = 7; ySize = 6; xLeft = 0.5; yTop = 0.5;
FigDirect=[DropboxFolder filesep 'Figures'];

%% Look at negative co-variance as a fx of allele distance 
 
% Load your constructs - listed as they are in DataStatus.xlsx
ConstructList= {'KrDist';'KrProx';'KrBothSep';'KrDistDuplicN';'KrProxDuplic';'KrBoth';'KrDist32C';'KrProx32C';'KrBothSep32C';'KrBoth32C';'Kr2xProx32C';'KrDist17C';'Kr2xDist32C';'KrBoth17C';'Kr2xDist17C';'KrProx17C';'Kr2xDistdBcd';'Kr2xProxdHb';'Kr2xProx17C';'KrSEdBcd';'KrInvSE';'KrSEdHb';'KrEndogDist';'KrEndogDist32C';'KrEndogDist17C';'Kr2xDistLessMS2';'4F_Kr4_Kr4';'KrSEChrom3';'Kr2xDistChrom3';'Kr2xProxChrom3';'Kr4_Chr2_Chr3';'Kr3_Chr2_Chr3'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
     
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');

for cc = 1:length(ConstructList)
    Data= LoadMS2SetsCS(ConstructList{cc});
    NEmbryos = length(Data);
    EmbNucInfo=[];
    AllEmbDistVect=[];
    AllEmbNuclei =[];
    EmbMeanDist = [];
    for ee = 1:NEmbryos
        PrefixName=Data(ee).Prefix;
        CompPars=[DropboxFolder filesep PrefixName filesep 'SpotCorrelationAdj.mat'];
        load(CompPars);
        EmbNuclei=[];
        EmbDistVect=[];
        EmbAPVect =[];
        for n = 1:length(SpotDiff)
            if isfield(SpotDiff,'Frame2')
                if ~isnan([SpotDiff(n).Frame1]) %& ~isnan([SpotDiff(n).Frame2])
                    if ~isnan([SpotDiff(n).Frame2])
                        Frames_1 = [SpotDiff(n).Frame1];
                        Frames_2 = [SpotDiff(n).Frame2];
                        %Record AP info 11/25/19 
                        APInfo = SpotDiff(n).APBin;
                
                        CommonFrames = intersect(Frames_1, Frames_2); %Find the frames where both spots exist
                        IndexFrames_1= nan(1,length(CommonFrames));
                        IndexFrames_2 = nan(1,length(CommonFrames)); % set up array to fill with indices of frames where both spots exist
                        for f = 1:length(CommonFrames)
                            IndexFrames_1(f) = find(Frames_1 == CommonFrames(f));
                            IndexFrames_2(f) = find(Frames_2 == CommonFrames(f));
                        end 
                        X_vals_1 = [SpotDiff(n).xPos1(IndexFrames_1)];
                        Y_vals_1 = [SpotDiff(n).yPos1(IndexFrames_1)];
                        X_vals_2 = [SpotDiff(n).xPos2(IndexFrames_2)];
                        Y_vals_2 = [SpotDiff(n).yPos2(IndexFrames_2)]; 
                        DistanceVector = nan(1, length(CommonFrames));
                        DistanceVectorAll = nan(1, length(SpotDiff(n).nc14Frame1));
                        for f = 1:length(CommonFrames)
                            DistanceVector(f) = sqrt((X_vals_2(f) - X_vals_1(f))^2 + (Y_vals_2(f) - Y_vals_1(f))^2);
                            %Make vector with all time points, but only
                            %fill in where both spots exist
                            DistanceVectorAll(CommonFrames(f)) = sqrt((X_vals_2(f) - X_vals_1(f))^2 + (Y_vals_2(f) - Y_vals_1(f))^2);
                        end                                               
        
                    else
                        continue
                    end
                    else
                        continue
            
                    end
                    else
                        continue
                    end
            
            SpotDist(cc).Embryo(ee).NucleusInfo(n).DistanceVector = [DistanceVector];
            SpotDist(cc).Embryo(ee).NucleusInfo(n).AllPtsDistanceVector = [DistanceVectorAll];
            SpotDist(cc).Embryo(ee).MeanDist(n) = nanmean(DistanceVector);
            SpotDist(cc).Embryo(ee).APbin(n) = APInfo;
            %SpotDist(cc).Embryo(ee).Nucleus(n) = SpotDiff(n).Nucleus;
            EmbNuclei = [EmbNuclei, SpotDiff(n).Nucleus];
            %EmbMeanDist = [EmbMeanDist, nanmean(DistanceVector)];
            EmbNucInfo = [EmbNucInfo, SpotDiff(n).Nucleus];
            EmbDistVect = [EmbDistVect, nanmean(DistanceVector)];
            EmbAPVect = [EmbAPVect, APInfo];
        end
        SpotDist(cc).Embryo(ee).AllMeanDist = [EmbDistVect] %Vector of all mean allele distances of embryo
        SpotDist(cc).Embryo(ee).Nuclei = EmbNuclei; %Vector of all nuclei #'s in same order 
        SpotDist(cc).Embryo(ee).APBins =EmbAPVect;
        SpotDist(cc).Embryo(ee).Construct = ConstructList{cc};
        
    end
end
save([DropboxFolder filesep 'Constructs' filesep 'SpotDistance'],'SpotDist');

%% Connect negative covariance with allele distance of each nucleus 
% Load your constructs - listed as they are in DataStatus.xlsx
ConstructList= {'KrDist';'KrProx';'KrBothSep';'KrDistDuplicN';'KrProxDuplic';'KrBoth';'KrDist32C';'KrProx32C';'KrBothSep32C';'KrBoth32C';'Kr2xProx32C';'KrDist17C';'Kr2xDist32C';'KrBoth17C';'Kr2xDist17C';'KrProx17C';'Kr2xDistdBcd';'Kr2xProxdHb';'Kr2xProx17C';'KrSEdBcd';'KrInvSE';'KrSEdHb';'KrEndogDist';'KrEndogDist32C';'KrEndogDist17C';'Kr2xDistLessMS2';'4F_Kr4_Kr4';'KrSEChrom3';'Kr2xDistChrom3';'Kr2xProxChrom3';'Kr4_Chr2_Chr3';'Kr3_Chr2_Chr3'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
     
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');

load([DropboxFolder filesep 'Constructs' filesep 'TotalNoiseData.mat']);
for cc = 1:length(ConstructList)
    Data= LoadMS2SetsCS(ConstructList{cc});
    NEmbryos = length(Data);
    for pp = 1:length(WholeNoise)
        if strcmp(WholeNoise(pp).ConstructName,ConstructList{cc}) %see which row in WholeNoise structure corresponds to this construct
            Noiseset = pp;
            break
        end
    end
    EmbNegCovar =[];
    EmbAlleleDist =[];
    EmbAPBinInfo=[];
    for ee = 1:NEmbryos

        NegCovar=[];
        AlleleDist=[];
        APBinInfo =[];
        for nn = 1:length(SpotDist(cc).Embryo(ee).Nuclei)
            temp = find([WholeNoise(pp).Embryo(ee).Nucleus] == SpotDist(cc).Embryo(ee).Nuclei(nn));
            NegCovar = [NegCovar, WholeNoise(pp).Embryo(ee).CoVarNoise(temp)]; %Vector of values for a single embryo
            AlleleDist = [AlleleDist, SpotDist(cc).Embryo(ee).AllMeanDist(nn)];
            %Add AP bin info 11/25/19
            APBinInfo = [APBinInfo, SpotDist(cc).Embryo(ee).APbin(nn)];
        end
        SpotDist(cc).Embryo(ee).NegCovar = NegCovar;
        SpotDist(cc).Embryo(ee).AlleleDist = AlleleDist;
        EmbNegCovar = [EmbNegCovar, NegCovar];
        EmbAlleleDist = [EmbAlleleDist, AlleleDist];
        EmbAPBinInfo = [EmbAPBinInfo, APBinInfo];
        clear NegCovar AlleleDist APBinInfo
    end
    SpotDist(cc).NegCovar = EmbNegCovar;
    SpotDist(cc).AlleleDist = EmbAlleleDist;
    SpotDist(cc).APBin = EmbAPBinInfo;
    SpotDist(cc).Construct = ConstructList{cc};
   
end

save([DropboxFolder filesep 'Constructs' filesep 'SpotDistance'],'SpotDist');

%% Plotting
% Plot 1/2 violin plots of reporter distance broken into nuclei w vs w/o
% negative covariance for the 3 main reporter constructs 
figure
PixelSizeZoom = 0.2436; % Value from bfMetaData.mat -- this will be microscope specific!
PlaceCounter=1;
for cc = [4:6] %,29:31]%length(ConstructList)]%length(ConstructList)
    
    
    
    Data= LoadMS2SetsCS(ConstructList{cc});
    NEmbryos = length(Data);
    CompareMatrix = []; %reset for each construct looking at 
    
    NegIndex = find([SpotDist(cc).NegCovar] < 0); %indices of nuclei w negative covariance
    PosIndex = find([SpotDist(cc).NegCovar] >= 0); %indices of nuclei w covariance >= 0
    
    TotalDataPts_Neg=[SpotDist(cc).AlleleDist(NegIndex)]';
    TotalDataPts_Pos = [SpotDist(cc).AlleleDist(PosIndex)]';
    
    PointsUsed_Neg=length(TotalDataPts_Neg);
    PointsUsed_Pos = length(TotalDataPts_Pos);
    %plot in actual distance units by multiplying by size of individual pixel
    %in um
    % Try doing split violin plot 7/13/20
    hold on 
    violinPlot([TotalDataPts_Pos].*PixelSizeZoom,'xValues',[PlaceCounter],'color',Colors(cc).Color,'histOpt',1.1,'divFactor',1.5,'widthDiv',[2,1],'histOri','left','showMM',0);
    violinPlot([TotalDataPts_Neg].*PixelSizeZoom,'xValues',[PlaceCounter],'color',[Colors(cc).Color .* 0.5],'histOpt',1.1,'divFactor',1.5,'widthDiv',[2,2],'histOri','right','showMM',0);
    
    hold on
    y_1=nanmedian(([SpotDist(cc).AlleleDist(PosIndex)]).*PixelSizeZoom);
    q=line([PlaceCounter-0.5,PlaceCounter],[y_1,y_1]); q.Color='k'; q.LineWidth=2;
    
    
    y_2=nanmedian(([SpotDist(cc).AlleleDist(NegIndex)]).*PixelSizeZoom);
    q=line([PlaceCounter,PlaceCounter+0.5],[y_2,y_2]); q.Color='k'; q.LineWidth=2;

    %Figure out ylimits based on 75th percentile 
    UpLim(cc) = prctile(TotalDataPts_Pos,99);
    %Make sure negative cov nuclei don't have higher 75th percentile in
    %case negative-closer isn't true for some constructs 
    if (prctile(TotalDataPts_Neg,70)) > UpLim(cc)
        UpLim(cc) = prctile(TotalDataPts_Neg,99);
    end
    %ylim([0 (nanmax(UpLim).*PixelSizeZoom)]);
    ylim([0 6]);
    ylabel('projected reporter distance (um)');
    xticks([]);
    set(gca, 'YColor','k');
    set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
    set(gcf,'PaperUnits','inches');
    set(gcf,'PaperPosition',[xLeft yTop 7 6]);
    %Kolmogorov-Smirnov test
    [h(cc),p(cc),ks2stat(cc)]=kstest2([SpotDist(cc).AlleleDist(NegIndex)],[SpotDist(cc).AlleleDist(PosIndex)]);
    
    %should be fine to do KS test bc only comparing two groups, need KW test w multiple comparision corrections when comparing mult groups 
    
    PlaceCounter =PlaceCounter+1.5;
    %title({num2str(p(cc)) 'p-val'});
    %xticklabels({'positive co-variance', 'negative co-variance'})
    
end

%% Violin plot compare all avg distances btwn E's on Chr2
PixelSizeZoom = 0.2436; % Value from bfMetaData.mat 
FigDirect=[DropboxFolder filesep 'Figures'];

%Compare on Chr2
figure
PlaceCounter =1;
for cc = [4:6]
 
    
    Data= LoadMS2SetsCS(ConstructList{cc});
    NEmbryos = length(Data);
    CompareMatrix = []; %reset for each construct looking at 
    
    TotalDataPts_1=[SpotDist(cc).AlleleDist]';
    
    PointsUsed_1=length(TotalDataPts_1);
    
    % 5/8/20 add mult by pixelsize to have in um units 
    v=violinPlot([TotalDataPts_1].*PixelSizeZoom,'xValues',PlaceCounter,'color',Grey,'histOpt',1.1,'divFactor',1.5,'showMM',0);
    hold on
    %Overlay the total noise points of the individual nuclei
    scatter((PlaceCounter.*ones(PointsUsed_1,1)),(([SpotDist(cc).AlleleDist].*PixelSizeZoom)),5,'k','MarkerFaceColor',Colors(cc).Color,'jitter','on','jitterAmount',0.3);
    y=nanmedian((SpotDist(cc).AlleleDist(:)).*PixelSizeZoom);
    q=line([PlaceCounter-0.5,PlaceCounter+0.5],[y,y]); q.Color='k'; q.LineWidth=2;
    
    PlaceCounter=PlaceCounter+1;
end   
ylim([-0.5 6]);
title('chr2');    %xlabel('allele distance (AU)');
ylabel('projected reporter distance (um)');
xticks([1,2,3]);
saveas(gcf, [FigDirect filesep 'CompAlleleDist2xEs_Ch2','.pdf'],'pdf');

%Kruskal-wallis test w Bonferroni correction
 %Initialize matrix to hold values and group name
 AlleleDistKS = nan((length(SpotDist(4).AlleleDist)+length(SpotDist(5).AlleleDist)+length(SpotDist(6).AlleleDist)),1);
 AlleleGroups = nan((length(SpotDist(4).AlleleDist)+length(SpotDist(5).AlleleDist)+length(SpotDist(6).AlleleDist)),1);
 % populate with values 
 DLength = length(SpotDist(4).AlleleDist); PLength = length(SpotDist(5).AlleleDist); SELength = length(SpotDist(6).AlleleDist);
 AlleleDistKS([1:DLength]) = SpotDist(4).AlleleDist(:);
 AlleleDistKS([DLength+1:DLength+PLength]) = SpotDist(5).AlleleDist(:);
 AlleleDistKS([DLength+PLength+2:DLength+PLength+SELength+1]) = SpotDist(6).AlleleDist(:);
 
 AlleleGroups([1:DLength]) = 4;
 AlleleGroups([DLength+1:DLength+PLength]) = 5;
 AlleleGroups([DLength+PLength+2:DLength+PLength+SELength+1]) = 6;
[p,tbl,stats]=kruskalwallis(AlleleDistKS,AlleleGroups);
[c,m]=multcompare(stats,'CType','bonferroni'); 
 
%% Plot distribution and compare of total distribution of allele distances btwn reporter constructs 

EDist_SE=[];
EDist_DD=[];
EDist_DP=[];
EDist_Group =[];

for ee = 1:length(SpotDist(4).Embryo)
    EDist_DD = [EDist_DD, [SpotDist(4).Embryo(ee).NucleusInfo(:).DistanceVector]];
    EDist_Group = [EDist_Group, ones(1,length([SpotDist(4).Embryo(ee).NucleusInfo(:).DistanceVector]))];
end

for ee = 1:length(SpotDist(5).Embryo)
    EDist_DP = [EDist_DP, [SpotDist(5).Embryo(ee).NucleusInfo(:).DistanceVector]];
    EDist_Group = [EDist_Group,ones(1,length([SpotDist(5).Embryo(ee).NucleusInfo(:).DistanceVector])).*2];
end

for ee = 1:length(SpotDist(6).Embryo)
    EDist_SE = [EDist_SE, [SpotDist(6).Embryo(ee).NucleusInfo(:).DistanceVector]];
    EDist_Group = [EDist_Group,ones(1,length([SpotDist(6).Embryo(ee).NucleusInfo(:).DistanceVector])).*3];
end

EDist_All3 = [EDist_DD, EDist_DP, EDist_SE];

%KS test w Bonferoni mult comparison correction
[p,tbl,stats]=kruskalwallis(EDist_All3,EDist_Group);
[c,m]=multcompare(stats,'CType','bonferroni'); 

% plotting
%2xDistal
figure
PlaceCounter=1;
PointsUsed = length(EDist_DD); 
v=violinPlot([EDist_DD]'.*PixelSizeZoom,'xValues',PlaceCounter,'color',Colors(4).Color,'histOpt',1.1,'divFactor',1.5,'showMM',0);
hold on 
y=nanmedian([EDist_DD].*PixelSizeZoom);
q=line([PlaceCounter-0.5,PlaceCounter+0.5],[y,y]); q.Color='k'; q.LineWidth=2;

%2xProximal
PlaceCounter = 2; 
PointsUsed = length(EDist_DP); 
v=violinPlot([EDist_DP]'.*PixelSizeZoom,'xValues',PlaceCounter,'color',Colors(5).Color,'histOpt',1.1,'divFactor',1.5,'showMM',0);
y=nanmedian([EDist_DP].*PixelSizeZoom);
q=line([PlaceCounter-0.5,PlaceCounter+0.5],[y,y]); q.Color='k'; q.LineWidth=2;

%shadow pair
PlaceCounter = 3; 
PointsUsed = length(EDist_SE); 
v=violinPlot([EDist_SE]'.*PixelSizeZoom,'xValues',PlaceCounter,'color',Colors(6).Color,'histOpt',1.1,'divFactor',1.5,'showMM',0);
y=nanmedian([EDist_SE].*PixelSizeZoom);
q=line([PlaceCounter-0.5,PlaceCounter+0.5],[y,y]); q.Color='k'; q.LineWidth=2;

% indicate model cut off of hub size
q=line([0,4],[0.512,0.512]); q.Color=[0.5 0.5 0.5]; q.LineWidth=2, q.LineStyle='--';

ylabel('projected reporter distance (um)');
ylim([-0.5 7]);
xticks([]);
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);

saveas(gcf, [FigDirect filesep 'AllReporterDistDistribution','.pdf'],'pdf');
