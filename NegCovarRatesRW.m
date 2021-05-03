%% calculate coefficient of variation for each construct 
%load constructs from MovieDatabase.csv
ConstructList= {'KrDist';'KrProx';'KrBothSep';'KrDistDuplicN';'KrProxDuplic';'KrBoth';'KrDist32C';'KrProx32C';'KrBothSep32C';'KrBoth32C';'Kr2xProx32C';'KrDist17C';'Kr2xDist32C';'KrBoth17C';'Kr2xDist17C';'KrProx17C';'Kr2xDistdBcd';'Kr2xProxdHb';'Kr2xProx17C';'KrSEdBcd';'KrInvSE';'KrSEdHb';'KrEndogDist';'KrEndogDist32C';'KrEndogDist17C';'Kr2xDistLessMS2';'4F_Kr4_Kr4';'KrSEChrom3';'Kr2xDistChrom3';'Kr2xProxChrom3';'Kr4_Chr2_Chr3';'Kr3_Chr2_Chr3'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location

%Ask if only want nc14 info   for right now just going to use nc14
ncUse=input('Want to only use nc14?','s');
NumberNuclei=nan(20,41,[length(ConstructList)]);
AvoidNegs=input('Exclude negative co-variance?','s'); %Testing effect of removing info of nuclei w negative co-variance values 1/18/19

% go through each embryo of each construct
for cc=1:length(ConstructList)
     Data= LoadMS2SetsCS(ConstructList{cc});
    Datalength(cc)=length(Data);
    NEmbryos = length(Data);
    APbinID=[Data(1).APbinID];
    Label = ConstructList(cc);
    APTable=[];
    AllCorrTotmRNA=[];
    ConTimeAvg=[];
    AllNucTimeTable=[];
    AllIntraNoise=[];
    AllCoVarNoise=[];
    AllCorrSpots=[];
    AllTotalNoise=[];
    ConMeanTable=[];
    AllTotalmRNA=[];
    AllAvgProds=[];
    AllVarProds=[];
    AllNucInfo=[]; %Added 3/9/20 
    for ee=1:NEmbryos
        
        PrefixName=Data(ee).Prefix;
        if ncUse=='y'
        Filename=[DropboxFolder filesep PrefixName filesep 'SpotCorrelationAdj.mat'];
        else
            Filename=[DropboxFolder filesep PrefixName filesep 'SpotCorrelation.mat'];
        end
        load(Filename);
       
        APstuff=[SpotDiff.APBin];
        %Initialize so the same location in each parameter corresponds to
        %the same transcription spot 
        BothTotmRNA=nan(60,41);
        IntraNoise=nan(60, 41);
        CoVarNoise=nan(60,41);
        TotalNoiseVal=nan(60,41);  %bumped up to 60 1/14/19
        CorrSpots=nan(60,41);
        NucleusInfo = nan(60,41); %added to keep track of which nucleus data corresponds to 8/28/19 - using in other scripts

        for aa=1:length(APbinID)
            APsubset=[];
            APsubset=SpotDiff(APstuff==APbinID(aa));
            Spotstuff=[APsubset.SpotOne];
            APSpotsubset=APsubset(~isempty(Spotstuff));
            APNucs(aa)=length(APsubset);
        end
        MostAPNucs=max(APNucs);
        MostAPNucs=2*MostAPNucs;
        TimeTable=nan(MostAPNucs,41);
        MeanzTable=nan(MostAPNucs,41);
        TotalmRNA=nan(MostAPNucs,41);
        AvgProdAlleleOne=nan(MostAPNucs,41); %added 2/4/19 for mean-var relationships
        VarProdAlleleOne=nan(MostAPNucs,41);
        for aa=1:length(APbinID)
            APsubset=[];
%           
            APsubset=SpotDiff(APstuff==APbinID(aa));
            if ~isempty(APsubset)
               
            for bb=1:length(APsubset)
                AvgSpotOne=[];
                AvgSpotTwo=[];
                DiffVal=[];
                SpotCorr=[];
                SquaredSum=[];
                MultVal=[];
                NumberNuclei(ee,aa,cc)=length(APsubset);
                TimeLength=[];
                %1/22/19 changed to smoothed traces
                TimeLength=sum(~isnan([APsubset(bb).SmoothSpotOne]));   %Only want to look at frames where nucleus exisits (i.e 0 or number values)
                AvgFluo=nanmean([APsubset(bb).SmoothSpotOne]);
                VarFluo=nanstd([APsubset(bb).SmoothSpotOne]);
                
            TimeTable(bb,aa)=(VarFluo/AvgFluo); %CV values 
            MeanzTable(bb,aa)=AvgFluo; %save for mean-noise relationships
            
            NucleusInfo(bb,aa) = APsubset(bb).Nucleus; %record which nucleus data corresponds to 8/28/19

            %save total mRNA produced by each spot to allow to compare
            %noise-expression later 
            if ~isempty(APsubset(bb).TotalmRNAOne)
            TotalmRNA(bb,aa)=APsubset(bb).TotalmRNAOne;
            else
                TotalmRNA(bb,aa)=nan;
            end
            for ss=1:length(APsubset(bb).SmoothSpotOne)
                if ~isfield(APsubset,'SpotTwo')
                    break
                else
                if length([APsubset(bb).SmoothSpotOne]) ~= length([APsubset(bb).SmoothSpotTwo])
                    break
                else
                    %Perform calculations needed for
                    %inter-allele/covariance formulas 
                    DiffVal(ss)=((APsubset(bb).SmoothSpotOne(ss) - APsubset(bb).SmoothSpotTwo(ss))^2);
                    SquaredSum(ss)=(((APsubset(bb).SmoothSpotOne(ss)^2))+(APsubset(bb).SmoothSpotTwo(ss))^2);
                    MultVal(ss)=(APsubset(bb).SmoothSpotOne(ss) * (APsubset(bb).SmoothSpotTwo(ss)));
                    %Also save the individual spot values for correlation calculations
                    SpotCorr(ss,1)=APsubset(bb).SmoothSpotOne(ss); 
                    SpotCorr(ss,2)=APsubset(bb).SmoothSpotTwo(ss);
                end
                end
            end
            %Getting the mean values needed for inter-allele/covariance
            AvgDiffVal=nanmean(DiffVal);
            AvgSqrSum=nanmean(SquaredSum);
            AvgSpotOne=nanmean(APsubset(bb).SmoothSpotOne);
            VarSpotOne=var(APsubset(bb).SmoothSpotOne);
            AvgMultVal=nanmean(MultVal);
            %Calculate correlation of activity of the two alleles
            IndSpotCorr=corrcoef(SpotCorr,'Rows','complete');
            if ~isfield(APsubset,'SpotTwo')
                break
            else
                if isempty(APsubset(bb).TotalmRNATwo) | (isempty(APsubset(bb).TotalmRNAOne))
                    BothTotmRNA(bb,aa)=nan;
                    NucleusInfo(bb,aa) = APsubset(bb).Nucleus;
                else
                    BothTotmRNA(bb,aa)=(APsubset(bb).TotalmRNAOne+APsubset(bb).TotalmRNATwo); %Total mRNA produced in a nucleus by both alleles
                    NucleusInfo(bb,aa) = APsubset(bb).Nucleus;
                end
                AvgSpotTwo=nanmean(APsubset(bb).SmoothSpotTwo);
                %Calculate inter-allele noise and covariance (squared
                %values)
                 IntraNoise(bb,aa)=(AvgDiffVal/((2*(AvgSpotOne*AvgSpotTwo))));
                 CoVarNoise(bb,aa)=(((AvgMultVal) - ((AvgSpotOne)*(AvgSpotTwo)))/((AvgSpotOne) * (AvgSpotTwo))); 
                 TotalNoiseVal(bb,aa)=((AvgSqrSum-(2*(AvgSpotOne)*(AvgSpotTwo)))/(2*(AvgSpotOne)*(AvgSpotTwo))); %(<m^2 + p^2> - 2<m><p>)/(2<m><p>)
                 %Save mean and variance of the alleles
                 AvgProdAlleleOne(bb,aa)=AvgSpotOne;
                 VarProdAlleleOne(bb,aa)=VarSpotOne;
                 
                 %If chose to avoid negative co-variance this is where
                 %that happens 
                 if (AvoidNegs=='y') & (CoVarNoise(bb,aa) < 0)
                     CoVarNoise(bb,aa)=nan;
                     IntraNoise(bb,aa)=nan;
                     TotalNoiseVal(bb,aa)=nan;
                 end

              end
            if length(IndSpotCorr) >1
                CorrSpots(bb,aa)=IndSpotCorr(1,2);
            else
                CorrSpots(bb,aa)=nan;
            end
                
            end
            for bb=1:length(APsubset) %Save the mean and variance of the second alleles
                if ~isfield(APsubset,'SpotTwo')
                    break
                elseif ~isempty(APsubset(bb).SmoothSpotTwo)
                    TimeLength2=[];
                TimeLength2=(sum(~isnan(APsubset(bb).SmoothSpotTwo)));
                AvgFluo2=nanmean([APsubset(bb).SmoothSpotTwo]);
                VarFluo2=nanstd([APsubset(bb).SmoothSpotTwo]);
                
                TimeTable(bb+(length(APsubset)),aa)=(VarFluo2/AvgFluo2); %CV values
                if ~isempty(APsubset(bb).TotalmRNATwo);
                TotalmRNA(bb+(length(APsubset)),aa)=APsubset(bb).TotalmRNATwo
                else
                    TotalmRNA(bb+(length(APsubset)),aa)=nan;
                end
                MeanzTable(bb+(length(APsubset)),aa)=AvgFluo2;
                AvgProdAlleleOne(bb+(length(APsubset)),aa)=AvgFluo2;
                VarProdAlleleOne(bb+(length(APsubset)),aa)=VarFluo2;
                end
            end
             TotalmRNA(TotalmRNA==0)=nan;
                
                %TotalNoiseVal(TotalNoiseVal==0)=nan; removed 6/14/19
                
            end
        end
        
        CoVarNoise=real(CoVarNoise); %values were including tiny imaginary portions so get rid of those
        WholeNoise(cc).Embryo(ee).NoiseVals=IntraNoise;
        WholeNoise(cc).Embryo(ee).CoVarNoise=CoVarNoise;
        WholeNoise(cc).Embryo(ee).TotalNoiseCalc=TotalNoiseVal;
        WholeNoise(cc).Embryo(ee).SpotCorr=CorrSpots;
        WholeNoise(cc).Embryo(ee).TimeTable=TimeTable;
        WholeNoise(cc).Embryo(ee).TotalmRNA=TotalmRNA;
        WholeNoise(cc).Embryo(ee).MeanFluo=MeanzTable;
        WholeNoise(cc).Embryo(ee).EmbryoAvg=nanmean(TimeTable);
        WholeNoise(cc).Embryo(ee).Nucleus = NucleusInfo; %save nucleus info for comparison to allele distance
        
        AllCorrTotmRNA=[AllCorrTotmRNA; BothTotmRNA];
        AllNucTimeTable=[AllNucTimeTable;TimeTable];
        ConTimeAvg=[ConTimeAvg;(nanmean(TimeTable))];
        AllIntraNoise=[AllIntraNoise; IntraNoise];
        AllCoVarNoise=[AllCoVarNoise; CoVarNoise];
        AllTotalNoise=[AllTotalNoise;TotalNoiseVal];
        AllCorrSpots=[AllCorrSpots; CorrSpots];
        ConMeanTable=[ConMeanTable; MeanzTable];
        AllTotalmRNA=[AllTotalmRNA; TotalmRNA];
        AllAvgProds=[AllAvgProds; AvgProdAlleleOne]; %added 2/4/19 for mean-var relationships
        AllVarProds=[AllVarProds; VarProdAlleleOne];
        AllNucInfo =[AllNucInfo; NucleusInfo];
    end
    WholeNoise(cc).ConstructAvgNoise=nanmean(ConTimeAvg);
    WholeNoise(cc).CorrTotNoise=AllCorrTotmRNA;
    WholeNoise(cc).AllNucsCV=AllNucTimeTable;
    WholeNoise(cc).AvgCVAllNucs=nanmean([WholeNoise(cc).AllNucsCV]);
    WholeNoise(cc).SDAllNucs=nanstd([WholeNoise(cc).AllNucsCV]);
    WholeNoise(cc).SEAllNucs=(WholeNoise(cc).SDAllNucs)/(sqrt(length(WholeNoise(cc).AllNucsCV)));
    WholeNoise(cc).AllIntraNoise=AllIntraNoise;
    WholeNoise(cc).SDAllIntraNoise=nanstd(AllIntraNoise);
    WholeNoise(cc).SEAllIntraNoise=(nanstd(AllIntraNoise)/(sqrt(length(AllIntraNoise(~isnan(AllIntraNoise))))));
    WholeNoise(cc).Conf95InterNoise=(WholeNoise(cc).SEAllIntraNoise).*1.96;
    WholeNoise(cc).AllCoVarNoise=AllCoVarNoise;
    WholeNoise(cc).SDAllCoVarNoise=nanstd(AllCoVarNoise);
    WholeNoise(cc).SEAllCoVarNoise=(WholeNoise(cc).SDAllCoVarNoise)/(sqrt(length(AllCoVarNoise(~isnan(AllCoVarNoise)))));
    WholeNoise(cc).Conf95CoVar=(WholeNoise(cc).SEAllCoVarNoise).*1.96;
    WholeNoise(cc).TotalNoise=AllTotalNoise;
    WholeNoise(cc).SDTotalNoise=nanstd(AllTotalNoise);
    WholeNoise(cc).SETotalNoise=(WholeNoise(cc).SDTotalNoise)/(sqrt(length(AllTotalNoise(~isnan(AllTotalNoise)))));
    WholeNoise(cc).Conf95TotalNoise=(WholeNoise(cc).SETotalNoise).*1.96;
    WholeNoise(cc).AllCorrSpots=AllCorrSpots;
    WholeNoise(cc).SumSystNoise=[AllIntraNoise + AllCoVarNoise];
    WholeNoise(cc).AllAvgFluoVal=[ConMeanTable];
    WholeNoise(cc).AllTotalmRNA=[AllTotalmRNA];
    WholeNoise(cc).AllAvgProds=[AllAvgProds];  %added 2/4/19 RW for mean-var relationships
    WholeNoise(cc).AllVarProds=[AllVarProds];
    WholeNoise(cc).ConstructName=ConstructList{cc}; %Save associated construct name 6/19/19 RW
    WholeNoise(cc).NucInfo = AllNucInfo; %Added 3/9/20 
end

% Get rid of places where only have one data point for a whole AP bin
for cc=1:length(ConstructList)
    for aa=1:length(APbinID)
        if sum(~isnan(WholeNoise(cc).AllNucsCV(:,aa))) ==1
            WholeNoise(cc).AllNucsCV([find(~isnan(WholeNoise(cc).AllNucsCV(:,aa)))],aa)=nan;
            WholeNoise(cc).AvgCVAllNucs(aa)=nan;
            WholeNoise(cc).TotalNoise(:,aa)=nan;
            WholeNoise(cc).AllIntraNoise(:,aa)=nan;
            WholeNoise(cc).AllCoVarNoise(:,aa)=nan;
        end
        NumberNuclei(cc,aa)=sum(~isnan(WholeNoise(cc).AllNucsCV(:,aa)));
    end
end

%% Plotting info 
DistalColor=[1 64 172]./255;
Distal32CColor=[118 180 238] ./ 255;
DistalEmptyColor=[8 210 238] ./ 255; 
DoubDistColor=[73 184 253] ./ 255;
ProxColor=[238 123 23]./255;
ProxEmptyColor=[251 250 50] ./255;
Proximal32CColor=[251 150 10] ./ 255;
DoubProxColor=[215 183 58] ./ 255;
DoubProxEmptyColor=[251 220 50] ./ 255;
BothSepColor=[94 250 81] ./ 255;
BothSep32CColor=[150 250 81] ./255;
BothColor=[52 119 71]./255;
Both32CColor=[149 188 114] ./ 255;
BothEmptyColor=[12 250 100] ./ 255;
DoubProx32CColor=[200 150 100] ./ 255;
DoubDist4FColor = 'k'; 

Colors(1).Color=DistalColor;
Colors(2).Color=ProxColor;
Colors(3).Color=BothSepColor;
Colors(4).Color=DoubDistColor;
Colors(5).Color=DoubProxColor;
Colors(6).Color=BothColor;
Colors(7).Color=DistalColor;
Colors(8).Color=ProxColor;
Colors(9).Color=BothSepColor;
Colors(10).Color=BothColor;
Colors(11).Color=DoubProxColor;
Colors(12).Color=DistalColor;
Colors(13).Color=DoubDistColor;
Colors(14).Color=BothColor;
Colors(15).Color=DoubDistColor;
Colors(16).Color=ProxColor;
Colors(17).Color=DoubDistColor;
Colors(18).Color=DoubProxColor;
Colors(19).Color=DoubProxColor;
Colors(20).Color=BothColor;
Colors(21).Color=BothColor;
Colors(22).Color=BothColor;
Colors(23).Color=DistalColor;
Colors(24).Color=DistalColor;
Colors(25).Color=DistalColor;
Colors(26).Color = DoubDistColor;
Colors(27).Color = DoubDist4FColor;
Colors(28).Color = BothColor;
Colors(29).Color = DoubDistColor;
Colors(30).Color = DoubProxColor;
Colors(31).Color = DoubDistColor;
Colors(32).Color = BothColor;

% font/sizing info

fontsize=10;
fontname='Arial';
x_width=3; y_width=2.25;
x_widthsplit=1.5; y_widthsplit=1.125;
xSize = 7; ySize = 6; xLeft = 0.5; yTop = 0.5;

EggLength=APbinID .* 100;
% Set where to save figures to 
FigDirect=[DropboxFolder filesep 'Figures'];
% Save structure needed to create figures
if AvoidNegs == 'n'
save([DropboxFolder filesep 'Constructs' filesep 'TotalNoiseData'],'WholeNoise');
end

%% Make reference data structure

for cc=1:length(ConstructList)
    APAvgs=[nanmean(WholeNoise(cc).AllTotalmRNA)];
    APbinToUse=find(nanmean(WholeNoise(cc).AllTotalmRNA)==max(nanmean(WholeNoise(cc).AllTotalmRNA)));
    NormAPAvgs=[APAvgs./APAvgs(APbinToUse)];
    NoiseContributions(cc,1)=nanmean(WholeNoise(cc).AllIntraNoise(:,APbinToUse));
    NoiseContributions(cc,2)=nanmean(WholeNoise(cc).AllCoVarNoise(:,APbinToUse));
    NoiseContributions(cc,3)=nanmean(WholeNoise(cc).TotalNoise(:,APbinToUse));
    NoiseContributions(cc,4)=APbinToUse;
    NoiseContributions(cc,5)=WholeNoise(cc).Conf95InterNoise(APbinToUse);
    NoiseContributions(cc,6)=WholeNoise(cc).Conf95CoVar(APbinToUse);
    NoiseContributions(cc,7)=nanmedian(WholeNoise(cc).TotalNoise(:,APbinToUse));
    NoiseContributions(cc,8)=((NoiseContributions(cc,3)-NoiseContributions(cc,7))/(NoiseContributions(cc,3)))*100; % % mean is larger than median
    %Median values for intra,co-variance,total noise
    NoiseContributions(cc,9)=nanmedian(WholeNoise(cc).AllIntraNoise(:,APbinToUse));
    NoiseContributions(cc,10)=nanmedian(WholeNoise(cc).AllCoVarNoise(:,APbinToUse));
    NoiseContributions(cc,11)=nanmedian(WholeNoise(cc).TotalNoise(:,APbinToUse));
    %figuring out AP bin of 1/2 expression (anterior side)
    for aa=1:APbinToUse
        if (NormAPAvgs(aa) < 0.6) & (NormAPAvgs(aa+1) > 0.5)
            NoiseContributions(cc,12)=aa;
            break
        end
    end
    %Posterior end
    for aa=(APbinToUse+1):(length(APbinID)-1)
        if (NormAPAvgs(aa) < 0.6) & (NormAPAvgs(aa+1) < 0.5)
            NoiseContributions(cc,13)=aa;
            break
        end
    end
end
%% Percent nuclei w neg covariance f(x) of EL
for cc=[1:6,23]%length(ConstructList) %[1:6,23]%length(ConstructList)
    ContoUse=find(strcmp([ExpressionBins(2,:)],ConstructList{cc}));
    if length(ContoUse) >1
        ContoUse=ContoUse(1);
    end
    figure
    [ind1, ind2]=find(WholeNoise(cc).AllCoVarNoise<0);
    %get avg expression in each AP bin
    AvgAPExp = nanmean(WholeNoise(cc).AllTotalmRNA);
    HalfExp = prctile(AvgAPExp(NoiseContributions(cc,4)),75);
    %get ant half
    AntHalf = [AvgAPExp(1:NoiseContributions(cc,4)-1)] - HalfExp;
    [~,IdxA] = nanmin(abs(AntHalf));
    AntHalf = IdxA;  %started at 1 so no need to adjust number
    % posterior half
    PostHalf = [AvgAPExp(NoiseContributions(cc,4)+1:end)]-HalfExp;
    [~,IdxP] = nanmin(abs(PostHalf));
    PostHalf = IdxP+NoiseContributions(cc,4); % account for starting at peak 
    
for aa=1:length(APbinID)
    NegFract(aa)=length(find(ind2==aa));
    NegFract(aa)=NegFract(aa)/(sum(~isnan(WholeNoise(cc).AllCoVarNoise(:,aa))));
    %get rid of AP bins with 1 data point
    if sum(~isnan(WholeNoise(cc).AllCoVarNoise(:,aa))) < 2
        NegFract(aa) = nan;
    end
end
plot(EggLength, (NegFract.*100),'Color',Colors(cc).Color, 'LineWidth',2.5);
ylabel('% of nuclei with negative co-variance');
xlabel('% egg length');
xlim([0 100]);
ylim([0 35]);
hold on 
% mark peak expression
h=plot(EggLength(NoiseContributions(cc,4)),(NegFract(NoiseContributions(cc,4))*100),'o','LineWidth',2.5);
%Also want to mark 50% max exp boundaries
h_ant = plot(EggLength(AntHalf),(NegFract(AntHalf)*100),'*','LineWidth',2.5);
%h_ant7 = plot(EggLength(ExpressionBins{5,ContoUse}),(NegFract(ExpressionBins{5,ContoUse})*100),'v');
h_post= plot(EggLength(PostHalf),(NegFract(PostHalf)*100),'*');
%h_post7 = plot(EggLength(ExpressionBins{6,ContoUse}),(NegFract(ExpressionBins{6,ContoUse})*100),'v');
%h=line([EggLength(NoiseContributions(cc,4)) EggLength(NoiseContributions(cc,4))], [0 20]);
h.LineWidth=2; h.Color='k';
h_ant.LineWidth =2; h_ant.Color='k';  h_post.LineWidth=2; h_post.Color='k';
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 x_width y_width]);
set(gca,'Box','off', 'FontName',fontname,'FontSize',fontsize);
title(ConstructList{cc});
 if exist ('h_ant','var') 
     print([FigDirect filesep ConstructList{cc} 'NegCovarvEL_Bounds'],'-dsvg');
 else
     print( [FigDirect filesep ConstructList{cc} 'NegCovarvEL'],'-dsvg');
 end
end