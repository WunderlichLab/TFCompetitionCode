

%% Determine and graph total mRNA produced during nc14 by construct

% Load your constructs - listed as they are in DataStatus.xlsx
ConstructList= {'KrDist','KrProx','KrBothSep','KrDistEmpty','KrProxEmpty','KrDistDuplicN','KrProxDuplic','Kr2xProxEmpty','KrBoth','KrBothEmpty','KrDist32C','KrProx32C','KrBothSep32C','KrBoth32C','Kr2xProx32C','HbEmpty','KrDistDuplicN','KrDist17C','Kr2xDist32C','KrBoth17C','Kr2xDist17C','KrProx17C','Kr2xDistdBcd','Kr2xProxdHb','Kr2xProx17C','Kr2xDistEmpty','Kr2xDistEmpty32C','KrSEdBcd','KrInvSE','Kr2xDistLessMS2','KrSEdHb','KrEndogDist','KrEndogDist32C','KrEndogDist17C','KrSEChrom3','Kr2xDistChrom3','Kr2xProxChrom3','Kr2xDistChrom3_Empty','Kr4_1xBcd','KrSEChr3_Empty','Kr3_1xBcd','Kr4_1xHb','Kr5_1xHb','Kr3_1xHb','Kr4_1xZld', 'Kr5_1xStat92E','Kr3_1xZld','Kr3_1xStat92E','Kr4_Chr2_Chr3','Kr4_6xBcd','Kr3_Chr3_1xBcd','Kr5_1xBcd','KrBoth_BcdGFP','KrDist_BcdGFP','Kr3_Chr2_Chr3','Kr3_6xBcd','HbP2_HbP2','Kr4_1xTwi','HbP2_6xBcd','Kr1_6xBcd','Kr4_3xBcd','4F_Kr4_Kr4','Kr2_Kr2NoProm','Kr3_3xBcd','4F_Kr4_0', 'HbP2_1xBcd'} %{'KrDist','KrProx','KrBothSep', 'KrDistDuplicN', 'KrProxDuplic', 'KrBoth'};% %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, ~, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location

%Ask if only want nc14 info - respond 'y' 
ncUse=input('Want to only use nc14? y/n','s');
SlopeUse=input('Want to use Slope calculations? y/n','s');

%Organize all integrated fluorescence recordings for each construct
AvgmRNAProd=[];
for cc=1:length(ConstructList)
    Data= LoadMS2SetsCS(ConstructList{cc});
    Datalength(cc)=length(Data);
    NEmbryos = length(Data);
    %deal with if this is using the new pipeline 4/20/20 RW
    if contains(ConstructList{cc},'BcdGFP')
        APbinID=[Data(1).Particles.APbinID];
    else
        APbinID=[Data(1).APbinID];
    end
    Label = ConstructList(cc);
    AvgmRNAProdCon=[];
    firsttime=1;
    Timez=[];
    ConmRNAProdAllAP=[];
    ConmRNAProdSE=[];
    ConProdSD=[];
    EmbsArray=[];
     mRNAProdAllAP=[];
     mRNAProdErrorAllAP=[];
% Load the data for each embryo
    for ee=1:NEmbryos
        %Deal with this is from the new pipeline. 4/20/20
        if contains(ConstructList{cc},'BcdGFP')
            PrefixName=Data(ee).Particles.Prefix;
        else
            PrefixName=Data(ee).Prefix;
        end
        if SlopeUse=='y'
            filename=[DropboxFolder filesep PrefixName filesep 'BurstPropertiesSlope.mat'];
        elseif ncUse=='y'
        filename=[DropboxFolder filesep PrefixName filesep 'BurstPropertiesnc14.mat'];
        else
            filename=[DropboxFolder filesep PrefixName filesep 'BurstProperties.mat'];
        end 
        load(filename);
        CompPars=[DropboxFolder filesep PrefixName filesep 'CompiledParticles.mat'];
        load(CompPars);
        
        NumberBursts(ee,cc)=length([BurstProperties.Duration]);
        
        %seperate out by AP bin
        for aa=1:length(APbinID)
            ProdAP=[];
            ProdErrorAP=[];
            ProductionAP=find([BurstProperties.APBin]==APbinID(aa)); %find the nuclei with bursts 
            if isempty(ProductionAP)
                ProdAP=[ProdAP; nan];
            else
            for bb=1:length(ProductionAP)
                ProdAP=[ProdAP;BurstProperties(ProductionAP(bb)).TotalmRNA];  %put all mRNA outputs at a given AP value in a column going down
                ProdErrorAP=[ProdErrorAP;BurstProperties(ProductionAP(bb)).TotalmRNAError];           
            end
            end
            %make array with AP bins across as columns each page an embryo
            for bb=1:length(ProdAP)
                mRNAProdAllAP(bb,aa,ee)=ProdAP(bb);
            end
            for bb=1:length(ProdErrorAP)
                mRNAProdErrorAllAP(bb,aa,ee)=ProdErrorAP(bb);
            end
            mRNAProdAllAP(mRNAProdAllAP==0)=nan;  %Remove 0's recorded in BurstProperties to indicate nucleus exist without that spot of transcription active
            mRNAProdErrorAllAP(mRNAProdErrorAllAP==0)=nan;
            
            ProductionSD(ee,aa,cc)=nanstd(mRNAProdAllAP(:,aa,ee));
        ProductionSE(ee,aa,cc)=ProductionSD(ee,aa,cc)/sqrt(sum(~isnan(ProdAP))); %n=# of nuclei that produce mRNA
            clear ProductionAP
        end
        %Compile all raw duration data for a construct in one long column
        %to put in structure
        AvgProdAllAP(cc).nc14Time(ee)=ElapsedTime(end)-ElapsedTime(nc14);
        
         AvgProdAllAP(cc).EmbryosProd(ee).MeanProd=nanmean(mRNAProdAllAP(:,:,ee));
        AvgProdAllAP(cc).EmbryosProd(ee).SE=ProductionSE(ee,:,cc);
         AvgProdAllAP(cc).EmbryosProd(ee).SD=ProductionSD(ee,:,cc);
         AvgProdAllAP(cc).EmbryosProd(ee).AllProds=[mRNAProdAllAP(:,:,ee)];
         Timez=[Timez,(ElapsedTime(end)-ElapsedTime(nc14))];
    end
    %Combine the data from all embryos of a construct
    for bb=1:size(mRNAProdAllAP,3)
            ConmRNAProdAllAP=[ConmRNAProdAllAP; mRNAProdAllAP(:,:,bb)];
        end
        for bb=1:size(ProductionSD,3)
            ConProdSD=[ConProdSD; ProductionSD(:,:,bb)];
        end
        for bb=1:size(ProductionSE,3)
            ConmRNAProdSE=[ConmRNAProdSE;ProductionSE(:,:,bb)];
        end

        AvgProdAllAP(cc).nc14Time=mean(Timez);
        AvgProdAllAP(cc).AvgProd=nanmean(ConmRNAProdAllAP,1);  %Avg mRNA production of all embryos of a construct by AP position
        AvgProdAllAP(cc).ConSD=nanmean(ConProdSD,1);
        AvgProdAllAP(cc).ConSE=nanmean(ConmRNAProdSE,1);
        AvgProdAllAP(cc).AllProds=[ConmRNAProdAllAP];
        AvgProdAllAP(cc).AllSD=nanstd([AvgProdAllAP(cc).AllProds]);
        for aa=1:length(APbinID)
            %Perform SE calc by AP bin as there are different # of data
            %points in different AP bins
            AvgProdAllAP(cc).AllSE(aa)=(([AvgProdAllAP(cc).AllSD(aa)])./sqrt(sum(~isnan(AvgProdAllAP(cc).AllProds(:,aa)))));
        end
        AvgProdAllAP(cc).All95Conf=(AvgProdAllAP(cc).AllSE).*1.95;
        AvgProdAllAP(cc).EmbryoMeans=[EmbsArray];
        clear mRNAProdAllAP mRNAProdErrorAllAP;
end

% Get rid of single values for an AP bin
for cc=1:length(ConstructList)
    for aa=1:length(APbinID)
        if sum(~isnan(AvgProdAllAP(cc).AllProds(:,aa)))==1
            AvgProdAllAP(cc).AllSD(aa)=nan;
            AvgProdAllAP(cc).AvgProd(aa)=nan;
            AvgProdAllAP(cc).AllSE(aa)=nan;
            AvgProdAllAP(cc).All95Conf(aa)=nan;
        end
    end
end


    

%% Plotting 
EggLength=APbinID.*100;
FractInfo=[DropboxFolder filesep 'Constructs' filesep 'FractON.mat'];

% 3/3/2019 updated/corrected Frnap value to be 379 instead of 377
% MS2Conversion =F1  which for Wunderlich Nikon scope is = 1338

MS2Conversion=(379*((1.275+4.021)/1.5)); %Frnap * Telongation ((length of MS2 +rest of transcript) / elongation rate)
HbMS2Conversion=(379*((1.275+4.021)/1.5));

%Figuring out issues w laser 
MS2Conversion_New=(480*((1.275+4.021)/1.5));
for cc=1:length(ConstructList)
    AvgProdAllAP(cc).AvgmRNACounts=[[AvgProdAllAP(cc).AvgProd]./MS2Conversion];
end

% Set construct colors
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


Colors(1).Color=DistalColor; 
Colors(2).Color=ProxColor; 
Colors(3).Color=BothSepColor; 
Colors(4).Color=DistalColor;
Colors(5).Color=ProxColor;
Colors(6).Color=DoubDistColor; 
Colors(7).Color=DoubProxColor;
Colors(8).Color=DoubProxColor;
Colors(9).Color=BothColor;
Colors(10).Color=BothColor;
Colors(11).Color=DistalColor;
Colors(12).Color=ProxColor;
Colors(13).Color=BothSepColor;
Colors(14).Color=BothColor;
Colors(15).Color=DoubProxColor;
Colors(16).Color='k';
Colors(17).Color=DoubDistColor;
Colors(18).Color=DistalColor;
Colors(19).Color=DoubDistColor;
Colors(20).Color=BothColor;
Colors(21).Color=DoubDistColor;
Colors(22).Color=ProxColor;
Colors(23).Color=DoubDistColor;
Colors(24).Color=DoubProxColor;
Colors(25).Color=DoubProxColor;
Colors(26).Color=DoubDistColor;
Colors(27).Color=DoubDistColor;
Colors(28).Color=BothColor;
Colors(29).Color=BothColor;
Colors(30).Color=DoubDistColor;
Colors(31).Color=BothColor;
Colors(32).Color=DistalColor;
Colors(33).Color=DistalColor;
Colors(34).Color=DistalColor;
Colors(35).Color=BothColor;
Colors(36).Color = DoubDistColor;
Colors(37).Color = DoubProxColor;
%Colors(38).Color = DoubDistColor;
Colors(38).Color = DoubDistColor;
Colors(39).Color = DoubDistColor;
Colors(40).Color = BothColor;
Colors(41).Color=BothColor;
Colors(42).Color = DoubDistColor;
Colors(43).Color = DoubProxColor;
Colors(44).Color=BothColor;
Colors(49).Color=DoubDistColor;


% Set font/display parameters
FontUsed=input('Want larger font?','s');
if FontUsed=='y'
    fontsize=15;
else
fontsize=10;
end
fontname='Arial';
x_width=3; y_width=2.25;
x_widthsplit=1.5; y_widthsplit=1.125;
xSize = 7; ySize = 6; xLeft = 0.5; yTop = 0.5;

FigDirect=[DropboxFolder filesep 'Figures'];
save([DropboxFolder filesep 'Constructs' filesep 'AllTotalmRNAProd'],'AvgProdAllAP');
%% AP bins of max expression each construct
% Record the AP bin of max expression for each construct to have for other
% calculations (ie noise)
ExpressionBins=cell(2,(length(ConstructList)));
for cc=1:length(ConstructList)
    AvgProdAP=[nanmean(AvgProdAllAP(cc).AllProds)];
    APBintoUse=find((nanmean(AvgProdAllAP(cc).AllProds)==max(nanmean(AvgProdAllAP(cc).AllProds))));
    ExpressionBins{1,cc}=APBintoUse;
    ExpressionBins{2,cc}=ConstructList{cc};
end
save([DropboxFolder filesep 'Constructs' filesep 'ConstructExpressionBins'],'ExpressionBins');

%% Create figures for core constructs
for cc=[32,2,6,7,9];
    %figure
 errorbar(EggLength,AvgProdAllAP(cc).AvgProd,AvgProdAllAP(cc).All95Conf,'Color',Colors(cc).Color,'LineWidth',2.5);
 hold on 
 yyaxis right
 plot(EggLength, (AvgProdAllAP(cc).AvgProd./MS2Conversion), 'Color', Colors(cc).Color,'LineWidth',2.5);
 set(gca, 'FontSize', fontsize, 'FontName', fontname,'YColor','k');
 xlabel('% egg length')
 ylabel('total transcripts produced')
 %title('Avg mRNA production per nucleus') 
 ylim([0 950000./MS2Conversion]);
 xlim([0 100]);
 yyaxis left 
 ylabel('integrated fluorescence (AU)');
 ylim([0 950000]);
 print( [FigDirect filesep ConstructList{cc} 'TotalmRNA'],'-dpdf');
end

%% Singles/duplicated vs shadow enhancer pair
 figure 
 hold on 
 errorbar(EggLength, AvgProdAllAP(17).AvgProd, AvgProdAllAP(17).All95Conf,'Color',Colors(17).Color,'LineWidth',5.5);
 errorbar(EggLength, (AvgProdAllAP(7).AvgProd),(AvgProdAllAP(7).All95Conf),'Color',Colors(7).Color,'LineWidth',5.5);
 errorbar(EggLength, AvgProdAllAP(9).AvgProd,AvgProdAllAP(9).All95Conf,'Color',Colors(9).Color,'LineWidth',5.5);
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
 set(gcf,'PaperUnits','inches');
 set(gcf,'PaperPosition',[xLeft yTop 7 6]);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 xlim([0 100]);
 ylim([0 1000000]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(17).AvgProd./MS2Conversion),(AvgProdAllAP(17).All95Conf./MS2Conversion),'Color',Colors(17).Color,'LineWidth',5.5);
 errorbar(EggLength, (AvgProdAllAP(7).AvgProd./MS2Conversion),(AvgProdAllAP(7).All95Conf./MS2Conversion),'Color',Colors(7).Color,'LineWidth',5.5);
 errorbar(EggLength, (AvgProdAllAP(9).AvgProd./MS2Conversion),(AvgProdAllAP(9).All95Conf./MS2Conversion),'Color',Colors(9).Color,'LineWidth',5.5);
ylabel('total transcripts produced');
ylim([0 (1000000/MS2Conversion)]);
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'DuplicatedvsBothTotalmRNA','.pdf'],'pdf');

% 2x vs SE limit to 20-80% egg length 
figure 
 hold on 
 errorbar(EggLength, AvgProdAllAP(17).AvgProd, AvgProdAllAP(17).All95Conf,'Color',Colors(17).Color,'LineWidth',5.5);
 errorbar(EggLength, (AvgProdAllAP(7).AvgProd),(AvgProdAllAP(7).All95Conf),'Color',Colors(7).Color,'LineWidth',5.5);
 errorbar(EggLength, AvgProdAllAP(9).AvgProd,AvgProdAllAP(9).All95Conf,'Color',Colors(9).Color,'LineWidth',5.5);
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
 set(gcf,'PaperUnits','inches');
 set(gcf,'PaperPosition',[xLeft yTop 7 6]);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 xlim([20 80]);
 ylim([0 1000000]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(17).AvgProd./MS2Conversion),(AvgProdAllAP(17).All95Conf./MS2Conversion),'Color',Colors(17).Color,'LineWidth',5.5);
 errorbar(EggLength, (AvgProdAllAP(7).AvgProd./MS2Conversion),(AvgProdAllAP(7).All95Conf./MS2Conversion),'Color',Colors(7).Color,'LineWidth',5.5);
 errorbar(EggLength, (AvgProdAllAP(9).AvgProd./MS2Conversion),(AvgProdAllAP(9).All95Conf./MS2Conversion),'Color',Colors(9).Color,'LineWidth',5.5);
ylabel('total transcripts produced');
ylim([0 (1000000/MS2Conversion)]);
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'DuplicatedvsBothTotalmRNA_2080Limit','.pdf'],'pdf');

 figure 
 hold on 
 errorbar(EggLength, AvgProdAllAP(32).AvgProd, AvgProdAllAP(32).All95Conf,'Color',Colors(32).Color,'LineWidth',2.5);
 errorbar(EggLength, (AvgProdAllAP(2).AvgProd),(AvgProdAllAP(2).All95Conf),'Color',Colors(2).Color,'LineWidth',2.5);
 errorbar(EggLength, AvgProdAllAP(9).AvgProd,AvgProdAllAP(9).All95Conf,'Color',Colors(9).Color,'LineWidth',2.5);
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 set(gcf,'PaperUnits','inches');
 set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 xlim([0 100]);
 ylim([0 1000000]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(32).AvgProd./MS2Conversion),(AvgProdAllAP(32).All95Conf./MS2Conversion),'Color',Colors(32).Color,'LineWidth',2.5);
 errorbar(EggLength, (AvgProdAllAP(2).AvgProd./MS2Conversion),(AvgProdAllAP(2).All95Conf./MS2Conversion),'Color',Colors(2).Color,'LineWidth',2.5);
 errorbar(EggLength, (AvgProdAllAP(9).AvgProd./MS2Conversion),(AvgProdAllAP(9).All95Conf./MS2Conversion),'Color',Colors(9).Color,'LineWidth',2.5);
ylabel('total transcripts produced');
ylim([0 (1000000/MS2Conversion)]);
%set(gca,'FontSize',fontsize,'FontName',fontname,'YColor','k');
%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 x_width y_width]);
set(gca,'Box','off', 'FontName',fontname,'FontSize',fontsize,'YColor','k');
 print( [FigDirect filesep 'SinglesvsBothTotalmRNA_ED'],'-dpdf');
 
 %Just singles
 figure 
 hold on 
 errorbar(EggLength, AvgProdAllAP(32).AvgProd, AvgProdAllAP(32).All95Conf,'Color',Colors(32).Color,'LineWidth',2.5);
 errorbar(EggLength, (AvgProdAllAP(2).AvgProd),(AvgProdAllAP(2).All95Conf),'Color',Colors(2).Color,'LineWidth',2.5);
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 set(gcf,'PaperUnits','inches');
 set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 xlim([0 100]);
 ylim([0 1000000]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(32).AvgProd./MS2Conversion),(AvgProdAllAP(32).All95Conf./MS2Conversion),'Color',Colors(32).Color,'LineWidth',2.5);
 errorbar(EggLength, (AvgProdAllAP(2).AvgProd./MS2Conversion),(AvgProdAllAP(2).All95Conf./MS2Conversion),'Color',Colors(2).Color,'LineWidth',2.5);
ylabel('total transcripts produced');
ylim([0 (1000000/MS2Conversion)]);
%set(gca,'FontSize',fontsize,'FontName',fontname,'YColor','k');
%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 x_width y_width]);
set(gca,'Box','off', 'FontName',fontname,'FontSize',fontsize,'YColor','k');
 print( [FigDirect filesep 'SinglesTotalmRNA_ED'],'-dpdf');
 
 %% Compare 2x vs SE 20-80% with indications of Kr expression limits 
 % 2x vs SE limit to 20-80% egg length 
 Cutoff = input('What percetage to use for cut off?'); %Where to put lines for Kr expression boundaries
 %Use the SE pair's expression to calculate the boundaries
 
 for aa = 1:length(APbinID)
     if AvgProdAllAP(9).AvgProd(aa) >= (nanmax([AvgProdAllAP(9).AvgProd]) *(Cutoff/100));
        AntBound = aa;
        break
     end
 end
 
 for aa = ExpressionBins{1,9}:length(APbinID)
     if AvgProdAllAP(9).AvgProd(aa) <= (nanmax([AvgProdAllAP(9).AvgProd]) *(Cutoff/100));
        PostBound = aa;
        break
     end
 end
 
figure 
 hold on 
 errorbar(EggLength, AvgProdAllAP(17).AvgProd, AvgProdAllAP(17).All95Conf,'Color',Colors(17).Color,'LineWidth',5.5);
 errorbar(EggLength, (AvgProdAllAP(7).AvgProd),(AvgProdAllAP(7).All95Conf),'Color',Colors(7).Color,'LineWidth',5.5);
 errorbar(EggLength, AvgProdAllAP(9).AvgProd,AvgProdAllAP(9).All95Conf,'Color',Colors(9).Color,'LineWidth',5.5);
 %plot the boundaries of Kr pattern
 plot([EggLength(AntBound) EggLength(AntBound)], [0 AvgProdAllAP(9).AvgProd(AntBound)],'Color','k','LineWidth',2.5);
plot([EggLength(PostBound) EggLength(PostBound)], [0 AvgProdAllAP(9).AvgProd(PostBound)],'Color','k','LineWidth',2.5);
 set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
 set(gcf,'PaperUnits','inches');
 set(gcf,'PaperPosition',[xLeft yTop 7 6]);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 xlim([20 80]);
 ylim([0 1000000]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(17).AvgProd./MS2Conversion),(AvgProdAllAP(17).All95Conf./MS2Conversion),'Color',Colors(17).Color,'LineWidth',5.5);
 errorbar(EggLength, (AvgProdAllAP(7).AvgProd./MS2Conversion),(AvgProdAllAP(7).All95Conf./MS2Conversion),'Color',Colors(7).Color,'LineWidth',5.5);
 errorbar(EggLength, (AvgProdAllAP(9).AvgProd./MS2Conversion),(AvgProdAllAP(9).All95Conf./MS2Conversion),'Color',Colors(9).Color,'LineWidth',5.5);
ylabel('total transcripts produced');
ylim([0 (1000000/MS2Conversion)]);
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep ['DupvsBothTotalmRNA_' num2str(Cutoff) 'Thresh_2080Limit'],'.pdf'],'pdf');
 %% Compare Endogenous Distal spacing
figure
  %yyaxis left
 d= errorbar(EggLength,AvgProdAllAP(1).AvgProd,AvgProdAllAP(1).All95Conf,'Color',Colors(1).Color,'LineStyle',':','LineWidth',5);
 d.Bar.LineStyle ='dotted';
 hold on 
 errorbar(EggLength,AvgProdAllAP(32).AvgProd, AvgProdAllAP(32).All95Conf,'Color',Colors(32).Color,'LineWidth',5);
 yyaxis right
 plot(EggLength, (AvgProdAllAP(1).AvgProd./MS2Conversion), 'Color', Colors(1).Color,'LineStyle',':','LineWidth',5);
 plot(EggLength, (AvgProdAllAP(32).AvgProd./MS2Conversion), 'Color', Colors(32).Color,'LineWidth',5);
 set(gca, 'FontSize', fontsize, 'FontName', fontname,'YColor','k');
 set(gcf,'PaperUnits', 'inches');
set(gcf,'PaperUnits','inches','PaperPosition',[xLeft yTop xSize ySize]);
 xlabel('% egg length')
 ylabel('total transcripts produced')
 %title('Avg mRNA production per nucleus') 
 ylim([0 900000./MS2Conversion]);
 xlim([0 100]);
 yyaxis left 
 ylabel('integrated fluorescence (AU)');
 ylim([0 900000]);
print( [FigDirect filesep 'EndogvsOrigDistTotalmRNA'],'-dpdf');

% Altered color
figure
  %yyaxis left
 d= errorbar(EggLength,AvgProdAllAP(1).AvgProd,AvgProdAllAP(1).All95Conf,'Color',PDistalColor,'LineStyle',':','LineWidth',5);
 %d.Bar.LineStyle ='dotted';
 hold on 
 errorbar(EggLength,AvgProdAllAP(32).AvgProd, AvgProdAllAP(32).All95Conf,'Color',Colors(32).Color,'LineWidth',5);
 yyaxis right
 plot(EggLength, (AvgProdAllAP(1).AvgProd./MS2Conversion), 'Color', PDistalColor,'LineStyle',':','LineWidth',5);
 plot(EggLength, (AvgProdAllAP(32).AvgProd./MS2Conversion), 'Color', Colors(32).Color,'LineWidth',5);
 set(gca, 'FontSize', fontsize, 'FontName', fontname,'YColor','k');
 set(gcf,'PaperUnits', 'inches');
set(gcf,'PaperUnits','inches','PaperPosition',[xLeft yTop xSize ySize]);
 xlabel('% egg length')
 ylabel('total transcripts produced')
 %title('Avg mRNA production per nucleus') 
 ylim([0 900000./MS2Conversion]);
 xlim([0 100]);
 yyaxis left 
 ylabel('integrated fluorescence (AU)');
 ylim([0 900000]);
print( [FigDirect filesep 'EndogvsOrigDistTotalmRNA_AltColor'],'-dpdf');

%% Compare all at once
 figure 
 hold on 
 errorbar(EggLength, AvgProdAllAP(17).AvgProd, AvgProdAllAP(17).All95Conf,'Color',Colors(17).Color,'LineWidth',2.5);
 errorbar(EggLength, (AvgProdAllAP(7).AvgProd),(AvgProdAllAP(7).All95Conf),'Color',Colors(7).Color,'LineWidth',2.5);
 errorbar(EggLength, AvgProdAllAP(9).AvgProd,AvgProdAllAP(9).All95Conf,'Color',Colors(9).Color,'LineWidth',2.5);
 errorbar(EggLength, AvgProdAllAP(32).AvgProd, AvgProdAllAP(32).All95Conf,'Color',Colors(32).Color,'LineWidth',2.5);
 errorbar(EggLength, (AvgProdAllAP(2).AvgProd),(AvgProdAllAP(2).All95Conf),'Color',Colors(2).Color,'LineWidth',2.5);
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 set(gcf,'PaperUnits','inches');
 set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 xlim([0 100]);
 ylim([0 1000000]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(17).AvgProd./MS2Conversion),(AvgProdAllAP(17).All95Conf./MS2Conversion),'Color',Colors(17).Color,'LineWidth',2.5);
 errorbar(EggLength, (AvgProdAllAP(7).AvgProd./MS2Conversion),(AvgProdAllAP(7).All95Conf./MS2Conversion),'Color',Colors(7).Color,'LineWidth',2.5);
 errorbar(EggLength, (AvgProdAllAP(9).AvgProd./MS2Conversion),(AvgProdAllAP(9).All95Conf./MS2Conversion),'Color',Colors(9).Color,'LineWidth',2.5);
 errorbar(EggLength, (AvgProdAllAP(32).AvgProd./MS2Conversion),(AvgProdAllAP(32).All95Conf./MS2Conversion),'Color',Colors(32).Color,'LineWidth',2.5);
 errorbar(EggLength, (AvgProdAllAP(2).AvgProd./MS2Conversion),(AvgProdAllAP(2).All95Conf./MS2Conversion),'Color',Colors(2).Color,'LineWidth',2.5);
 ylabel('total transcripts produced');
ylim([0 (1000000/MS2Conversion)]);
%set(gca,'FontSize',fontsize,'FontName',fontname,'YColor','k');
%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 x_width y_width]);
set(gca,'Box','off', 'FontName',fontname,'FontSize',fontsize,'YColor','k');
%% Compare 2xDistal with lost MS2 loops
figure
 errorbar(EggLength,AvgProdAllAP(6).AvgProd,AvgProdAllAP(6).All95Conf,'Color',Colors(6).Color,'LineWidth',2.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(30).AvgProd, AvgProdAllAP(30).All95Conf,'Color',Colors(30).Color,'LineWidth',2.5,'LineStyle',':');
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 %title('Avg mRNA production per nucleus') 
 ylim([0 1000000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(6).AvgProd./MS2Conversion),(AvgProdAllAP(6).All95Conf./MS2Conversion),'Color',Colors(6).Color,'LineWidth',2.5,'LineStyle','-');
  errorbar(EggLength, (AvgProdAllAP(30).AvgProd./MS2Conversion),(AvgProdAllAP(30).All95Conf./MS2Conversion),'Color',Colors(30).Color,'LineWidth',2.5,'LineStyle',':');
ylabel('total transcripts produced');
ylim([0 (1000000/MS2Conversion)]);

%% Assess additivity of single enhancers vs duplicated 
% Distal Enhancer 
figure
errorbar(EggLength,[AvgProdAllAP(1).AvgProd].*2,AvgProdAllAP(1).All95Conf,'Color',Colors(6).Color,'LineWidth',5.5,'LineStyle',':');
 hold on 
 errorbar(EggLength,AvgProdAllAP(6).AvgProd, AvgProdAllAP(6).All95Conf,'Color',Colors(6).Color,'LineWidth',5.5);
 set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
 set(gcf,'PaperUnits','inches');
 set(gcf,'PaperPosition',[xLeft yTop 7 6]);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 %title('Avg mRNA production per nucleus') 
 ylim([0 1500000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (([AvgProdAllAP(1).AvgProd].*2)./MS2Conversion),(AvgProdAllAP(1).All95Conf./MS2Conversion),'Color',Colors(6).Color,'LineWidth',5.5,'LineStyle',':');
  errorbar(EggLength, (AvgProdAllAP(6).AvgProd./MS2Conversion),(AvgProdAllAP(6).All95Conf./MS2Conversion),'Color',Colors(6).Color,'LineWidth',5.5,'LineStyle','-');
ylabel('total transcripts produced');
ylim([0 (1500000/MS2Conversion)]);
set(gca,'Box','off', 'FontName',fontname,'FontSize',fontsize,'YColor','k');
saveas(gcf, [FigDirect filesep 'DupDistvsPredictedAdditive','.pdf'],'pdf');


%Proximal enhancer
figure
errorbar(EggLength,[AvgProdAllAP(2).AvgProd].*2,AvgProdAllAP(2).All95Conf,'Color',Colors(7).Color,'LineWidth',5.5,'LineStyle',':');
 hold on 
 errorbar(EggLength,AvgProdAllAP(7).AvgProd, AvgProdAllAP(7).All95Conf,'Color',Colors(7).Color,'LineWidth',5.5);
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
 set(gcf,'PaperUnits','inches');
 set(gcf,'PaperPosition',[xLeft yTop 7 6]); xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 %title('Avg mRNA production per nucleus') 
 ylim([0 1500000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (([AvgProdAllAP(2).AvgProd].*2)./MS2Conversion),(AvgProdAllAP(2).All95Conf./MS2Conversion),'Color',Colors(7).Color,'LineWidth',5.5,'LineStyle',':');
  errorbar(EggLength, (AvgProdAllAP(7).AvgProd./MS2Conversion),(AvgProdAllAP(7).All95Conf./MS2Conversion),'Color',Colors(7).Color,'LineWidth',5.5,'LineStyle','-');
ylabel('total transcripts produced');
ylim([0 (1500000/MS2Conversion)]);
set(gca,'Box','off', 'FontName',fontname,'FontSize',fontsize,'YColor','k');
saveas(gcf, [FigDirect filesep 'DupProxvsPredictedAdditive','.pdf'],'pdf');


%% Compare chromosome III constructs
figure 
Gray = [0.5 0.5 0.5];
 hold on 
 errorbar(EggLength, AvgProdAllAP(9).AvgProd, AvgProdAllAP(9).All95Conf,'Color',Colors(9).Color,'LineWidth',2.5);
 errorbar(EggLength, AvgProdAllAP(35).AvgProd,AvgProdAllAP(35).All95Conf,'Color',Colors(35).Color,'LineWidth',2.5,'LineStyle','--');
  errorbar(EggLength, AvgProdAllAP(55).AvgProd,AvgProdAllAP(55).All95Conf,'Color','r','LineWidth',2.5,'LineStyle','--'); 
 set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
 set(gcf,'PaperUnits','inches');
 set(gcf,'PaperPosition',[xLeft yTop 7 6]);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 xlim([0 100]);
 ylim([0 1000000]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(9).AvgProd./MS2Conversion),(AvgProdAllAP(9).All95Conf./MS2Conversion),'Color',Colors(9).Color,'LineWidth',2.5);
 errorbar(EggLength, (AvgProdAllAP(35).AvgProd./MS2Conversion),(AvgProdAllAP(35).All95Conf./MS2Conversion),'Color',Colors(35).Color,'LineWidth',2.5,'LineStyle','--');
 errorbar(EggLength, (AvgProdAllAP(55).AvgProd./MS2Conversion),(AvgProdAllAP(55).All95Conf./MS2Conversion),'Color','r','LineWidth',2.5,'LineStyle','--');
 ylabel('total transcripts produced');
ylim([0 (1000000/MS2Conversion)]);
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'Chr2v3_Both_TotmRNA','.svg'],'svg');

%Hemizygotes
figure 
 hold on 
 errorbar(EggLength, AvgProdAllAP(35).AvgProd, AvgProdAllAP(35).All95Conf,'Color',Colors(35).Color,'LineWidth',2.5);
 errorbar(EggLength, AvgProdAllAP(40).AvgProd,AvgProdAllAP(40).All95Conf,'Color',Colors(40).Color,'LineWidth',2.5,'LineStyle',':');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
 set(gcf,'PaperUnits','inches');
 set(gcf,'PaperPosition',[xLeft yTop 7 6]);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 xlim([0 100]);
 ylim([0 1200000]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(35).AvgProd./MS2Conversion),(AvgProdAllAP(35).All95Conf./MS2Conversion),'Color',Colors(35).Color,'LineWidth',2.5);
 errorbar(EggLength, (AvgProdAllAP(40).AvgProd./MS2Conversion),(AvgProdAllAP(40).All95Conf./MS2Conversion),'Color',Colors(40).Color,'LineWidth',2.5,'LineStyle',':');
ylabel('total transcripts produced');
ylim([0 (1200000/MS2Conversion)]);
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);

%% 2xDistal Ch II vs III
figure 
 hold on 
 errorbar(EggLength, AvgProdAllAP(6).AvgProd, AvgProdAllAP(6).All95Conf,'Color',Colors(6).Color,'LineWidth',2.5);
 errorbar(EggLength, AvgProdAllAP(36).AvgProd,AvgProdAllAP(36).All95Conf,'Color',Colors(36).Color,'LineWidth',2.5,'LineStyle','--');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
 set(gcf,'PaperUnits','inches');
 set(gcf,'PaperPosition',[xLeft yTop 7 6]);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 xlim([0 100]);
 ylim([0 1000000]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(6).AvgProd./MS2Conversion),(AvgProdAllAP(6).All95Conf./MS2Conversion),'Color',Colors(6).Color,'LineWidth',2.5);
 errorbar(EggLength, (AvgProdAllAP(36).AvgProd./MS2Conversion),(AvgProdAllAP(36).All95Conf./MS2Conversion),'Color',Colors(36).Color,'LineWidth',2.5,'LineStyle','--');
ylabel('total transcripts produced');
ylim([0 (1000000/MS2Conversion)]);
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'Chr2v3_2xDist_TotmRNA','.svg'],'svg');

%Compare hemizygotes
figure 
 hold on 
 errorbar(EggLength, AvgProdAllAP(36).AvgProd, AvgProdAllAP(36).All95Conf,'Color',Colors(36).Color,'LineWidth',2.5);
 errorbar(EggLength, AvgProdAllAP(38).AvgProd,AvgProdAllAP(38).All95Conf,'Color',Colors(38).Color,'LineWidth',2.5,'LineStyle',':');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
 set(gcf,'PaperUnits','inches');
 set(gcf,'PaperPosition',[xLeft yTop 7 6]);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 xlim([0 100]);
 ylim([0 1200000]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(36).AvgProd./MS2Conversion),(AvgProdAllAP(36).All95Conf./MS2Conversion),'Color',Colors(36).Color,'LineWidth',2.5);
 errorbar(EggLength, (AvgProdAllAP(38).AvgProd./MS2Conversion),(AvgProdAllAP(38).All95Conf./MS2Conversion),'Color',Colors(38).Color,'LineWidth',2.5,'LineStyle',':');
ylabel('total transcripts produced');
ylim([0 (1200000/MS2Conversion)]);
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
%% 2xProximal Ch II vs III
figure 
 hold on 
 errorbar(EggLength, AvgProdAllAP(7).AvgProd, AvgProdAllAP(7).All95Conf,'Color',Colors(7).Color,'LineWidth',2.5);
 errorbar(EggLength, AvgProdAllAP(37).AvgProd,AvgProdAllAP(37).All95Conf,'Color',Colors(37).Color,'LineWidth',2.5,'LineStyle','--');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
 set(gcf,'PaperUnits','inches');
 set(gcf,'PaperPosition',[xLeft yTop 7 6]);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 xlim([0 100]);
 ylim([0 1000000]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(7).AvgProd./MS2Conversion),(AvgProdAllAP(7).All95Conf./MS2Conversion),'Color',Colors(7).Color,'LineWidth',2.5);
 errorbar(EggLength, (AvgProdAllAP(37).AvgProd./MS2Conversion),(AvgProdAllAP(37).All95Conf./MS2Conversion),'Color',Colors(37).Color,'LineWidth',2.5,'LineStyle','--');
ylabel('total transcripts produced');
ylim([0 (1000000/MS2Conversion)]);
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'Chr2v3_2xProx_TotmRNA','.svg'],'svg');

%% 
figure 
 hold on 
 errorbar(EggLength, AvgProdAllAP(4).AvgProd, AvgProdAllAP(4).All95Conf,'Color',Colors(4).Color,'LineWidth',2.5,'LineStyle','--');
 errorbar(EggLength, (AvgProdAllAP(5).AvgProd),(AvgProdAllAP(5).All95Conf),'Color',Colors(5).Color,'LineWidth',2.5,'LineStyle','--');
  errorbar(EggLength, (AvgProdAllAP(3).AvgProd),(AvgProdAllAP(3).All95Conf),'Color',Colors(3).Color,'LineWidth',2.5);
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 set(gcf,'PaperUnits','inches');
 set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 xlim([0 100]);
 ylim([0 1000000]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(4).AvgProd./MS2Conversion),(AvgProdAllAP(4).All95Conf./MS2Conversion),'Color',Colors(4).Color,'LineWidth',2.5);
 errorbar(EggLength, (AvgProdAllAP(5).AvgProd./MS2Conversion),(AvgProdAllAP(5).All95Conf./MS2Conversion),'Color',Colors(5).Color,'LineWidth',2.5);
 errorbar(EggLength, (AvgProdAllAP(3).AvgProd./MS2Conversion),(AvgProdAllAP(3).All95Conf./MS2Conversion),'Color',Colors(3).Color,'LineWidth',2.5);
 ylabel('total transcripts produced');
ylim([0 (1000000/MS2Conversion)]);
%set(gca,'FontSize',fontsize,'FontName',fontname,'YColor','k');
%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 x_width y_width]);
set(gca,'Box','off', 'FontName',fontname,'FontSize',fontsize,'YColor','k');
%% Find AP bins of 50%, 70% expression 
Cutoff = input('What fraction max expression to use?');
for cc = 1:length(ConstructList)
    CutoffExp = (max(AvgProdAllAP(cc).AvgProd) * Cutoff);
    for aa = 1:ExpressionBins{1,cc}%length(APbinID)  
        if AvgProdAllAP(cc).AvgProd(aa) >= CutoffExp;
            AntBound(cc) = aa-1;
            break 
        end
%         if (AvgProdAllAP(cc).AvgProd(aa)) < CutoffExp;
%             continue
%         elseif (AvgProdAllAP(cc).AvgProd(aa)) > CutoffExp;
%             AntBound(cc) = aa;
%         else
%             AntBound(cc) = aa-1;
%         
%         end
    end

    for bb = ExpressionBins{1,cc}:length(APbinID)
        if AvgProdAllAP(cc).AvgProd(bb) <= CutoffExp;
            PostBound(cc) = bb;
            break
        end
%         if AvgProdAllAP(cc).AvgProd(bb) > (max(AvgProdAllAP(cc).AvgProd)*Cutoff);
%             continue
%         elseif AvgProdAllAP(cc).AvgProd(bb) <= (max(AvgProdAllAP(cc).AvgProd)*Cutoff);
%             PostBound(cc) = bb;
%         end
    end
    if Cutoff == 0.5
        ExpressionBins{3,cc} = AntBound(cc);
        ExpressionBins{4,cc}= PostBound(cc);
    end
    if Cutoff == 0.7
        ExpressionBins{5,cc}=AntBound(cc);
        ExpressionBins{6,cc}=PostBound(cc);
    end
    
end


 %%
 %Compare single and doubles
 figure 
 errorbar(EggLength, AvgProdAllAP(1).AvgProd, AvgProdAllAP(1).All95Conf,'Color',Colors(1).Color,'LineWidth',2.5);
 hold on 
 errorbar(EggLength, AvgProdAllAP(17).AvgProd, AvgProdAllAP(17).All95Conf,'Color',Colors(17).Color,'LineWidth',2.5);
 legend('Distal', '2x Distal');
 xlabel('% Egg length')
 ylabel('integrated fluorescence')
 title('Avg mRNA production per nucleus')
 set(gca,'FontSize',fontsize,'FontName',fontname);
 ylim([0 700000]);
 xlim([0 100]);
 
 figure
 errorbar(EggLength,AvgProdAllAP(2).AvgProd,AvgProdAllAP(2).All95Conf,'Color',Colors(2).Color,'LineWidth',1.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(7).AvgProd, AvgProdAllAP(7).All95Conf,'Color',Colors(7).Color,'LineWidth',1.5);
 legend('Proximal',  '2x Proximal');
 xlabel('% Egg length')
 ylabel('integrated fluorescence')
 title('Avg mRNA production per nucleus') 
 ylim([0 700000]);
 xlim([0 100]);
 
 figure
 yyaxis left
 errorbar(EggLength,AvgProdAllAP(1).AvgProd,AvgProdAllAP(1).All95Conf,'Color',Colors(1).Color,'LineWidth',1.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(2).AvgProd, AvgProdAllAP(2).All95Conf,'Color',Colors(2).Color,'LineWidth',1.5);
 errorbar(EggLength,AvgProdAllAP(3).AvgProd, AvgProdAllAP(3).All95Conf,'Color',Colors(3).Color,'LineWidth',1.5);
 %legend('Distal','Proximal','Both Sep');
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% Egg length')
 ylabel('integrated fluorescence')
 %title('Avg mRNA production per nucleus') 
 ylim([0 700000]);
 xlim([0 100]);
 yyaxis right 
 errorbar(EggLength,(AvgProdAllAP(1).AvgProd./(MS2Conversion)),(AvgProdAllAP(1).All95Conf./(MS2Conversion)),'Color',Colors(1).Color,'LineWidth',1.5,'LineStyle','-');
 hold on 
 errorbar(EggLength,(AvgProdAllAP(2).AvgProd./(MS2Conversion)), (AvgProdAllAP(2).All95Conf./(MS2Conversion)),'Color',Colors(2).Color,'LineWidth',1.5,'LineStyle','-');
 errorbar(EggLength,(AvgProdAllAP(3).AvgProd./(MS2Conversion)), (AvgProdAllAP(3).All95Conf./(MS2Conversion)),'Color',Colors(3).Color,'LineWidth',1.5,'LineStyle','-');
 ylabel('Total mRNA molecules');
 ylim([0 (700000/MS2Conversion)]);
 
  figure
  %yyaxis left
 errorbar(EggLength,AvgProdAllAP(1).AvgProd,AvgProdAllAP(1).All95Conf,'Color',Colors(1).Color,'LineWidth',2.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(2).AvgProd, AvgProdAllAP(2).All95Conf,'Color',Colors(2).Color,'LineWidth',2.5);
 errorbar(EggLength,AvgProdAllAP(9).AvgProd, AvgProdAllAP(9).All95Conf,'Color',Colors(9).Color,'LineWidth',2.5);
 yyaxis right
 plot(EggLength, (AvgProdAllAP(1).AvgProd./MS2Conversion), 'Color', Colors(1).Color);
 plot(EggLength, (AvgProdAllAP(2).AvgProd./MS2Conversion), 'Color', Colors(2).Color);
 plot(EggLength, (AvgProdAllAP(9).AvgProd./MS2Conversion), 'Color', Colors(9).Color);
 set(gca,'Box','off', 'FontSize', fontsize, 'FontName', fontname,'YColor','k');
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 x_width y_width]);
 xlabel('% egg length')
 ylabel('total transcripts produced')
 %title('Avg mRNA production per nucleus') 
 ylim([0 900000./MS2Conversion]);
 xlim([0 100]);
 yyaxis left 
 ylabel('integrated fluorescence (AU)');
 ylim([0 900000]);
print( [FigDirect filesep 'SinglesvBothTotalmRNA'],'-dsvg');
 
  figure
 errorbar(EggLength,AvgProdAllAP(1).AvgProd,AvgProdAllAP(1).All95Conf,'Color',Colors(1).Color,'LineWidth',2.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(2).AvgProd, AvgProdAllAP(2).All95Conf,'Color',Colors(2).Color,'LineWidth',2.5);
 %legend('Distal','Proximal','Both Sep');
 yyaxis right
 plot(EggLength, (AvgProdAllAP(1).AvgProd./MS2Conversion), 'Color', Colors(1).Color);
 plot(EggLength, (AvgProdAllAP(2).AvgProd./MS2Conversion), 'Color', Colors(2).Color);
 ylabel('total transcripts produced');
 set(gca, 'FontSize', fontsize, 'FontName', fontname,'YColor','k');
 ylim([0 900000./MS2Conversion]);
 xlim([0 100]);
 yyaxis left
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 %title('Avg mRNA production per nucleus') 
 ylim([0 900000]);
 print( '-painters',[FigDirect filesep 'DistProxTotalmRNA'],'-dsvg');

 %% single allele 
  figure
 errorbar(EggLength,AvgProdAllAP(3).AvgProd,AvgProdAllAP(3).All95Conf,'Color',Colors(3).Color,'LineWidth',2.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(4).AvgProd, AvgProdAllAP(4).All95Conf,'Color',Colors(4).Color,'LineWidth',2.5,'LineStyle',':');
  errorbar(EggLength,AvgProdAllAP(5).AvgProd, AvgProdAllAP(5).All95Conf,'Color',Colors(5).Color,'LineWidth',2.5,'LineStyle',':');
 legend('Both Sep','1x Distal','1x Proximal');
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 %title('Avg mRNA production per nucleus') 
 ylim([0 1400000]);
 xlim([0 100]);
 
 figure
 errorbar(EggLength,AvgProdAllAP(1).AvgProd,AvgProdAllAP(1).All95Conf,'Color',Colors(1).Color,'LineWidth',2.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(4).AvgProd, AvgProdAllAP(4).All95Conf,'Color',Colors(4).Color,'LineWidth',2.5,'LineStyle',':');
 %legend('Distal','1x Distal');
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 %title('Avg mRNA production per nucleus') 
 ylim([0 1400000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(1).AvgProd./MS2Conversion),(AvgProdAllAP(1).All95Conf./MS2Conversion),'Color',Colors(1).Color,'LineWidth',2.5,'LineStyle','-');
  errorbar(EggLength, (AvgProdAllAP(4).AvgProd./MS2Conversion),(AvgProdAllAP(4).All95Conf./MS2Conversion),'Color',Colors(4).Color,'LineWidth',2.5,'LineStyle',':');
ylabel('total transcripts produced');
ylim([0 (1400000/MS2Conversion)]);
set(gca, 'YColor','k');
HemiPct = (AvgProdAllAP(1).AvgProd(ExpressionBins{1,1}) / AvgProdAllAP(4).AvgProd(ExpressionBins{1,4})); %Homo/Hemi
title([num2str(HemiPct),'% of hemi']);
print( '-painters',[FigDirect filesep 'HemiDistTotalmRNA'],'-dpdf');
 
  figure
 errorbar(EggLength,AvgProdAllAP(2).AvgProd,AvgProdAllAP(2).All95Conf,'Color',Colors(2).Color,'LineWidth',2.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(5).AvgProd, AvgProdAllAP(5).All95Conf,'Color',Colors(5).Color,'LineWidth',2.5,'LineStyle',':');
 %legend('Proximal','1x Proximal');
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 %title('Avg mRNA production per nucleus') 
 ylim([0 1400000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(2).AvgProd./MS2Conversion),(AvgProdAllAP(2).All95Conf./MS2Conversion),'Color',Colors(2).Color,'LineWidth',2.5,'LineStyle','-');
  errorbar(EggLength, (AvgProdAllAP(5).AvgProd./MS2Conversion),(AvgProdAllAP(5).All95Conf./MS2Conversion),'Color',Colors(5).Color,'LineWidth',2.5,'LineStyle',':');
ylabel('total transcripts produced');
ylim([0 (1400000/MS2Conversion)]);
set(gca,'YColor','k');
HemiPct = (AvgProdAllAP(2).AvgProd(ExpressionBins{1,2}) / AvgProdAllAP(5).AvgProd(ExpressionBins{1,5})); %Homo/Hemi
title([num2str(HemiPct),'% of hemi']);
print( '-painters',[FigDirect filesep 'HemiProxTotalmRNA'],'-dpdf');

figure
 errorbar(EggLength,AvgProdAllAP(6).AvgProd,AvgProdAllAP(6).All95Conf,'Color',Colors(6).Color,'LineWidth',2.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(26).AvgProd, AvgProdAllAP(26).All95Conf,'Color',Colors(26).Color,'LineWidth',2.5,'LineStyle',':');
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 %title('Avg mRNA production per nucleus') 
 ylim([0 1400000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(6).AvgProd./MS2Conversion),(AvgProdAllAP(6).All95Conf./MS2Conversion),'Color',Colors(6).Color,'LineWidth',2.5,'LineStyle','-');
  errorbar(EggLength, (AvgProdAllAP(26).AvgProd./MS2Conversion),(AvgProdAllAP(26).All95Conf./MS2Conversion),'Color',Colors(26).Color,'LineWidth',2.5,'LineStyle',':');
ylabel('total transcripts produced');
ylim([0 (1400000/MS2Conversion)]);
set(gca, 'YColor','k');
print( '-painters',[FigDirect filesep 'Hemi2xDistTotalmRNA'],'-dpdf');

 figure
 errorbar(EggLength,AvgProdAllAP(9).AvgProd,AvgProdAllAP(9).All95Conf,'Color',Colors(9).Color,'LineWidth',2.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(10).AvgProd, AvgProdAllAP(10).All95Conf,'Color',Colors(10).Color,'LineWidth',2.5,'LineStyle',':');
 %legend('Both','1x Both');
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 %title('Avg mRNA production per nucleus') 
 ylim([0 1400000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(9).AvgProd./MS2Conversion),(AvgProdAllAP(9).All95Conf./MS2Conversion),'Color',Colors(9).Color,'LineWidth',2.5,'LineStyle','-');
  errorbar(EggLength, (AvgProdAllAP(10).AvgProd./MS2Conversion),(AvgProdAllAP(10).All95Conf./MS2Conversion),'Color',Colors(10).Color,'LineWidth',2.5,'LineStyle',':');
ylabel('total transcripts produced');
ylim([0 (1400000/MS2Conversion)]);
set(gca,'YColor','k');
print( '-painters',[FigDirect filesep 'HemiBothTotalmRNA'],'-dpdf');

figure
 errorbar(EggLength,AvgProdAllAP(7).AvgProd,AvgProdAllAP(7).All95Conf,'Color',Colors(7).Color,'LineWidth',2.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(8).AvgProd,AvgProdAllAP(8).All95Conf,'Color',Colors(8).Color,'LineWidth',2.5,'LineStyle',':');
 %legend('2xProximal','Single 2xProximal');
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 ylim([0 1400000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(7).AvgProd./MS2Conversion),(AvgProdAllAP(7).All95Conf./MS2Conversion),'Color',Colors(7).Color,'LineWidth',2.5,'LineStyle','-');
  errorbar(EggLength, (AvgProdAllAP(8).AvgProd./MS2Conversion),(AvgProdAllAP(8).All95Conf./MS2Conversion),'Color',Colors(8).Color,'LineWidth',2.5,'LineStyle',':');
ylabel('total transcripts produced');
ylim([0 (1400000/MS2Conversion)]);
set(gca,'YColor','k');
print( '-painters',[FigDirect filesep 'Hemi2xProxTotalmRNA'],'-dpdf');

%% No promoter/MS2 constructs
figure 
errorbar(EggLength, AvgProdAllAP(1).AvgProd, AvgProdAllAP(1).All95Conf, 'Color', Colors(1).Color, 'LineWidth', 3.5);
hold on 
errorbar(EggLength, AvgProdAllAP(4).AvgProd, AvgProdAllAP(4).All95Conf, 'Color', Colors(1).Color, 'LineStyle','--', 'LineWidth',3.5);
errorbar(EggLength, AvgProdAllAP(63).AvgProd, AvgProdAllAP(63).All95Conf, 'Color', 'k', 'LineWidth', 3.5);
ylim([0 1500000]);
xlim([20 80]);
yyaxis right 
errorbar(EggLength, (AvgProdAllAP(1).AvgProd./MS2Conversion), (AvgProdAllAP(1).All95Conf./MS2Conversion), 'Color', Colors(1).Color, 'LineWidth', 3.5);
hold on 
errorbar(EggLength, (AvgProdAllAP(4).AvgProd./MS2Conversion), (AvgProdAllAP(4).All95Conf./MS2Conversion), 'Color', Colors(1).Color, 'LineStyle','--', 'LineWidth',3.5);
errorbar(EggLength, (AvgProdAllAP(63).AvgProd./MS2Conversion), (AvgProdAllAP(63).All95Conf./MS2Conversion), 'Color', 'k', 'LineWidth', 3.5);
ylim([0 1500000/MS2Conversion]);
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'Dist_noProm_APLine','.pdf'],'pdf');


% Compare to hemi vs homozygotes at peak expression bins 
figure 
x_sizeNarrow =1.3;
y_sizeNarrow =3.97;
figure
lineProps.col{1} = [Colors(1).Color];
mseb([-9:10], [ones(1,20).*AvgProdAllAP(4).AvgProd(ExpressionBins{1,4})], ones(1,20).*AvgProdAllAP(4).All95Conf(ExpressionBins{1,4}), lineProps,0)
hold on
mseb([-9:10], [ones(1,20).*AvgProdAllAP(1).AvgProd(ExpressionBins{1,1})], ones(1,20).*AvgProdAllAP(1).All95Conf(ExpressionBins{1,1}), lineProps,0)
%h2 = line([0,10],[AvgProdAllAP(6).AvgProd(ExpressionBins{1,6}),AvgProdAllAP(6).AvgProd(ExpressionBins{1,6})]);
%h2.Color = Colors(6).Color; h2.LineStyle ='--'; h2.LineWidth =2.5;
errorbar(0.2, AvgProdAllAP(63).AvgProd(ExpressionBins{1,63}),AvgProdAllAP(63).All95Conf(ExpressionBins{1,63}),'o','Color','k','LineWidth',2.5);
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop x_sizeNarrow y_sizeNarrow]);
ylabel('integrated fluorescence (AU)');
xticks([0.2]);
xticklabels('No prom/MS2');
xlim([0 0.4]);
ylim([3e5, 1.5e6]);
saveas(gcf, [FigDirect filesep 'Distal_CompNoProm_Shade','.pdf'],'pdf');
% calc pval of noprom vs hemizygotes
[h,p] = kstest2([AvgProdAllAP(63).AllProds(:,ExpressionBins{1,63})],[AvgProdAllAP(4).AllProds(:,ExpressionBins{1,4})]);

% also do KW test - need in one matrix so make same length
KWTestMat = nan(300,2);
KWTestMat([1:length(AvgProdAllAP(63).AllProds(:,ExpressionBins{1,63}))],1) = [AvgProdAllAP(63).AllProds(:,ExpressionBins{1,63})];
KWTestMat([1:length(AvgProdAllAP(4).AllProds(:,ExpressionBins{1,4}))],2) = [AvgProdAllAP(4).AllProds(:,ExpressionBins{1,4})];
[p,tbl,stats] = kruskalwallis(KWTestMat);
multcompare(stats);

%% Compare old vs new MS2 conversion 
cc = 56; 
%cc2 = 4;
CutOffDate = input('which embryo is first of new data?');
NewAvg = []; NewAvgSE=[];
OldAvg=[]; OldAvgSE=[];
for ee = 1:CutOffDate-1
    OldAvg = [OldAvg; AvgProdAllAP(cc).EmbryosProd(ee).MeanProd];
    OldAvgSE = [OldAvgSE; AvgProdAllAP(cc).EmbryosProd(ee).SE];
end
for ee = CutOffDate:length(AvgProdAllAP(cc).EmbryosProd)
    NewAvg=[NewAvg; AvgProdAllAP(cc).EmbryosProd(ee).MeanProd];
    NewAvgSE = [NewAvgSE; AvgProdAllAP(cc).EmbryosProd(ee).SE];
end

figure
errorbar(EggLength, nanmean(OldAvg), nanmean(OldAvgSE), 'LineWidth', 2.5, 'Color', 'k');
hold on 
errorbar(EggLength, nanmean(NewAvg), nanmean(NewAvgSE), 'LineWidth', 2.5, 'Color', 'r');
xlim([20 80]);
ylim([0 1500000]);
ylabel('integrated fluorescence (AU)');
xlabel('% egg length');
legend({'old data','new data'});
yyaxis right 
errorbar(EggLength, (nanmean(OldAvg)./MS2Conversion), (nanmean(OldAvgSE)./MS2Conversion), 'LineWidth', 2.5, 'Color', 'k');
hold on 
errorbar(EggLength, (nanmean(NewAvg)./MS2Conversion), (nanmean(NewAvgSE)./MS2Conversion), 'LineWidth', 2.5, 'Color', 'r');
ylim([0 1500000/MS2Conversion]);
ylabel('total transcripts produced');
title('old MS2 conversion');

% Try w new MS2 conversion on the new data 
figure 
errorbar(EggLength, (nanmean(OldAvg)./MS2Conversion), (nanmean(OldAvgSE)./MS2Conversion), 'LineWidth', 2.5, 'Color', 'k');
hold on 
errorbar(EggLength, (nanmean(NewAvg)./MS2Conversion_New), (nanmean(NewAvgSE)./MS2Conversion_New), 'LineWidth', 2.5, 'Color', 'r');
ylabel('total transcripts produced');
title('old MS2 conversion');
xlim([20 80]);
ylim([0 1500000/MS2Conversion]);
title('new MS2 conversion on new data');
legend({'old data','new data'});



%% Hemi vs homo chr 2 vs 3 
% Duplicated distal 
figure
% Chr 2
XError = AvgProdAllAP(26).All95Conf(ExpressionBins{1,26});
YError = AvgProdAllAP(6).All95Conf(ExpressionBins{1,6});    
errorbar(AvgProdAllAP(26).AvgProd(ExpressionBins{1,26}), AvgProdAllAP(6).AvgProd(ExpressionBins{1,6}),YError, YError, XError, XError,'o','Color',Colors(6).Color,'LineWidth',2.5);
%plot(AvgProdAllAP(26).AvgProd(ExpressionBins{1,26}), AvgProdAllAP(6).AvgProd(ExpressionBins{1,6}),'o', 'Color',Colors(6).Color, 'LineWidth',2.5);
hold on 
% Chr 3
XError = AvgProdAllAP(38).All95Conf(ExpressionBins{1,38});
YError = AvgProdAllAP(36).All95Conf(ExpressionBins{1,36});    
b = errorbar(AvgProdAllAP(38).AvgProd(ExpressionBins{1,38}), AvgProdAllAP(36).AvgProd(ExpressionBins{1,36}),YError, YError, XError, XError,'o','Color',Colors(6).Color,'LineWidth',2.5);
b.Bar.LineStyle = 'dotted'

xlim([0.7e6 1.4e6]);
ylim([0.7e6 1.4e6]); 
hline = refline(1,0);
hline.LineWidth =2.5;
hline.Color ='k'
% add dotted lines of expected homozygous values
    %Chr 2
h = line([AvgProdAllAP(26).AvgProd(ExpressionBins{1,26}), AvgProdAllAP(26).AvgProd(ExpressionBins{1,26})],[0.7e6, AvgProdAllAP(26).AvgProd(ExpressionBins{1,26})]);
h.LineStyle ='--'; h.Color = 'k'; h.LineWidth=2;
h2 = line([0.7e6, AvgProdAllAP(26).AvgProd(ExpressionBins{1,26})],[AvgProdAllAP(26).AvgProd(ExpressionBins{1,26}), AvgProdAllAP(26).AvgProd(ExpressionBins{1,26})]);
h2.Color ='k'; h2.LineStyle ='--'; h2.LineWidth =2;
    %Chr 3
h = line([AvgProdAllAP(38).AvgProd(ExpressionBins{1,38}), AvgProdAllAP(38).AvgProd(ExpressionBins{1,38})],[0.7e6, AvgProdAllAP(38).AvgProd(ExpressionBins{1,38})]);
h.LineStyle ='--'; h.Color = 'k'; h.LineWidth=2;
h2 = line([0.7e6, AvgProdAllAP(38).AvgProd(ExpressionBins{1,38})],[AvgProdAllAP(38).AvgProd(ExpressionBins{1,38}), AvgProdAllAP(38).AvgProd(ExpressionBins{1,38})]);
h2.Color ='k'; h2.LineStyle ='--'; h2.LineWidth =2;

ylabel('homozygous average total mRNA per allele');
xlabel('hemizygous average total mRNA per allele');
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'DupDist_HemiHomo_Chr2_3','.pdf'],'pdf');

% Bar graph for inset in figure
figure 
Chr2_RelExp = AvgProdAllAP(6).AvgProd(ExpressionBins{1,6}) / AvgProdAllAP(26).AvgProd(ExpressionBins{1,26});
Chr3_RelExp = AvgProdAllAP(36).AvgProd(ExpressionBins{1,36}) / AvgProdAllAP(38).AvgProd(ExpressionBins{1,38});
RelExpLvls = [Chr2_RelExp*100, Chr3_RelExp*100];
bar(RelExpLvls,'FaceColor', Colors(6).Color);
hold on 
 %propogate dat error
Chr2FractError = (AvgProdAllAP(6).All95Conf(ExpressionBins{1,6}) / AvgProdAllAP(6).AvgProd(ExpressionBins{1,6}))^2;
Chr2FractError = Chr2FractError + (AvgProdAllAP(26).All95Conf(ExpressionBins{1,26}) / AvgProdAllAP(26).AvgProd(ExpressionBins{1,26}))^2;
Chr2FractError = sqrt(Chr2FractError) * RelExpLvls(1); %mult and division gives relative std

Chr3FractError = (AvgProdAllAP(36).All95Conf(ExpressionBins{1,36}) / AvgProdAllAP(36).AvgProd(ExpressionBins{1,36}))^2;
Chr3FractError = Chr3FractError + (AvgProdAllAP(38).All95Conf(ExpressionBins{1,38}) / AvgProdAllAP(38).AvgProd(ExpressionBins{1,38}))^2;
Chr3FractError = sqrt(Chr3FractError) * RelExpLvls(2);
errorbar(1, RelExpLvls(1), Chr2FractError, 'Color', 'k', 'LineWidth', 2.5);
errorbar(2, RelExpLvls(2), Chr3FractError, 'Color', 'k', 'LineWidth', 2.5);
ylabel('% homozygous expression');
xticks([1,2]);
xticklabels({'Chr2', 'Chr3'});
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'DupDist_RelHomoLvls_Chr2v3','.svg'],'svg');

% shadow pair 
figure
% Chr 2
XError = AvgProdAllAP(10).All95Conf(ExpressionBins{1,10});
YError = AvgProdAllAP(9).All95Conf(ExpressionBins{1,9});    
errorbar(AvgProdAllAP(10).AvgProd(ExpressionBins{1,10}), AvgProdAllAP(9).AvgProd(ExpressionBins{1,9}),YError, YError, XError, XError,'o','Color',Colors(9).Color,'LineWidth',2.5);
%plot(AvgProdAllAP(26).AvgProd(ExpressionBins{1,26}), AvgProdAllAP(6).AvgProd(ExpressionBins{1,6}),'o', 'Color',Colors(6).Color, 'LineWidth',2.5);
hold on 
% Chr 3
XError = AvgProdAllAP(40).All95Conf(ExpressionBins{1,40});
YError = AvgProdAllAP(35).All95Conf(ExpressionBins{1,35});    
b = errorbar(AvgProdAllAP(40).AvgProd(ExpressionBins{1,40}), AvgProdAllAP(35).AvgProd(ExpressionBins{1,35}),YError, YError, XError, XError,'o','Color',Colors(9).Color,'LineWidth',2.5);
b.Bar.LineStyle = 'dotted'

xlim([0.7e6 1.6e6]);
ylim([0.7e6 1.6e6]); 
hline = refline(1,0);
hline.LineWidth =2.5;
hline.Color ='k'
ylabel('homozygous average total mRNA per allele');
xlabel('hemizygous average total mRNA per allele');
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'SE_HemiHomo_Chr2_3','.pdf'],'pdf');

% Bar graph for inset in figure
figure 
Chr2_RelExp = AvgProdAllAP(9).AvgProd(ExpressionBins{1,9}) / AvgProdAllAP(10).AvgProd(ExpressionBins{1,10});
Chr3_RelExp = AvgProdAllAP(35).AvgProd(ExpressionBins{1,35}) / AvgProdAllAP(40).AvgProd(ExpressionBins{1,40});
RelExpLvls = [Chr2_RelExp*100, Chr3_RelExp*100];
bar(RelExpLvls,'FaceColor', Colors(9).Color);
hold on 
 %propogate dat error
Chr2FractError = (AvgProdAllAP(9).All95Conf(ExpressionBins{1,9}) / AvgProdAllAP(9).AvgProd(ExpressionBins{1,9}))^2;
Chr2FractError = Chr2FractError + (AvgProdAllAP(10).All95Conf(ExpressionBins{1,10}) / AvgProdAllAP(10).AvgProd(ExpressionBins{1,10}))^2;
Chr2FractError = sqrt(Chr2FractError) * RelExpLvls(1); %mult and division gives relative std

Chr3FractError = (AvgProdAllAP(35).All95Conf(ExpressionBins{1,35}) / AvgProdAllAP(35).AvgProd(ExpressionBins{1,35}))^2;
Chr3FractError = Chr3FractError + (AvgProdAllAP(40).All95Conf(ExpressionBins{1,40}) / AvgProdAllAP(40).AvgProd(ExpressionBins{1,40}))^2;
Chr3FractError = sqrt(Chr3FractError) * RelExpLvls(2);
errorbar(1, RelExpLvls(1), Chr2FractError, 'Color', 'k', 'LineWidth', 2.5);
errorbar(2, RelExpLvls(2), Chr3FractError, 'Color', 'k', 'LineWidth', 2.5);
ylabel('% homozygous expression');
xticks([1,2]);
xticklabels({'Chr2', 'Chr3'});
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'SE_RelHomoLvls_Chr2v3','.svg'],'svg');
%% Compare hemizygotes to TFBS array hetreozygotes
figure
 errorbar(EggLength,AvgProdAllAP(6).AvgProd,AvgProdAllAP(6).All95Conf,'Color',Colors(6).Color,'LineWidth',3.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(26).AvgProd, AvgProdAllAP(26).All95Conf,'Color',Colors(26).Color,'LineWidth',3.5,'LineStyle',':');
 errorbar(EggLength, AvgProdAllAP(39).AvgProd, AvgProdAllAP(39).All95Conf, 'Color', 'k', 'LineWidth',3.5,'LineStyle',"--");
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 %title('Avg mRNA production per nucleus') 
 ylim([0 1400000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(6).AvgProd./MS2Conversion),(AvgProdAllAP(6).All95Conf./MS2Conversion),'Color',Colors(6).Color,'LineWidth',3.5,'LineStyle','-');
errorbar(EggLength, (AvgProdAllAP(26).AvgProd./MS2Conversion),(AvgProdAllAP(26).All95Conf./MS2Conversion),'Color',Colors(26).Color,'LineWidth',3.5,'LineStyle',':');
errorbar(EggLength, (AvgProdAllAP(39).AvgProd./MS2Conversion),(AvgProdAllAP(39).All95Conf./MS2Conversion),'Color','k','LineWidth',3.5,'LineStyle','--');
TFPct = (AvgProdAllAP(39).AvgProd(ExpressionBins{1,39}) / AvgProdAllAP(26).AvgProd(ExpressionBins{1,26})); %Homo/Hemi
ylabel('total transcripts produced');
ylim([0 (1400000/MS2Conversion)]);
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
% % higher hemi is than homo 
TFInc = (AvgProdAllAP(26).AvgProd(ExpressionBins{1,26}) - AvgProdAllAP(39).AvgProd(ExpressionBins{1,39}))/(AvgProdAllAP(26).AvgProd(ExpressionBins{1,26}));
title([num2str(HemiPct*100),'% of hemi', num2str(TFInc),'%higher than TF']);
saveas(gcf, [FigDirect filesep 'TotProdKr2xDist_1xBcd','.pdf'],'pdf');
%print( '-painters',[FigDirect filesep 'Hemi2xDistTotalmRNA'],'-dpdf');

% 2xDist w 1xHb array
figure
 errorbar(EggLength,AvgProdAllAP(6).AvgProd,AvgProdAllAP(6).All95Conf,'Color',Colors(6).Color,'LineWidth',3.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(26).AvgProd, AvgProdAllAP(26).All95Conf,'Color',Colors(26).Color,'LineWidth',3.5,'LineStyle',':');
 errorbar(EggLength, AvgProdAllAP(42).AvgProd, AvgProdAllAP(42).All95Conf, 'Color', 'k', 'LineWidth',3.5,'LineStyle',"--");
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 %title('Avg mRNA production per nucleus') 
 ylim([0 1400000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(6).AvgProd./MS2Conversion),(AvgProdAllAP(6).All95Conf./MS2Conversion),'Color',Colors(6).Color,'LineWidth',3.5,'LineStyle','-');
errorbar(EggLength, (AvgProdAllAP(26).AvgProd./MS2Conversion),(AvgProdAllAP(26).All95Conf./MS2Conversion),'Color',Colors(26).Color,'LineWidth',3.5,'LineStyle',':');
errorbar(EggLength, (AvgProdAllAP(42).AvgProd./MS2Conversion),(AvgProdAllAP(42).All95Conf./MS2Conversion),'Color','k','LineWidth',3.5,'LineStyle','--');
ylabel('total transcripts produced');
ylim([0 (1400000/MS2Conversion)]);
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'TotProdKr2xDist_1xHb','.pdf'],'pdf');

%2xDist w 1xZld array
figure
 errorbar(EggLength,AvgProdAllAP(6).AvgProd,AvgProdAllAP(6).All95Conf,'Color',Colors(6).Color,'LineWidth',3.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(26).AvgProd, AvgProdAllAP(26).All95Conf,'Color',Colors(26).Color,'LineWidth',3.5,'LineStyle',':');
 errorbar(EggLength, AvgProdAllAP(45).AvgProd, AvgProdAllAP(45).All95Conf, 'Color', 'k', 'LineWidth',3.5,'LineStyle',"--");
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 %title('Avg mRNA production per nucleus') 
 ylim([0 1400000]);
 xlim([20 80]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(6).AvgProd./MS2Conversion),(AvgProdAllAP(6).All95Conf./MS2Conversion),'Color',Colors(6).Color,'LineWidth',3.5,'LineStyle','-');
errorbar(EggLength, (AvgProdAllAP(26).AvgProd./MS2Conversion),(AvgProdAllAP(26).All95Conf./MS2Conversion),'Color',Colors(26).Color,'LineWidth',3.5,'LineStyle',':');
errorbar(EggLength, (AvgProdAllAP(45).AvgProd./MS2Conversion),(AvgProdAllAP(45).All95Conf./MS2Conversion),'Color','k','LineWidth',3.5,'LineStyle','--');
ylabel('total transcripts produced');
ylim([0 (1400000/MS2Conversion)]);
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'TotProdKr2xDist_1xZld','.pdf'],'pdf');

% 2xDist w Twi array
figure
 errorbar(EggLength,AvgProdAllAP(6).AvgProd,AvgProdAllAP(6).All95Conf,'Color',Colors(6).Color,'LineWidth',3.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(26).AvgProd, AvgProdAllAP(26).All95Conf,'Color',Colors(26).Color,'LineWidth',3.5,'LineStyle',':');
 errorbar(EggLength, AvgProdAllAP(58).AvgProd, AvgProdAllAP(58).All95Conf, 'Color', 'k', 'LineWidth',3.5,'LineStyle',"--");
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 %title('Avg mRNA production per nucleus') 
 ylim([0 1400000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(6).AvgProd./MS2Conversion),(AvgProdAllAP(6).All95Conf./MS2Conversion),'Color',Colors(6).Color,'LineWidth',3.5,'LineStyle','-');
errorbar(EggLength, (AvgProdAllAP(26).AvgProd./MS2Conversion),(AvgProdAllAP(26).All95Conf./MS2Conversion),'Color',Colors(26).Color,'LineWidth',3.5,'LineStyle',':');
errorbar(EggLength, (AvgProdAllAP(58).AvgProd./MS2Conversion),(AvgProdAllAP(58).All95Conf./MS2Conversion),'Color','k','LineWidth',3.5,'LineStyle','--');
ylabel('total transcripts produced');
ylim([0 (1400000/MS2Conversion)]);
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
%saveas(gcf, [FigDirect filesep 'TotProdKr2xDist_1xZld','.pdf'],'pdf');

% SE w 1xBcd array
figure
 errorbar(EggLength,AvgProdAllAP(9).AvgProd,AvgProdAllAP(9).All95Conf,'Color',Colors(9).Color,'LineWidth',3.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(10).AvgProd, AvgProdAllAP(10).All95Conf,'Color',Colors(10).Color,'LineWidth',3.5,'LineStyle',':');
 errorbar(EggLength, AvgProdAllAP(41).AvgProd, AvgProdAllAP(41).All95Conf, 'Color', 'k', 'LineWidth',3.5,'LineStyle',"--");
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 %title('Avg mRNA production per nucleus') 
 ylim([0 1400000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(9).AvgProd./MS2Conversion),(AvgProdAllAP(9).All95Conf./MS2Conversion),'Color',Colors(9).Color,'LineWidth',3.5,'LineStyle','-');
errorbar(EggLength, (AvgProdAllAP(10).AvgProd./MS2Conversion),(AvgProdAllAP(10).All95Conf./MS2Conversion),'Color',Colors(10).Color,'LineWidth',3.5,'LineStyle',':');
errorbar(EggLength, (AvgProdAllAP(41).AvgProd./MS2Conversion),(AvgProdAllAP(41).All95Conf./MS2Conversion),'Color','k','LineWidth',3.5,'LineStyle','--');
ylabel('total transcripts produced');
ylim([0 (1400000/MS2Conversion)]);
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
TFPct = (AvgProdAllAP(41).AvgProd(ExpressionBins{1,41}) / AvgProdAllAP(10).AvgProd(ExpressionBins{1,10})); %Homo/Hemi
title([num2str(HemiPct),'% of hemi']);
saveas(gcf, [FigDirect filesep 'TotProdKrSE_1xBcd','.pdf'],'pdf');

% Both w 1xHb array

figure
 errorbar(EggLength,AvgProdAllAP(9).AvgProd,AvgProdAllAP(9).All95Conf,'Color',Colors(9).Color,'LineWidth',3.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(10).AvgProd, AvgProdAllAP(10).All95Conf,'Color',Colors(10).Color,'LineWidth',3.5,'LineStyle',':');
 errorbar(EggLength, AvgProdAllAP(44).AvgProd, AvgProdAllAP(44).All95Conf, 'Color', 'k', 'LineWidth',3.5,'LineStyle',"--");
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 %title('Avg mRNA production per nucleus') 
 ylim([0 1400000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(9).AvgProd./MS2Conversion),(AvgProdAllAP(9).All95Conf./MS2Conversion),'Color',Colors(9).Color,'LineWidth',3.5,'LineStyle','-');
errorbar(EggLength, (AvgProdAllAP(10).AvgProd./MS2Conversion),(AvgProdAllAP(10).All95Conf./MS2Conversion),'Color',Colors(10).Color,'LineWidth',3.5,'LineStyle',':');
errorbar(EggLength, (AvgProdAllAP(44).AvgProd./MS2Conversion),(AvgProdAllAP(44).All95Conf./MS2Conversion),'Color','k','LineWidth',3.5,'LineStyle','--');
ylabel('total transcripts produced');
ylim([0 (1400000/MS2Conversion)]);
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'TotProdKrSE_1xHb','.pdf'],'pdf');

% Both w 1xZld array

figure
 errorbar(EggLength,AvgProdAllAP(9).AvgProd,AvgProdAllAP(9).All95Conf,'Color',Colors(9).Color,'LineWidth',3.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(10).AvgProd, AvgProdAllAP(10).All95Conf,'Color',Colors(10).Color,'LineWidth',3.5,'LineStyle',':');
 errorbar(EggLength, AvgProdAllAP(47).AvgProd, AvgProdAllAP(47).All95Conf, 'Color', 'k', 'LineWidth',3.5,'LineStyle',"--");
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 %title('Avg mRNA production per nucleus') 
 ylim([0 1400000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(9).AvgProd./MS2Conversion),(AvgProdAllAP(9).All95Conf./MS2Conversion),'Color',Colors(9).Color,'LineWidth',3.5,'LineStyle','-');
errorbar(EggLength, (AvgProdAllAP(10).AvgProd./MS2Conversion),(AvgProdAllAP(10).All95Conf./MS2Conversion),'Color',Colors(10).Color,'LineWidth',3.5,'LineStyle',':');
errorbar(EggLength, (AvgProdAllAP(47).AvgProd./MS2Conversion),(AvgProdAllAP(47).All95Conf./MS2Conversion),'Color','k','LineWidth',3.5,'LineStyle','--');
ylabel('total transcripts produced');
ylim([0 (1400000/MS2Conversion)]);
xlim([20 80]);
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'TotProdKrSE_1xZld','.pdf'],'pdf');

% Both w 1xStat92E array

figure
 errorbar(EggLength,AvgProdAllAP(9).AvgProd,AvgProdAllAP(9).All95Conf,'Color',Colors(9).Color,'LineWidth',3.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(10).AvgProd, AvgProdAllAP(10).All95Conf,'Color',Colors(10).Color,'LineWidth',3.5,'LineStyle',':');
 errorbar(EggLength, AvgProdAllAP(48).AvgProd, AvgProdAllAP(48).All95Conf, 'Color', 'k', 'LineWidth',3.5,'LineStyle',"--");
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 %title('Avg mRNA production per nucleus') 
 ylim([0 1400000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(9).AvgProd./MS2Conversion),(AvgProdAllAP(9).All95Conf./MS2Conversion),'Color',Colors(9).Color,'LineWidth',3.5,'LineStyle','-');
errorbar(EggLength, (AvgProdAllAP(10).AvgProd./MS2Conversion),(AvgProdAllAP(10).All95Conf./MS2Conversion),'Color',Colors(10).Color,'LineWidth',3.5,'LineStyle',':');
errorbar(EggLength, (AvgProdAllAP(48).AvgProd./MS2Conversion),(AvgProdAllAP(48).All95Conf./MS2Conversion),'Color','k','LineWidth',3.5,'LineStyle','--');
ylabel('total transcripts produced');
ylim([0 (1400000/MS2Conversion)]);
xlim([20 80])
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'TotProdKrSE_1xStat92E','.pdf'],'pdf');

% 2xProx w 1xHb array
figure
 errorbar(EggLength,AvgProdAllAP(7).AvgProd,AvgProdAllAP(7).All95Conf,'Color',Colors(7).Color,'LineWidth',3.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(8).AvgProd, AvgProdAllAP(8).All95Conf,'Color',Colors(8).Color,'LineWidth',3.5,'LineStyle',':');
 errorbar(EggLength, AvgProdAllAP(43).AvgProd, AvgProdAllAP(43).All95Conf, 'Color', 'k', 'LineWidth',3.5,'LineStyle',"--");
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 %title('Avg mRNA production per nucleus') 
 ylim([0 1000000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(7).AvgProd./MS2Conversion),(AvgProdAllAP(7).All95Conf./MS2Conversion),'Color',Colors(7).Color,'LineWidth',3.5,'LineStyle','-');
errorbar(EggLength, (AvgProdAllAP(8).AvgProd./MS2Conversion),(AvgProdAllAP(8).All95Conf./MS2Conversion),'Color',Colors(8).Color,'LineWidth',3.5,'LineStyle',':');
errorbar(EggLength, (AvgProdAllAP(43).AvgProd./MS2Conversion),(AvgProdAllAP(43).All95Conf./MS2Conversion),'Color','k','LineWidth',3.5,'LineStyle','--');
ylabel('total transcripts produced');
ylim([0 (1000000/MS2Conversion)]);
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
%Calculate % of hemi expression happens in homo 
HemiPct = (AvgProdAllAP(8).AvgProd(ExpressionBins{1,8}) / AvgProdAllAP(7).AvgProd(ExpressionBins{1,8})); %Homo/Hemi
title([num2str(HemiPct),'% of hemi']);
saveas(gcf, [FigDirect filesep 'TotProdKr2xProx_1xHb','.pdf'],'pdf');

% 2xProx w 1xStat92E array
figure
 errorbar(EggLength,AvgProdAllAP(7).AvgProd,AvgProdAllAP(7).All95Conf,'Color',Colors(7).Color,'LineWidth',3.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(8).AvgProd, AvgProdAllAP(8).All95Conf,'Color',Colors(8).Color,'LineWidth',3.5,'LineStyle',':');
 errorbar(EggLength, AvgProdAllAP(46).AvgProd, AvgProdAllAP(46).All95Conf, 'Color', 'k', 'LineWidth',3.5,'LineStyle',"--");
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 %title('Avg mRNA production per nucleus') 
 ylim([0 1000000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(7).AvgProd./MS2Conversion),(AvgProdAllAP(7).All95Conf./MS2Conversion),'Color',Colors(7).Color,'LineWidth',3.5,'LineStyle','-');
errorbar(EggLength, (AvgProdAllAP(8).AvgProd./MS2Conversion),(AvgProdAllAP(8).All95Conf./MS2Conversion),'Color',Colors(8).Color,'LineWidth',3.5,'LineStyle',':');
errorbar(EggLength, (AvgProdAllAP(46).AvgProd./MS2Conversion),(AvgProdAllAP(46).All95Conf./MS2Conversion),'Color','k','LineWidth',3.5,'LineStyle','--');
ylabel('total transcripts produced');
ylim([0 (1000000/MS2Conversion)]);
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'TotProdKr2xProx_1xStat92E','.pdf'],'pdf');

% 2xProx w 1xBcd array
figure
 errorbar(EggLength,AvgProdAllAP(7).AvgProd,AvgProdAllAP(7).All95Conf,'Color',Colors(7).Color,'LineWidth',3.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(8).AvgProd, AvgProdAllAP(8).All95Conf,'Color',Colors(8).Color,'LineWidth',3.5,'LineStyle',':');
 errorbar(EggLength, AvgProdAllAP(52).AvgProd, AvgProdAllAP(52).All95Conf, 'Color', 'k', 'LineWidth',3.5,'LineStyle',"--");
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 %title('Avg mRNA production per nucleus') 
 ylim([0 1000000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(7).AvgProd./MS2Conversion),(AvgProdAllAP(7).All95Conf./MS2Conversion),'Color',Colors(7).Color,'LineWidth',3.5,'LineStyle','-');
errorbar(EggLength, (AvgProdAllAP(8).AvgProd./MS2Conversion),(AvgProdAllAP(8).All95Conf./MS2Conversion),'Color',Colors(8).Color,'LineWidth',3.5,'LineStyle',':');
errorbar(EggLength, (AvgProdAllAP(52).AvgProd./MS2Conversion),(AvgProdAllAP(52).All95Conf./MS2Conversion),'Color','k','LineWidth',3.5,'LineStyle','--');
ylabel('total transcripts produced');
ylim([0 (1000000/MS2Conversion)]);
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'TotProdKr2xProx_1xBcd','.pdf'],'pdf');

% Proximal w 6xBcd array
figure
 errorbar(EggLength,AvgProdAllAP(2).AvgProd,AvgProdAllAP(2).All95Conf,'Color',Colors(2).Color,'LineWidth',3.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(5).AvgProd, AvgProdAllAP(5).All95Conf,'Color',Colors(5).Color,'LineWidth',3.5,'LineStyle',':');
 errorbar(EggLength, AvgProdAllAP(60).AvgProd, AvgProdAllAP(60).All95Conf, 'Color', 'k', 'LineWidth',3.5,'LineStyle',"--");
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 %title('Avg mRNA production per nucleus') 
 ylim([0 1000000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(2).AvgProd./MS2Conversion),(AvgProdAllAP(2).All95Conf./MS2Conversion),'Color',Colors(2).Color,'LineWidth',3.5,'LineStyle','-');
errorbar(EggLength, (AvgProdAllAP(5).AvgProd./MS2Conversion),(AvgProdAllAP(5).All95Conf./MS2Conversion),'Color',Colors(5).Color,'LineWidth',3.5,'LineStyle',':');
errorbar(EggLength, (AvgProdAllAP(60).AvgProd./MS2Conversion),(AvgProdAllAP(60).All95Conf./MS2Conversion),'Color','k','LineWidth',3.5,'LineStyle','--');
ylabel('total transcripts produced');
ylim([0 (1000000/MS2Conversion)]);
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
%% 6x TFBS arrays compared to hemi vs homozygous
ArrayColor = [179 16 182] ./255;
Array3Color =[155 16 152] ./255;
Array6Color = [115 16 107]./255;
figure
 %errorbar(EggLength,AvgProdAllAP(6).AvgProd,AvgProdAllAP(6).All95Conf,'Color',Colors(6).Color,'LineWidth',3.5);
 lineProps.col{1} = 'k' %Colors(6).Color
 mseb(EggLength, AvgProdAllAP(6).AvgProd, AvgProdAllAP(6).All95Conf, lineProps,0)
 hold on 
 lineProps.col{1} = Colors(26).Color
 %lineProps.style{1} = ':';
 mseb(EggLength, AvgProdAllAP(26).AvgProd, AvgProdAllAP(26).All95Conf, lineProps,0)
 %errorbar(EggLength,AvgProdAllAP(26).AvgProd, AvgProdAllAP(26).All95Conf,'Color',Colors(26).Color,'LineWidth',3.5,'LineStyle',':');
 lineProps.col{1} = ArrayColor;
 %lineProps.style{1} = '--';
 %mseb(EggLength, AvgProdAllAP(39).AvgProd, AvgProdAllAP(39).All95Conf, lineProps,0)
 %errorbar(EggLength, AvgProdAllAP(39).AvgProd, AvgProdAllAP(39).All95Conf, 'Color', 'k', 'LineWidth',3.5,'LineStyle',"--");
 
 mseb(EggLength, AvgProdAllAP(50).AvgProd, AvgProdAllAP(50).All95Conf, lineProps,0)
 %errorbar(EggLength, AvgProdAllAP(50).AvgProd, AvgProdAllAP(50).All95Conf, 'Color', grey, 'LineWidth',3.5,'LineStyle',"--");
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 %title('Avg mRNA production per nucleus') 
 ylim([0 1400000]);
 xlim([20 80]);
 yyaxis right
 lineProps.col{1} = 'k' %Colors(6).Color
 %lineProps.style{1} = '-';
 mseb(EggLength, (AvgProdAllAP(6).AvgProd./MS2Conversion), (AvgProdAllAP(6).All95Conf./MS2Conversion), lineProps,0)
 %errorbar(EggLength, (AvgProdAllAP(6).AvgProd./MS2Conversion),(AvgProdAllAP(6).All95Conf./MS2Conversion),'Color',Colors(6).Color,'LineWidth',3.5,'LineStyle','-');
lineProps.col{1} = Colors(26).Color
%lineProps.style{1} = ':';
 mseb(EggLength, (AvgProdAllAP(26).AvgProd./MS2Conversion), (AvgProdAllAP(26).All95Conf./MS2Conversion), lineProps,0)
 %errorbar(EggLength, (AvgProdAllAP(26).AvgProd./MS2Conversion),(AvgProdAllAP(26).All95Conf./MS2Conversion),'Color',Colors(26).Color,'LineWidth',3.5,'LineStyle',':');
lineProps.col{1} = ArrayColor;
%lineProps.style{1} = '--';
 %mseb(EggLength, (AvgProdAllAP(39).AvgProd./MS2Conversion), (AvgProdAllAP(39).All95Conf./MS2Conversion), lineProps,0)
 %errorbar(EggLength, (AvgProdAllAP(39).AvgProd./MS2Conversion),(AvgProdAllAP(39).All95Conf./MS2Conversion),'Color','k','LineWidth',3.5,'LineStyle','--');
lineProps.col{1} = Array6Color;
 mseb(EggLength, (AvgProdAllAP(50).AvgProd./MS2Conversion), (AvgProdAllAP(50).All95Conf./MS2Conversion), lineProps,0)
 %errorbar(EggLength, (AvgProdAllAP(50).AvgProd./MS2Conversion),(AvgProdAllAP(50).All95Conf./MS2Conversion),'Color',grey,'LineWidth',3.5,'LineStyle','--');
ylabel('total transcripts produced');
ylim([0 (1400000/MS2Conversion)]);
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'TotProdKr2xDist_6xBcd_Shade','.pdf'],'pdf');

% 2xDistal 3xBcd
figure
 %errorbar(EggLength,AvgProdAllAP(6).AvgProd,AvgProdAllAP(6).All95Conf,'Color',Colors(6).Color,'LineWidth',3.5);
 lineProps.col{1} = Colors(6).Color
 mseb(EggLength, AvgProdAllAP(6).AvgProd, AvgProdAllAP(6).All95Conf, lineProps,0)
 hold on 
 lineProps.col{1} = Colors(26).Color
 %lineProps.style{1} = ':';
 mseb(EggLength, AvgProdAllAP(26).AvgProd, AvgProdAllAP(26).All95Conf, lineProps,0)
 %errorbar(EggLength,AvgProdAllAP(26).AvgProd, AvgProdAllAP(26).All95Conf,'Color',Colors(26).Color,'LineWidth',3.5,'LineStyle',':');
 lineProps.col{1} = ArrayColor;
 %lineProps.style{1} = '--';
 %mseb(EggLength, AvgProdAllAP(39).AvgProd, AvgProdAllAP(39).All95Conf, lineProps,0)
 %errorbar(EggLength, AvgProdAllAP(39).AvgProd, AvgProdAllAP(39).All95Conf, 'Color', 'k', 'LineWidth',3.5,'LineStyle',"--");
 
 %mseb(EggLength, AvgProdAllAP(50).AvgProd, AvgProdAllAP(50).All95Conf, lineProps,0)
 mseb(EggLength, AvgProdAllAP(61).AvgProd, AvgProdAllAP(61).All95Conf, lineProps,0)
 %errorbar(EggLength, AvgProdAllAP(50).AvgProd, AvgProdAllAP(50).All95Conf, 'Color', grey, 'LineWidth',3.5,'LineStyle',"--");
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 %title('Avg mRNA production per nucleus') 
 ylim([0 1400000]);
 xlim([20 80]);
 yyaxis right
 lineProps.col{1} = Colors(6).Color
 %lineProps.style{1} = '-';
 mseb(EggLength, (AvgProdAllAP(6).AvgProd./MS2Conversion), (AvgProdAllAP(6).All95Conf./MS2Conversion), lineProps,0)
 %errorbar(EggLength, (AvgProdAllAP(6).AvgProd./MS2Conversion),(AvgProdAllAP(6).All95Conf./MS2Conversion),'Color',Colors(6).Color,'LineWidth',3.5,'LineStyle','-');
lineProps.col{1} = Colors(26).Color
%lineProps.style{1} = ':';
 mseb(EggLength, (AvgProdAllAP(26).AvgProd./MS2Conversion), (AvgProdAllAP(26).All95Conf./MS2Conversion), lineProps,0)
 %errorbar(EggLength, (AvgProdAllAP(26).AvgProd./MS2Conversion),(AvgProdAllAP(26).All95Conf./MS2Conversion),'Color',Colors(26).Color,'LineWidth',3.5,'LineStyle',':');
lineProps.col{1} = ArrayColor;
%lineProps.style{1} = '--';
 %mseb(EggLength, (AvgProdAllAP(39).AvgProd./MS2Conversion), (AvgProdAllAP(39).All95Conf./MS2Conversion), lineProps,0)
 %errorbar(EggLength, (AvgProdAllAP(39).AvgProd./MS2Conversion),(AvgProdAllAP(39).All95Conf./MS2Conversion),'Color','k','LineWidth',3.5,'LineStyle','--');
lineProps.col{1} = Array6Color;
 %mseb(EggLength, (AvgProdAllAP(50).AvgProd./MS2Conversion), (AvgProdAllAP(50).All95Conf./MS2Conversion), lineProps,0)
lineProps.col{1} = Array3Color;
 mseb(EggLength, (AvgProdAllAP(61).AvgProd./MS2Conversion), (AvgProdAllAP(61).All95Conf./MS2Conversion), lineProps,0)
 %errorbar(EggLength, (AvgProdAllAP(50).AvgProd./MS2Conversion),(AvgProdAllAP(50).All95Conf./MS2Conversion),'Color',grey,'LineWidth',3.5,'LineStyle','--');
ylabel('total transcripts produced');
ylim([0 (1400000/MS2Conversion)]);
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'TotProdKr2xDist_3xBcd_Shade','.pdf'],'pdf');

% KrSE w 6xBcd
ArrayColor = [179 16 182] ./255;
Array6Color = [115 16 107]./255;
figure

lineProps.col{1} = 'k'%[Colors(9).Color];
mseb(EggLength, AvgProdAllAP(9).AvgProd, AvgProdAllAP(9).All95Conf, lineProps,0)
  hold on 
  lineProps.col{1} = [Colors(9).Color];
mseb(EggLength, AvgProdAllAP(10).AvgProd, AvgProdAllAP(10).All95Conf, lineProps,0);
lineProps.col{1} = ArrayColor;
mseb(EggLength, AvgProdAllAP(41).AvgProd, AvgProdAllAP(41).All95Conf, lineProps,0);
lineProps.col{1}=Array6Color;
mseb(EggLength, AvgProdAllAP(56).AvgProd, AvgProdAllAP(56).All95Conf, lineProps,0);
set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 %title('Avg mRNA production per nucleus') 
 ylim([0 1400000]);
 xlim([0 100]);
 yyaxis right
 lineProps.col{1} = 'k'%Colors(9).Color
mseb(EggLength, (AvgProdAllAP(9).AvgProd./MS2Conversion), (AvgProdAllAP(9).All95Conf./MS2Conversion), lineProps,0)
  hold on 
  lineProps.col{1} = [Colors(9).Color];
mseb(EggLength, (AvgProdAllAP(10).AvgProd./MS2Conversion), (AvgProdAllAP(10).All95Conf./MS2Conversion), lineProps,0);
lineProps.col{1} = ArrayColor;
mseb(EggLength, (AvgProdAllAP(41).AvgProd./MS2Conversion), (AvgProdAllAP(41).All95Conf./MS2Conversion), lineProps,0);
lineProps.col{1}=Array6Color;
mseb(EggLength, (AvgProdAllAP(56).AvgProd./MS2Conversion), (AvgProdAllAP(56).All95Conf./MS2Conversion), lineProps,0);

ylabel('total transcripts produced');
ylim([0 (1400000/MS2Conversion)]);
xlim([20 80]);
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'TotProdKrSE_6xBcd_Shade','.pdf'],'pdf');

% KrSE w 3xBcd
ArrayColor = [179 16 182] ./255;
Array3Color =[155 16 152] ./255;
Array6Color = [115 16 107]./255;
figure

lineProps.col{1} = 'k'%[Colors(9).Color];
mseb(EggLength, AvgProdAllAP(9).AvgProd, AvgProdAllAP(9).All95Conf, lineProps,0)
  hold on 
  lineProps.col{1} = Colors(9).Color
mseb(EggLength, AvgProdAllAP(10).AvgProd, AvgProdAllAP(10).All95Conf, lineProps,0);
lineProps.col{1} = ArrayColor;
%mseb(EggLength, AvgProdAllAP(41).AvgProd, AvgProdAllAP(41).All95Conf, lineProps,0);
lineProps.col{1} = Array3Color;
mseb(EggLength, AvgProdAllAP(64).AvgProd, AvgProdAllAP(64).All95Conf, lineProps,0);
lineProps.col{1}=Array6Color;
%mseb(EggLength, AvgProdAllAP(56).AvgProd, AvgProdAllAP(56).All95Conf, lineProps,0);
set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 %title('Avg mRNA production per nucleus') 
 ylim([0 1500000]);
 xlim([0 100]);
 yyaxis right
 lineProps.col{1} = 'k'%Colors(9).Color
mseb(EggLength, (AvgProdAllAP(9).AvgProd./MS2Conversion), (AvgProdAllAP(9).All95Conf./MS2Conversion), lineProps,0)
  hold on 
  lineProps.col{1} = Colors(9).Color
mseb(EggLength, (AvgProdAllAP(10).AvgProd./MS2Conversion), (AvgProdAllAP(10).All95Conf./MS2Conversion), lineProps,0);
lineProps.col{1} = ArrayColor;
%mseb(EggLength, (AvgProdAllAP(41).AvgProd./MS2Conversion), (AvgProdAllAP(41).All95Conf./MS2Conversion), lineProps,0);
lineProps.col{1} = Array3Color;
mseb(EggLength, (AvgProdAllAP(64).AvgProd./MS2Conversion), (AvgProdAllAP(64).All95Conf./MS2Conversion), lineProps,0);
lineProps.col{1}=Array6Color;
%mseb(EggLength, (AvgProdAllAP(56).AvgProd./MS2Conversion), (AvgProdAllAP(56).All95Conf./MS2Conversion), lineProps,0);

ylabel('total transcripts produced');
ylim([0 (1500000/MS2Conversion)]);
xlim([20 80]);
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'TotProdKrSE_3xBcd_Shade','.pdf'],'pdf');

%% Homozygotes as fx of hemizygotes plots
%Distal 
figure
plot(AvgProdAllAP(4).AvgProd, AvgProdAllAP(1).AvgProd,'o', 'Color',Colors(1).Color, 'LineWidth',2.5);
hold on
hline = refline(1,0);
hline.LineWidth =2.5;
hline.Color ='k'
ylabel('homozygous average total mRNA per allele');
xlabel('hemizygous average total mRNA per allele');
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);

%Proximal
figure
plot(AvgProdAllAP(5).AvgProd, AvgProdAllAP(2).AvgProd,'o', 'Color',Colors(2).Color, 'LineWidth',2.5);
hold on
hline = refline(1,0);
hline.LineWidth =2.5;
hline.Color ='k'
ylabel('homozygous average total mRNA per allele');
xlabel('hemizygous average total mRNA per allele');
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);

%SE pair
figure
plot(AvgProdAllAP(10).AvgProd, AvgProdAllAP(9).AvgProd,'o', 'Color',Colors(9).Color, 'LineWidth',2.5);
hold on
hline = refline(1,0);
hline.LineWidth =2.5;
hline.Color ='k'
ylabel('homozygous average total mRNA per allele');
xlabel('hemizygous average total mRNA per allele');
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);

% 2xDistal
figure
plot(AvgProdAllAP(26).AvgProd, AvgProdAllAP(6).AvgProd,'o', 'Color',Colors(6).Color, 'LineWidth',2.5);
hold on
hline = refline(1,0);
hline.LineWidth =2.5;
hline.Color ='k'
ylabel('homozygous average total mRNA per allele');
xlabel('hemizygous average total mRNA per allele');
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);

%2xProx 
figure
plot(AvgProdAllAP(8).AvgProd, AvgProdAllAP(7).AvgProd,'o', 'Color',Colors(7).Color, 'LineWidth',2.5);
hold on
hline = refline(1,0);
hline.LineWidth =2.5;
hline.Color ='k'
ylabel('homozygous average total mRNA per allele');
xlabel('hemizygous average total mRNA per allele');
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);

% A single combined graph 
figure 
Constructs_Homo =[1,2,9,6,7];
Constructs_Hemi = [4,5,10,26,8];
APBins_Homo = [ExpressionBins{1,[1,2,9,6,7]}];
APBins_Hemi = [ExpressionBins{1,[4,5,10,26,8]}];
for ss = 1:length(APBins_Homo)
    XError = AvgProdAllAP(Constructs_Hemi(ss)).All95Conf(APBins_Hemi(ss));
    YError = AvgProdAllAP(Constructs_Homo(ss)).All95Conf(APBins_Homo(ss));
    
    errorbar(AvgProdAllAP(Constructs_Hemi(ss)).AvgProd(APBins_Hemi(ss)), AvgProdAllAP(Constructs_Homo(ss)).AvgProd(APBins_Homo(ss)),YError, YError, XError, XError,'o','Color',Colors(Constructs_Homo(ss)).Color,'LineWidth',2.5);
    %plot(AvgProdAllAP(Constructs_Hemi(ss)).AvgProd(APBins_Hemi(ss)), AvgProdAllAP(Constructs_Homo(ss)).AvgProd(APBins_Homo(ss)),'o','LineWidth',2.5,'Color',Colors(Constructs_Homo(ss)).Color);
    hold on 
end
XVals = [4e5, 6e5, 10e5, 14e5];
YVals = XVals;
plot(XVals, YVals,'Color','k','LineWidth',2.5,'LineStyle','--')
area(XVals, YVals,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.2);
%hline = refline(1,0);
%hline.Color = 'k'; hline.LineWidth = 2;
ylim([4e5 14e5]);
xlim([4e5 14e5]);
ylabel('homozygous total mRNA per allele');
xlabel('hemizygous total mRNA per allele');
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'HomoExp_Fx_HemiExp','.pdf'],'pdf');

%% Get % of hemizygous expression for all constructs 
DistHemiAP = ExpressionBins{1,4};
ProxHemiAP = ExpressionBins{1,5};
DDHemiAP = ExpressionBins{1,26};
DoubProxHemi = ExpressionBins{1,8};
SEHemiAP = ExpressionBins{1,10};
HbHemiAP = ExpressionBins{1,16};

DistComps = [1,63];
ProxComps = [2,60]; % homo and 6xBcd 
DDComps = [6,39,42,45,50,61];
DPComps = [7,43,46,52];
SEComps = [9,41,44,47,48,56,64];
HbP2Comps = [57,66,59];

% Initiate table to fill in 
HemiChangeTable = table();

% Distal enhancer
Counter=0;
for cc = DistComps
    Counter = Counter+1;
    DistChange(1,Counter) = (AvgProdAllAP(cc).AvgProd(ExpressionBins{1,cc}) - AvgProdAllAP(4).AvgProd(ExpressionBins{1,4}))/AvgProdAllAP(4).AvgProd(DistHemiAP);
    % Perc change assuming hemi is "final"
    DistChange(2,Counter) = (AvgProdAllAP(4).AvgProd(DistHemiAP))/AvgProdAllAP(cc).AvgProd(ExpressionBins{1,cc});

    % Save in an excel sheet for easy generation/reference
    HemiChangeTable.(ConstructList{cc}) = DistChange(:,Counter);
end

% Proximal enhancer
Counter = 0;
for cc = ProxComps
    Counter = Counter+1;
    ProxChange(1,Counter) = (AvgProdAllAP(cc).AvgProd(ExpressionBins{1,cc}) - AvgProdAllAP(5).AvgProd(ExpressionBins{1,5}))/AvgProdAllAP(5).AvgProd(ProxHemiAP);
    % Perc change assuming hemi is "final"
    ProxChange(2,Counter) = (AvgProdAllAP(5).AvgProd(ProxHemiAP))/AvgProdAllAP(cc).AvgProd(ExpressionBins{1,cc});

    % Save in an excel sheet for easy generation/reference
    HemiChangeTable.(ConstructList{cc}) = ProxChange(:,Counter);
end

% 2xDistal 
Counter = 0;
for cc = DDComps
    Counter = Counter+1;
    DDChange(1,Counter) = (AvgProdAllAP(cc).AvgProd(ExpressionBins{1,cc}) - AvgProdAllAP(26).AvgProd(ExpressionBins{1,26}))/AvgProdAllAP(26).AvgProd(DDHemiAP);
    % Perc change assuming hemi is "final"
    DDChange(2,Counter) = (AvgProdAllAP(26).AvgProd(DDHemiAP))/AvgProdAllAP(cc).AvgProd(ExpressionBins{1,cc});

    % Save in an excel sheet for easy generation/reference
    HemiChangeTable.(ConstructList{cc}) = DDChange(:,Counter);
end

% 2xProximal
Counter = 0;
for cc = DPComps
    Counter = Counter+1;
    DPChange(1,Counter) = (AvgProdAllAP(cc).AvgProd(ExpressionBins{1,cc})- AvgProdAllAP(8).AvgProd(ExpressionBins{1,8}))/AvgProdAllAP(8).AvgProd(DoubProxHemi);
    % Perc change assuming hemi is "final"
    DPChange(2,Counter) = (AvgProdAllAP(8).AvgProd(DoubProxHemi))/AvgProdAllAP(cc).AvgProd(ExpressionBins{1,cc});

    % Save in an excel sheet for easy generation/reference
    HemiChangeTable.(ConstructList{cc}) = DPChange(:,Counter);
end

% SE
Counter = 0;
for cc = SEComps
    Counter = Counter+1;
    SEChange(1,Counter) = (AvgProdAllAP(cc).AvgProd(ExpressionBins{1,cc}) - AvgProdAllAP(10).AvgProd(ExpressionBins{1,10}))/AvgProdAllAP(10).AvgProd(SEHemiAP);
    % Perc change assuming hemi is "final"
    SEChange(2,Counter) = (AvgProdAllAP(10).AvgProd(SEHemiAP))/AvgProdAllAP(cc).AvgProd(ExpressionBins{1,cc});

    % Save in an excel sheet for easy generation/reference
    HemiChangeTable.(ConstructList{cc}) = SEChange(:,Counter);
end

% HbP2
Counter = 0;
for cc = HbP2Comps
    Counter = Counter+1;
    HbChange(1,Counter) = (AvgProdAllAP(cc).AvgProd(ExpressionBins{1,cc}) - AvgProdAllAP(16).AvgProd(ExpressionBins{1,16}))/AvgProdAllAP(16).AvgProd(HbHemiAP);
    % Perc change assuming hemi is "final"
    HbChange(2,Counter) = (AvgProdAllAP(16).AvgProd(HbHemiAP))/AvgProdAllAP(cc).AvgProd(ExpressionBins{1,cc});

    % Save in an excel sheet for easy generation/reference
    HemiChangeTable.(ConstructList{cc}) = HbChange(:,Counter);
end
    
writetable(HemiChangeTable, ['C:/Users/Rachel2/Documents/livemRNA\Data\DynamicsResults' filesep 'Constructs' filesep 'ExpChangefromHemi.csv'] );

%% Test if variance is different in hemi vs homozygotes
    %Distal
    x = [AvgProdAllAP(1).AllProds(:,ExpressionBins{1,1})];
    y = [AvgProdAllAP(4).AllProds(:,ExpressionBins{1,4})];
[h,p,ci,stats] = vartest2(x,y);

% Proximal
x = [AvgProdAllAP(2).AllProds(:,ExpressionBins{1,2})];
y = [AvgProdAllAP(5).AllProds(:,ExpressionBins{1,5})];
[h,p,ci,stats] = vartest2(x,y);

% Shadow pair
x = [AvgProdAllAP(9).AllProds(:,ExpressionBins{1,9})];
y = [AvgProdAllAP(10).AllProds(:,ExpressionBins{1,10})];
[h,p,ci,stats] = vartest2(x,y);

% 2x Distal
x = [AvgProdAllAP(6).AllProds(:,ExpressionBins{1,6})];
y = [AvgProdAllAP(26).AllProds(:,ExpressionBins{1,26})];
[h,p,ci,stats] = vartest2(x,y);

% 2x Proximal
x = [AvgProdAllAP(7).AllProds(:,ExpressionBins{1,7})];
y = [AvgProdAllAP(8).AllProds(:,ExpressionBins{1,8})];
[h,p,ci,stats] = vartest2(x,y);
%% Arrays as fx of hemi production
% Graph of 2xDistal with TFBS arrays
figure
DupDist_TFs = [39,45,42,58];
APExpBins_DDTFs = [ExpressionBins{1,[DupDist_TFs]}];
AP_DDHemi = ExpressionBins{1,26};
    % plot dashed lines for hemi/homo expression levels
h = line([0,10],[AvgProdAllAP(26).AvgProd(AP_DDHemi),AvgProdAllAP(26).AvgProd(AP_DDHemi)]);
h.Color = Colors(26).Color; h.LineStyle ='--'; h.LineWidth =2.5;
hold on 
h2 = line([0,10],[AvgProdAllAP(6).AvgProd(ExpressionBins{1,6}),AvgProdAllAP(6).AvgProd(ExpressionBins{1,6})]);
h2.Color = Colors(6).Color; h2.LineStyle ='--'; h2.LineWidth =2.5;
ss=1;
for tt = 1:length(DupDist_TFs)
    errorbar(ss, AvgProdAllAP(DupDist_TFs(tt)).AvgProd(APExpBins_DDTFs(tt)),AvgProdAllAP(DupDist_TFs(tt)).All95Conf(APExpBins_DDTFs(tt)),'o','Color','k','LineWidth',2.5);
    ss = ss+2;
end
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
xticks([1 3 5 7 9]);
xticklabels({'1xBcd', '1xZld', '1xHb', '1xTwi'});
ylabel('peak total mRNA produced');
ylim([0.5e6, 1.5e6]);
xlim([0 8]);
saveas(gcf, [FigDirect filesep 'DupDist_TFBSArrays_HemiHomoLines','.pdf'],'pdf');

% Graph of shadow pair with TFBS arrays
figure
SE_TFs = [41,47,44,48];
APExpBins_SETFs = [ExpressionBins{1,[SE_TFs]}];
AP_SEHemi = ExpressionBins{1,10};
    % plot dashed lines for hemi/homo expression levels -- w shading for
    % error
lineProps.col{1} = [Colors(9).Color];
mseb([-9:10], [ones(1,20).*AvgProdAllAP(10).AvgProd(ExpressionBins{1,10})], ones(1,20).*AvgProdAllAP(10).All95Conf(ExpressionBins{1,10}), lineProps,0)
hold on 
mseb([-9:10], [ones(1,20).*AvgProdAllAP(9).AvgProd(ExpressionBins{1,9})], ones(1,20).*AvgProdAllAP(9).All95Conf(ExpressionBins{1,9}), lineProps,0)
%h2 = line([0,10],[AvgProdAllAP(9).AvgProd(ExpressionBins{1,9}),AvgProdAllAP(9).AvgProd(ExpressionBins{1,9})]);
%h2.Color = Colors(9).Color; h2.LineStyle ='--'; h2.LineWidth =2.5;
ss=1;
for tt = 1:length(SE_TFs)
    errorbar(ss, AvgProdAllAP(SE_TFs(tt)).AvgProd(APExpBins_SETFs(tt)),AvgProdAllAP(SE_TFs(tt)).All95Conf(APExpBins_SETFs(tt)),'o','Color','k','LineWidth',2.5);
    ss = ss+2;
end
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 5]);
xticks([1 3 5 7 9]);
xticklabels({'1xBcd', '1xZld', '1xHb', '1xStat92E'});
ylabel('peak total mRNA produced');
ylim([5e5, 1.5e6]);
xlim([0 8]);
saveas(gcf, [FigDirect filesep 'SE_TFBSArrays_HemiHomoLines_Shade','.pdf'],'pdf');

% Single graphs for 3 main constructs w 1xBcd
  %shadow pair -- w shading around hemi/homo lines for error
figure
lineProps.col{1} = [Colors(9).Color];
mseb([-9:10], [ones(1,20).*AvgProdAllAP(10).AvgProd(AP_SEHemi)], ones(1,20).*AvgProdAllAP(10).All95Conf(AP_SEHemi), lineProps,0)
%h = line([0,10],[AvgProdAllAP(10).AvgProd(AP_SEHemi),AvgProdAllAP(10).AvgProd(AP_SEHemi)]);
%h.Color = Colors(10).Color; h.LineStyle ='--'; h.LineWidth =2.5;
hold on 
mseb([-9:10], ones(1,20).*AvgProdAllAP(9).AvgProd(ExpressionBins{1,9}), ones(1,20).*AvgProdAllAP(9).All95Conf(ExpressionBins{1,9}), lineProps,0)
%h2 = line([0,10],[AvgProdAllAP(9).AvgProd(ExpressionBins{1,9}),AvgProdAllAP(9).AvgProd(ExpressionBins{1,9})]);
%h2.Color = Colors(9).Color; h2.LineStyle ='--'; h2.LineWidth =2.5;
errorbar(1, AvgProdAllAP(41).AvgProd(ExpressionBins{1,41}),AvgProdAllAP(41).All95Conf(ExpressionBins{1,41}),'o','Color','k','LineWidth',2.5);
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
ylabel('integrated fluorescence (AU)');
xticks([1]);
xticklabels('1x Bcd');
xlim([0 2]);
ylim([3e5, 1.5e6]);
saveas(gcf, [FigDirect filesep 'SE_1xBcd_HemiHomoLines_Shade','.pdf'],'pdf');

 %duplicated distal -- shaded error
x_sizeNarrow =1.3;
y_sizeNarrow =3.97;
figure
lineProps.col{1} = [Colors(6).Color];
mseb([-9:10], [ones(1,20).*AvgProdAllAP(26).AvgProd(AP_DDHemi)], ones(1,20).*AvgProdAllAP(26).All95Conf(AP_DDHemi), lineProps,0)
hold on
mseb([-9:10], [ones(1,20).*AvgProdAllAP(6).AvgProd(ExpressionBins{1,6})], ones(1,20).*AvgProdAllAP(6).All95Conf(ExpressionBins{1,6}), lineProps,0)
%h2 = line([0,10],[AvgProdAllAP(6).AvgProd(ExpressionBins{1,6}),AvgProdAllAP(6).AvgProd(ExpressionBins{1,6})]);
%h2.Color = Colors(6).Color; h2.LineStyle ='--'; h2.LineWidth =2.5;
errorbar(0.2, AvgProdAllAP(39).AvgProd(ExpressionBins{1,39}),AvgProdAllAP(39).All95Conf(ExpressionBins{1,39}),'o','Color','k','LineWidth',2.5);
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 2 5]);
ylabel('integrated fluorescence (AU)');
xticks([0.2]);
xticklabels('1x Bcd');
xlim([0 0.4]);
ylim([3e5, 1.5e6]);
saveas(gcf, [FigDirect filesep 'DupDist_1xBcd_HemiHomoLines_Shade','.pdf'],'pdf');

 %duplicated proximal
figure
lineProps.col{1} = [Colors(7).Color];
mseb([-9:10], [ones(1,20).*AvgProdAllAP(8).AvgProd(ExpressionBins{1,8})], ones(1,20).*AvgProdAllAP(8).All95Conf(ExpressionBins{1,8}), lineProps,0)
%h = line([0,10],[AvgProdAllAP(8).AvgProd(ExpressionBins{1,8}),AvgProdAllAP(8).AvgProd(ExpressionBins{1,8})]);
%h.Color = Colors(7).Color; h.LineStyle ='--'; h.LineWidth =2.5;
hold on 
mseb([-9:10], [ones(1,20).*AvgProdAllAP(7).AvgProd(ExpressionBins{1,7})], ones(1,20).*AvgProdAllAP(7).All95Conf(ExpressionBins{1,7}), lineProps,0)
%h2 = line([0,10],[AvgProdAllAP(7).AvgProd(ExpressionBins{1,7}),AvgProdAllAP(7).AvgProd(ExpressionBins{1,7})]);
%h2.Color = Colors(7).Color; h2.LineStyle ='--'; h2.LineWidth =2.5;
errorbar(0.2, AvgProdAllAP(52).AvgProd(ExpressionBins{1,52}),AvgProdAllAP(52).All95Conf(ExpressionBins{1,52}),'o','Color','k','LineWidth',2.5);
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 2 5]);
ylabel('integrated fluorescence (AU)');
xticks([0.2]);
xticklabels('1x Bcd');
xlim([0 0.4]);
ylim([3e5, 1.5e6]);
saveas(gcf, [FigDirect filesep 'DupProx_1xBcd_HemiHomoLines_Shade','.pdf'],'pdf');

 %HbP2 w 6xBcd -- Figure 4
x_sizeNarrow =1.3;
y_sizeNarrow =3.97;
figure
lineProps.col{1} = [0,0,0];
mseb([-9:10], [ones(1,20).*AvgProdAllAP(16).AvgProd(ExpressionBins{1,16})], ones(1,20).*AvgProdAllAP(16).All95Conf(ExpressionBins{1,16}), lineProps,0)
hold on
mseb([-9:10], [ones(1,20).*AvgProdAllAP(57).AvgProd(ExpressionBins{1,57})], ones(1,20).*AvgProdAllAP(57).All95Conf(ExpressionBins{1,57}), lineProps,0)
%h2 = line([0,10],[AvgProdAllAP(6).AvgProd(ExpressionBins{1,6}),AvgProdAllAP(6).AvgProd(ExpressionBins{1,6})]);
%h2.Color = Colors(6).Color; h2.LineStyle ='--'; h2.LineWidth =2.5;
errorbar(0.2, AvgProdAllAP(59).AvgProd(ExpressionBins{1,59}),AvgProdAllAP(59).All95Conf(ExpressionBins{1,59}),'o','Color','k','LineWidth',2.5);
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 2 5]);
ylabel('integrated fluorescence (AU)');
xticks([0.2]);
xticklabels('6x Bcd');
xlim([0 0.4]);
ylim([3e5 1.6e6]);
saveas(gcf, [FigDirect filesep 'HbP2_6xBcd_HemiHomoLines_Shade','.pdf'],'pdf');

% HbP2 w 1xBcd array - supplement (?)
x_sizeNarrow =1.3;
y_sizeNarrow =3.97;
figure
lineProps.col{1} = [0,0,0];
mseb([-9:10], [ones(1,20).*AvgProdAllAP(16).AvgProd(ExpressionBins{1,16})], ones(1,20).*AvgProdAllAP(16).All95Conf(ExpressionBins{1,16}), lineProps,0)
hold on
mseb([-9:10], [ones(1,20).*AvgProdAllAP(57).AvgProd(ExpressionBins{1,57})], ones(1,20).*AvgProdAllAP(57).All95Conf(ExpressionBins{1,57}), lineProps,0)
%h2 = line([0,10],[AvgProdAllAP(6).AvgProd(ExpressionBins{1,6}),AvgProdAllAP(6).AvgProd(ExpressionBins{1,6})]);
%h2.Color = Colors(6).Color; h2.LineStyle ='--'; h2.LineWidth =2.5;
errorbar(0.2, AvgProdAllAP(66).AvgProd(ExpressionBins{1,66}),AvgProdAllAP(66).All95Conf(ExpressionBins{1,66}),'o','Color','k','LineWidth',2.5);
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 2 5]);
ylabel('integrated fluorescence (AU)');
xticks([0.2]);
xticklabels('1x Bcd');
xlim([0 0.4]);
ylim([3e5 1.6e6]);
saveas(gcf, [FigDirect filesep 'HbP2_1xBcd_HemiHomoLines_Shade','.pdf'],'pdf');



% Dose-dependent graphs for Bcd binding site arrays 
    %shadow pair
DupDist_Bcd = [39,61,50];
SE_Bcd = [41,64,56];
Hb_Bcd = [66, 59];
APExpBins_DupDist = [ExpressionBins{1,[DupDist_Bcd]}];
APExpBins_SE = [ExpressionBins{1,[SE_Bcd]}];
APExpBins_Hb = [ExpressionBins{1,[Hb_Bcd]}];
AP_DDHemi = ExpressionBins{1,26};
AP_SEHemi = ExpressionBins{1,10};
AP_HbHemi = ExpressionBins{1,16};
%DupDist_Colors = [Colors(6).Color;[73 99 253]./253; ArrayColor; Array3Color; Array6Color];

%shadow pair data -- shade error bars
 figure
lineProps.col{1} = [Colors(9).Color];
mseb([-9:10], [ones(1,20).*AvgProdAllAP(10).AvgProd(ExpressionBins{1,10})], ones(1,20).*AvgProdAllAP(10).All95Conf(ExpressionBins{1,10}), lineProps,0)
hold on 
mseb([-9:10], [ones(1,20).*AvgProdAllAP(9).AvgProd(ExpressionBins{1,9})], ones(1,20).*AvgProdAllAP(9).All95Conf(ExpressionBins{1,9}), lineProps,0)
% h = line([0,10],[AvgProdAllAP(9).AvgProd(ExpressionBins{1,9}),AvgProdAllAP(9).AvgProd(ExpressionBins{1,9})]);
% h.Color = Colors(9).Color; h.LineStyle ='--'; h.LineWidth =2.5;
pp=1;
 for ss = 1:length(SE_Bcd)
     errorbar(pp,AvgProdAllAP(SE_Bcd(ss)).AvgProd(APExpBins_SE(ss)), AvgProdAllAP(SE_Bcd(ss)).All95Conf(APExpBins_SE(ss)),'o','Color','k','LineWidth',2.5);
     pp = pp+2;
     hold on
 end
%  XVals = [0.4e6,1e6,1.1e6,1.2e6,1.3e6,1.4e6,1.6e6];
% YVals = XVals;
% plot(XVals, YVals,'Color','k','LineWidth',2.5,'LineStyle','--')
xlim([0 6]);
ylim([3e5 1.6e6]);
ylabel('integrated fluorescence (AU)');
xlabel('# of Bcd binding sites')
xticks([1,3,5]);
xticklabels({'6','18', '36'});
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'SE_BcdArrays_HemiHomoLines','.pdf'],'pdf');

 %duplicated distal data -- shade errorbars
 figure
lineProps.col{1} = [Colors(6).Color];
mseb([-9:10], [ones(1,20).*AvgProdAllAP(26).AvgProd(ExpressionBins{1,26})], ones(1,20).*AvgProdAllAP(26).All95Conf(ExpressionBins{1,26}), lineProps,0)
hold on 
mseb([-9:10], [ones(1,20).*AvgProdAllAP(6).AvgProd(ExpressionBins{1,6})], ones(1,20).*AvgProdAllAP(6).All95Conf(ExpressionBins{1,6}), lineProps,0)
% h = line([0,10],[AvgProdAllAP(6).AvgProd(ExpressionBins{1,6}),AvgProdAllAP(6).AvgProd(ExpressionBins{1,6})]);
% h.Color = Colors(6).Color; h.LineStyle ='--'; h.LineWidth =2.5;
pp=1;
 for ss = 1:length(DupDist_Bcd)
     errorbar(pp,AvgProdAllAP(DupDist_Bcd(ss)).AvgProd(APExpBins_DupDist(ss)), AvgProdAllAP(DupDist_Bcd(ss)).All95Conf(APExpBins_DupDist(ss)),'o','Color','k','LineWidth',2.5);
     pp = pp+2;
     hold on
 end
%  XVals = [0.4e6,1e6,1.1e6,1.2e6,1.3e6,1.4e6,1.6e6];
% YVals = XVals;
% plot(XVals, YVals,'Color','k','LineWidth',2.5,'LineStyle','--')
xlim([0 6]);
ylim([3e5 1.6e6]);
ylabel('integrated fluorescence (AU)');
xlabel('# of Bcd binding arrays')
xticks([1,3,5]);
xticklabels({'1x','3x', '6x'});
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'DupDist_BcdArrays_HemiHomoLines_Shade','.pdf'],'pdf');

%HbP2 data -- shade error bars
 figure
lineProps.col{1} = ['k'];
mseb([-9:10], [ones(1,20).*AvgProdAllAP(16).AvgProd(ExpressionBins{1,16})], ones(1,20).*AvgProdAllAP(16).All95Conf(ExpressionBins{1,16}), lineProps,0)
hold on 
mseb([-9:10], [ones(1,20).*AvgProdAllAP(57).AvgProd(ExpressionBins{1,57})], ones(1,20).*AvgProdAllAP(57).All95Conf(ExpressionBins{1,57}), lineProps,0)
% h = line([0,10],[AvgProdAllAP(9).AvgProd(ExpressionBins{1,9}),AvgProdAllAP(9).AvgProd(ExpressionBins{1,9})]);
% h.Color = Colors(9).Color; h.LineStyle ='--'; h.LineWidth =2.5;
pp=1;
 for ss = 1:length(Hb_Bcd)
     errorbar(pp,AvgProdAllAP(Hb_Bcd(ss)).AvgProd(APExpBins_Hb(ss)), AvgProdAllAP(Hb_Bcd(ss)).All95Conf(APExpBins_Hb(ss)),'o','Color','k','LineWidth',2.5);
     pp = pp+2;
     hold on
 end
%  XVals = [0.4e6,1e6,1.1e6,1.2e6,1.3e6,1.4e6,1.6e6];
% YVals = XVals;
% plot(XVals, YVals,'Color','k','LineWidth',2.5,'LineStyle','--')
xlim([0 4]);
ylim([1e5 1.2e6]);
ylabel('integrated fluorescence (AU)');
xlabel('# of Bcd binding arrays')
xticks([1,3]);
xticklabels({'1x','6x'});
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 4 6]);
saveas(gcf, [FigDirect filesep 'HbP2_BcdArrays_HemiHomoLines','.pdf'],'pdf');

%% Supplemental TFBS array effect plots 

% 2xDistal w Zld
x_sizeNarrow =1.3;
y_sizeNarrow =3.97;
figure
lineProps.col{1} = Colors(6).Color;
mseb([-9:10], [ones(1,20).*AvgProdAllAP(6).AvgProd(ExpressionBins{1,6})], ones(1,20).*AvgProdAllAP(6).All95Conf(ExpressionBins{1,6}), lineProps,0)
hold on
mseb([-9:10], [ones(1,20).*AvgProdAllAP(26).AvgProd(ExpressionBins{1,26})], ones(1,20).*AvgProdAllAP(26).All95Conf(ExpressionBins{1,26}), lineProps,0)

errorbar(0.2, AvgProdAllAP(45).AvgProd(ExpressionBins{1,45}),AvgProdAllAP(45).All95Conf(ExpressionBins{1,45}),'o','Color','k','LineWidth',2.5);
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 2 5]);
ylabel('integrated fluorescence (AU)');
xticks([0.2]);
xticklabels('1x Zld');
xlim([0 0.4]);
ylim([3e5 1.6e6]);


% 2xDistal w Hb
figure
lineProps.col{1} = Colors(6).Color;
mseb([-9:10], [ones(1,20).*AvgProdAllAP(6).AvgProd(ExpressionBins{1,6})], ones(1,20).*AvgProdAllAP(6).All95Conf(ExpressionBins{1,6}), lineProps,0)
hold on
mseb([-9:10], [ones(1,20).*AvgProdAllAP(26).AvgProd(ExpressionBins{1,26})], ones(1,20).*AvgProdAllAP(26).All95Conf(ExpressionBins{1,26}), lineProps,0)

errorbar(0.2, AvgProdAllAP(42).AvgProd(ExpressionBins{1,42}),AvgProdAllAP(42).All95Conf(ExpressionBins{1,42}),'o','Color','k','LineWidth',2.5);
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 2 5]);
ylabel('integrated fluorescence (AU)');
xticks([0.2]);
xticklabels('1x Hb');
xlim([0 0.4]);
ylim([3e5 1.6e6]);

% 2xdistal w Stat92E (to be imaged)



% 2xProx w Stat92E
figure
lineProps.col{1} = Colors(7).Color;
mseb([-9:10], [ones(1,20).*AvgProdAllAP(7).AvgProd(ExpressionBins{1,7})], ones(1,20).*AvgProdAllAP(7).All95Conf(ExpressionBins{1,7}), lineProps,0)
hold on
mseb([-9:10], [ones(1,20).*AvgProdAllAP(8).AvgProd(ExpressionBins{1,8})], ones(1,20).*AvgProdAllAP(8).All95Conf(ExpressionBins{1,8}), lineProps,0)

errorbar(0.2, AvgProdAllAP(46).AvgProd(ExpressionBins{1,46}),AvgProdAllAP(46).All95Conf(ExpressionBins{1,46}),'o','Color','k','LineWidth',2.5);
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 2 5]);
ylabel('integrated fluorescence (AU)');
xticks([0.2]);
xticklabels('1x Stat92E');
xlim([0 0.4]);
ylim([3e5 1.6e6]);


% 2xProx w Hb
figure
lineProps.col{1} = Colors(7).Color;
mseb([-9:10], [ones(1,20).*AvgProdAllAP(7).AvgProd(ExpressionBins{1,7})], ones(1,20).*AvgProdAllAP(7).All95Conf(ExpressionBins{1,7}), lineProps,0)
hold on
mseb([-9:10], [ones(1,20).*AvgProdAllAP(8).AvgProd(ExpressionBins{1,8})], ones(1,20).*AvgProdAllAP(8).All95Conf(ExpressionBins{1,8}), lineProps,0)

errorbar(0.2, AvgProdAllAP(43).AvgProd(ExpressionBins{1,43}),AvgProdAllAP(43).All95Conf(ExpressionBins{1,43}),'o','Color','k','LineWidth',2.5);
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 2 5]);
ylabel('integrated fluorescence (AU)');
xticks([0.2]);
xticklabels('1x Stat92E');
xlim([0 0.4]);
ylim([3e5 1.6e6]);
%% Separate plots for boundary shifts 
% Set what percentile you want as exp cut offs 
ExpLim = 50;
SE_Bcd = [10,9,41,64,56];
APExpBins_SE = [ExpressionBins{1,[SE_Bcd]}];
DD_Bcd = [26,6,39,61,50];
APExpBins_DD = [ExpressionBins{1,[DD_Bcd]}];

figure
hold on 
% Find bounds - SE
% using threshold of 50% max exp of the homozygous construct
ExpThresh = (nanmax([AvgProdAllAP(9).AvgProd]))*0.5;
for bb = 1:length(SE_Bcd)
    ExpLvls = [AvgProdAllAP(SE_Bcd(bb)).AvgProd];
    % Ant bound
    Ant_Exp = ExpLvls(1:APExpBins_SE(bb)-1);
    Post_Exp = ExpLvls(APExpBins_SE(bb)+1:end);
    %gives the AP bin of left most that is at least % expression
    if bb == 2
        AntBound_SE(bb) = min(find(Ant_Exp >= ExpThresh));
    else
        AntBound_SE(bb) = min(find(Ant_Exp >= ExpThresh));
    end
    % AP bin right most that is at least % exp
    PostBound_SE(bb) = max(find(Post_Exp >= ExpThresh))+APExpBins_SE(bb);
    plot([EggLength(AntBound_SE(bb)) EggLength(PostBound_SE(bb))],[bb bb],'LineWidth',5,'Color',Colors(10).Color);
end

% Find bounds - 2xDist
% using threshold of 50% max exp of the homozygous construct
ExpThresh = (nanmax([AvgProdAllAP(6).AvgProd]))*0.5;
for bb = 1:length(DD_Bcd)
    ExpLvls = [AvgProdAllAP(DD_Bcd(bb)).AvgProd];
    % Ant bound
    Ant_Exp = ExpLvls(1:APExpBins_DD(bb)-1);
    Post_Exp = ExpLvls(APExpBins_DD(bb)+1:end);
    %gives the AP bin of left most that is at least % expression
    AntBound_DD(bb) = min(find(Ant_Exp >= ExpThresh));
    % AP bin right most that is at least % exp
    PostBound_DD(bb) = max(find(Post_Exp >= ExpThresh))+ APExpBins_DD(bb);
    plot([EggLength(AntBound_DD(bb)) EggLength(PostBound_DD(bb))],[bb+length(DD_Bcd) bb+length(DD_Bcd)],'LineWidth',5,'Color',Colors(6).Color);
end
xlim([20 80]);
ylim([0 bb+length(DD_Bcd)+1]);
xlabel('% egg length');
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);

ArrayColor = [179 16 182] ./255;
Array3Color =[155 16 152] ./255;
Array6Color = [115 16 107]./255;

%separate graphs for each construct
figure
SE_Colors = [0 0 0;Colors(9).Color; ArrayColor; Array3Color; Array6Color];
for bb = 1:length(SE_Bcd)
    plot([EggLength(AntBound_SE(bb)) EggLength(PostBound_SE(bb))],[bb bb],'LineWidth',5,'Color',SE_Colors(bb,:));
    hold on
end
ylim([0 5.5]);
xlim([20 80])
xlabel('% egg length');
yticks([1 2 3 4 5 6]);
yticklabels({'hemi', 'homo', '1xBcd', '3xBcd', '6xBcd'});
hline = line([EggLength(AntBound_SE(1)) EggLength(AntBound_SE(1))], [0 5.5]);
hline.Color ='k'; hline.LineStyle ='--';hline.LineWidth = 2;
hline_P = line([EggLength(PostBound_SE(1)) EggLength(PostBound_SE(1))], [0 5.5]);
hline_P.Color ='k'; hline_P.LineStyle ='--'; hline_P.LineWidth = 2;
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'SE_AntShifts','.pdf'],'pdf');


% Duplicated distal
figure
Array3Color =[155 16 152] ./255;
DD_Colors = [0 0 0;Colors(6).Color; ArrayColor; Array3Color; Array6Color];
for bb = 1:length(DD_Bcd)
    plot([EggLength(AntBound_DD(bb)) EggLength(PostBound_DD(bb))],[bb bb],'LineWidth',5,'Color',DD_Colors(bb,:));
    hold on
end
ylim([0 5.5]);
xlim([20 80])
xlabel('% egg length');
yticks([1 2 3 4 5 6]);
yticklabels({'hemi','homo', '1xBcd', '3xBcd', '6xBcd'});
hline = line([EggLength(AntBound_DD(1)) EggLength(AntBound_DD(1))], [0 5.5]);
hline.Color ='k'; hline.LineStyle ='--';hline.LineWidth = 2;
hline_P = line([EggLength(PostBound_DD(1)) EggLength(PostBound_DD(1))], [0 5.5]);
hline_P.Color ='k'; hline_P.LineStyle ='--'; hline_P.LineWidth = 2;
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep '2xDist_AntShifts','.pdf'],'pdf');

% Each construct, peak exp + boundaries
figure 
SE_Colors = [0 0 0;Colors(9).Color; ArrayColor; Array6Color];
for bb = 1:length(SE_Bcd)
    plot([EggLength(AntBound_SE(bb)) EggLength(PostBound_SE(bb))],[nanmax(AvgProdAllAP(SE_Bcd(bb)).AvgProd) nanmax(AvgProdAllAP(SE_Bcd(bb)).AvgProd)],'LineWidth',5,'Color',SE_Colors(bb,:));
    hold on
end
xlim([20 80]);
xlabel('% egg length');
ylabel('peak mRNA expression');
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);

%% Bootstrapping error bars for boundary shift
%use 1/2 max of homozygous construct
CutoffLim = nanmax([AvgProdAllAP(6).AvgProd]) *0.5;%input('What percentile to use?');
Ntimes = 1000; % set the number of times to try bootstrapping
counter=0;
for cc = [26,6,39,61,50]; % set constructs to do this for
    counter = counter+1;
    ConData = [AvgProdAllAP(cc).AllProds]; % sample from this
    [TrueAntBound, TruePostBound] = CalcExpBounds(ConData,CutoffLim); %calc from orig data
    BootBound_All=[];
    %bootstrap error in exp bound shifts
    for tt = 1:Ntimes
            y = datasample(ConData,size(ConData,1),1); %Randomly sample rows (columns are AP bins)
            [BootBound_Ant, BootBound_Post] = CalcExpBounds(y, CutoffLim); %The fx we want performed every time
            BootBound_Indiv = [BootBound_Ant, BootBound_Post];
            BootBound_All = [BootBound_All; BootBound_Indiv]; %Get all samplings into one matrix 
        
    end
     BootBoundStd(counter,:) = nanstd(BootBound_All); %This is the SE from bootstrapping!
     BootBoundUp95(counter,:) = prctile(BootBound_All,97.5); %Upper CI estimates
     BootBoundLow95(counter,:) = prctile(BootBound_All,2.5); %Lower CI estimates
     % is this the correct way to estimate CI ???
     BootBound95Try(counter,:) = BootBoundStd(counter,:) .* 1.96;
end

% Shadow pair
CutoffLim = nanmax([AvgProdAllAP(9).AvgProd]) *0.5;
counter=0;
for cc = [10,9,41,64,56]; % set constructs to do this for
    counter = counter+1;
    ConData = [AvgProdAllAP(cc).AllProds]; % sample from this
    [TrueAntBound, TruePostBound] = CalcExpBounds(ConData,CutoffLim); %calc from orig data
    BootBound_All=[];
    %bootstrap error in exp bound shifts
    for tt = 1:Ntimes
            y = datasample(ConData,size(ConData,1),1); %Randomly sample rows (columns are AP bins)
            [BootBound_Ant, BootBound_Post] = CalcExpBounds(y, CutoffLim); %The fx we want performed every time
            BootBound_Indiv = [BootBound_Ant, BootBound_Post];
            BootBound_All = [BootBound_All; BootBound_Indiv]; %Get all samplings into one matrix 
        
    end
     BootBoundStd_SE(counter,:) = nanstd(BootBound_All); %This is the SE from bootstrapping!
     BootBoundUp95_SE(counter,:) = prctile(BootBound_All,97.5); %Upper CI estimates
     BootBoundLow95_SE(counter,:) = prctile(BootBound_All,2.5); %Lower CI estimates
     % is this the correct way to estimate CI ???
     BootBound95Try_SE(counter,:) = BootBoundStd_SE(counter,:) .* 1.96;
end
%%
% Duplicated distal
x_sizeNarrow =3.7;
y_sizeNarrow =5.52; 
figure
Grey =[0.5 0.5 0.5];
Array3Color =[155 16 152] ./255;
% Set transparency level (0:1) for errorbars
alpha = 0.65; 
DD_Colors = [0 0 0;Colors(6).Color; ArrayColor; Array3Color; Array6Color];
for bb = 1:length(DD_Bcd)
    plot([EggLength(AntBound_DD(bb)) EggLength(PostBound_DD(bb))],[bb bb],'LineWidth',5,'Color',DD_Colors(bb,:));
    hold on
    h = errorbar(EggLength(AntBound_DD(bb)),bb,BootBoundStd(bb,1),'horizontal','Color', 'k','LineWidth',2.5);
      
    % Set transparency (undocumented)
    set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h.Line.ColorData(1:3); 255*alpha])
    
    h2= errorbar(EggLength(PostBound_DD(bb)),bb,BootBoundStd(bb,2),'horizontal','Color','k','LineWidth', 2.5);
    set([h2.Bar, h2.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h2.Line.ColorData(1:3); 255*alpha])

end
ylim([0 6]);
xlim([30 70])
xlabel('% egg length');
yticks([1 2 3 4 5 6]);
yticklabels({'hemi','homo', '1xBcd', '3xBcd', '6xBcd'});
hline = line([EggLength(AntBound_DD(2)) EggLength(AntBound_DD(2))], [0 6]);
hline.Color ='k'; hline.LineStyle ='--';hline.LineWidth = 2;
hline_P = line([EggLength(PostBound_DD(2)) EggLength(PostBound_DD(2))], [0 6]);
hline_P.Color ='k'; hline_P.LineStyle ='--'; hline_P.LineWidth = 2;
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop x_sizeNarrow y_sizeNarrow]);
saveas(gcf, [FigDirect filesep '2xDist_AntShifts_ErrorBars','.pdf'],'pdf');

% shadow pair
figure
Grey =[0.5 0.5 0.5];
Array3Color =[155 16 152] ./255;
SE_Colors = [Colors(9).Color; Colors(9).Color; ArrayColor; Array3Color; Array6Color];
for bb = 1:length(SE_Bcd)
    plot([EggLength(AntBound_SE(bb)) EggLength(PostBound_SE(bb))],[bb bb],'LineWidth',5,'Color',SE_Colors(bb,:));
    hold on
    errorbar(EggLength(AntBound_SE(bb)),bb,BootBoundStd(bb,1),'horizontal','Color', 'k','LineWidth',2.5);
    errorbar(EggLength(PostBound_SE(bb)),bb,BootBoundStd(bb,2),'horizontal','Color','k','LineWidth', 2.5);
end
ylim([0 6]);
xlim([30 70])
xlabel('% egg length');
yticks([1 2 3 4 5 6]);
yticklabels({'hemi','homo', '1xBcd', '3xBcd', '6xBcd'});
hline = line([EggLength(AntBound_SE(2)) EggLength(AntBound_SE(2))], [0 6]);
hline.Color ='k'; hline.LineStyle ='--';hline.LineWidth = 2;
hline_P = line([EggLength(PostBound_SE(2)) EggLength(PostBound_SE(2))], [0 6]);
hline_P.Color ='k'; hline_P.LineStyle ='--'; hline_P.LineWidth = 2;
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop x_sizeNarrow y_sizeNarrow]);
saveas(gcf, [FigDirect filesep 'SE_AntShifts_ErrorBars','.pdf'],'pdf');


%% 4F saturating MCP-GFP hemi vs homo
figure
errorbar(EggLength, AvgProdAllAP(62).AvgProd, AvgProdAllAP(62).All95Conf,'Color','k','LineStyle','--','LineWidth', 2.5);
hold on 
errorbar(EggLength, AvgProdAllAP(65).AvgProd, AvgProdAllAP(65).All95Conf,'Color','k','LineWidth', 2.5);
xlabel('% egg length');
ylabel('integrated fluorescence (AU)');
xlim([20 80]);
ylim([0 16e5]);
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
%% HbP2 - MS2 competition assessments  
figure 
errorbar(EggLength, AvgProdAllAP(16).AvgProd, AvgProdAllAP(16).All95Conf, 'Color', 'k','LineWidth',2.5);
hold on 
errorbar(EggLength, AvgProdAllAP(57).AvgProd, AvgProdAllAP(57).All95Conf, 'Color', grey,'LineWidth',2.5);
set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 %title('Avg mRNA production per nucleus') 
 ylim([0 1000000]);
 xlim([15 60]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(16).AvgProd./MS2Conversion),(AvgProdAllAP(16).All95Conf./MS2Conversion),'Color','k','LineWidth',3.5,'LineStyle','-');
errorbar(EggLength, (AvgProdAllAP(57).AvgProd./MS2Conversion),(AvgProdAllAP(57).All95Conf./MS2Conversion),'Color',grey,'LineWidth',3.5,'LineStyle',':');
ylabel('total transcripts produced');
ylim([0 (1000000/MS2Conversion)]);
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'HbP2_HemivHomo_TotOut','.pdf'],'pdf');

% HbP2 w 6xBcd
ArrayColor = [179 16 182] ./255;
Array6Color = [115 16 107]./255;
figure

lineProps.col{1} = [1 1 1];
mseb(EggLength, AvgProdAllAP(16).AvgProd, AvgProdAllAP(16).All95Conf, lineProps,0)
  hold on 
lineProps.col{1} = [0.5 0.5 0.5];
mseb(EggLength, AvgProdAllAP(57).AvgProd, AvgProdAllAP(57).All95Conf, lineProps,0);
%lineProps.col{1} = ArrayColor;
%mseb(EggLength, AvgProdAllAP(41).AvgProd, AvgProdAllAP(41).All95Conf, lineProps,0);
lineProps.col{1}=Array6Color;
mseb(EggLength, AvgProdAllAP(59).AvgProd, AvgProdAllAP(59).All95Conf, lineProps,0);
set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 %title('Avg mRNA production per nucleus') 
 ylim([0 1400000]);
 xlim([0 100]);
 yyaxis right
 lineProps.col{1} = Colors(16).Color
mseb(EggLength, (AvgProdAllAP(16).AvgProd./MS2Conversion), (AvgProdAllAP(16).All95Conf./MS2Conversion), lineProps,0)
  hold on 
lineProps.col{1} = [0.5 0.5 0.5];  
mseb(EggLength, (AvgProdAllAP(57).AvgProd./MS2Conversion), (AvgProdAllAP(57).All95Conf./MS2Conversion), lineProps,0);
%lineProps.col{1} = ArrayColor;
%mseb(EggLength, (AvgProdAllAP(41).AvgProd./MS2Conversion), (AvgProdAllAP(41).All95Conf./MS2Conversion), lineProps,0);
lineProps.col{1}=Array6Color;
mseb(EggLength, (AvgProdAllAP(59).AvgProd./MS2Conversion), (AvgProdAllAP(59).All95Conf./MS2Conversion), lineProps,0);

ylabel('total transcripts produced');
ylim([0 (1400000/MS2Conversion)]);
xlim([20 80]);
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'HbP2_6xBcdArray_TotOut','.pdf'],'pdf');

% Homo as fx of hemi graph
 figure
lineProps.col{1} = [0.5 0.5 0.5];
mseb([-9:10], [ones(1,20).*AvgProdAllAP(57).AvgProd(ExpressionBins{1,57})], ones(1,20).*AvgProdAllAP(57).All95Conf(ExpressionBins{1,57}), lineProps,0)
hold on 
mseb([-9:10], [ones(1,20).*AvgProdAllAP(16).AvgProd(ExpressionBins{1,16})], ones(1,20).*AvgProdAllAP(16).All95Conf(ExpressionBins{1,16}), lineProps,0)

errorbar(1, AvgProdAllAP(59).AvgProd(ExpressionBins{1,59}), AvgProdAllAP(59).All95Conf(ExpressionBins{1,59}),'o','LineWidth',2.5,'Color','k');

xlim([0 3]);
ylim([3e5 1.6e6]);
ylabel('integrated fluorescence (AU)');
xlabel('# of Bcd binding arrays')
xticks([1]);
xticklabels({'6x'});
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'HbP2_BcdArrays_HemiHomoLines_Shade','.pdf'],'pdf');



%% Compare production Chr2 vs Chr3 
% 2xDistal 
figure
 errorbar(EggLength,AvgProdAllAP(6).AvgProd,AvgProdAllAP(6).All95Conf,'Color',Colors(6).Color,'LineWidth',3.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(36).AvgProd, AvgProdAllAP(36).All95Conf,'Color',Colors(36).Color,'LineWidth',3.5,'LineStyle','--');
 %errorbar(EggLength, AvgProdAllAP(51).AvgProd, AvgProdAllAP(51).All95Conf, 'Color', 'k', 'LineWidth',3.5,'LineStyle',"--");
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 %title('Avg mRNA production per nucleus') 
 ylim([0 1400000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(6).AvgProd./MS2Conversion),(AvgProdAllAP(6).All95Conf./MS2Conversion),'Color',Colors(6).Color,'LineWidth',3.5,'LineStyle','-');
errorbar(EggLength, (AvgProdAllAP(36).AvgProd./MS2Conversion),(AvgProdAllAP(36).All95Conf./MS2Conversion),'Color',Colors(36).Color,'LineWidth',3.5,'LineStyle','--');
ylabel('total transcripts produced');
ylim([0 (1400000/MS2Conversion)]);
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
%saveas(gcf, [FigDirect filesep 'TotProdKrSE_Chr3_1xBcd','.pdf'],'pdf');

%Hemi vs Homo on ChrIII - 2xDistal
  figure
 errorbar(EggLength,AvgProdAllAP(6).AvgProd,AvgProdAllAP(6).All95Conf,'Color',Colors(6).Color,'LineWidth',3.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(38).AvgProd, AvgProdAllAP(38).All95Conf,'Color',Colors(38).Color,'LineWidth',3.5,'LineStyle','--');
 %errorbar(EggLength, AvgProdAllAP(51).AvgProd, AvgProdAllAP(51).All95Conf, 'Color', 'k', 'LineWidth',3.5,'LineStyle',"--");
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 %title('Avg mRNA production per nucleus') 
 ylim([0 1400000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(6).AvgProd./MS2Conversion),(AvgProdAllAP(6).All95Conf./MS2Conversion),'Color',Colors(6).Color,'LineWidth',3.5,'LineStyle','-');
errorbar(EggLength, (AvgProdAllAP(36).AvgProd./MS2Conversion),(AvgProdAllAP(36).All95Conf./MS2Conversion),'Color',Colors(36).Color,'LineWidth',3.5,'LineStyle','--');
ylabel('total transcripts produced');
ylim([0 (1400000/MS2Conversion)]);
title('Homo vs hemi ChrIII');
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);

%Shadow pair
figure
 errorbar(EggLength,AvgProdAllAP(9).AvgProd,AvgProdAllAP(9).All95Conf,'Color',Colors(9).Color,'LineWidth',3.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(35).AvgProd, AvgProdAllAP(35).All95Conf,'Color',Colors(35).Color,'LineWidth',3.5,'LineStyle','--');
 %errorbar(EggLength, AvgProdAllAP(51).AvgProd, AvgProdAllAP(51).All95Conf, 'Color', 'k', 'LineWidth',3.5,'LineStyle',"--");
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 %title('Avg mRNA production per nucleus') 
 ylim([0 1400000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(9).AvgProd./MS2Conversion),(AvgProdAllAP(9).All95Conf./MS2Conversion),'Color',Colors(9).Color,'LineWidth',3.5,'LineStyle','-');
errorbar(EggLength, (AvgProdAllAP(35).AvgProd./MS2Conversion),(AvgProdAllAP(35).All95Conf./MS2Conversion),'Color',Colors(35).Color,'LineWidth',3.5,'LineStyle','--');
ylabel('total transcripts produced');
ylim([0 (1400000/MS2Conversion)]);
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
%saveas(gcf, [FigDirect filesep 'TotProdKrSE_Chr3_1xBcd','.pdf'],'pdf');

  % Homo vs Hemi - Shadow pair
  figure
 errorbar(EggLength,AvgProdAllAP(9).AvgProd,AvgProdAllAP(9).All95Conf,'Color',Colors(9).Color,'LineWidth',3.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(40).AvgProd, AvgProdAllAP(40).All95Conf,'Color',Colors(40).Color,'LineWidth',3.5,'LineStyle','--');
 %errorbar(EggLength, AvgProdAllAP(51).AvgProd, AvgProdAllAP(51).All95Conf, 'Color', 'k', 'LineWidth',3.5,'LineStyle',"--");
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 %title('Avg mRNA production per nucleus') 
 ylim([0 1400000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(9).AvgProd./MS2Conversion),(AvgProdAllAP(9).All95Conf./MS2Conversion),'Color',Colors(9).Color,'LineWidth',3.5,'LineStyle','-');
errorbar(EggLength, (AvgProdAllAP(40).AvgProd./MS2Conversion),(AvgProdAllAP(40).All95Conf./MS2Conversion),'Color',Colors(40).Color,'LineWidth',3.5,'LineStyle','--');
ylabel('total transcripts produced');
ylim([0 (1400000/MS2Conversion)]);
title('Homo vs Hemi ChrIII')
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);

%2xProximal 
figure
 errorbar(EggLength,AvgProdAllAP(7).AvgProd,AvgProdAllAP(7).All95Conf,'Color',Colors(7).Color,'LineWidth',3.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(37).AvgProd, AvgProdAllAP(37).All95Conf,'Color',Colors(37).Color,'LineWidth',3.5,'LineStyle','--');
 %errorbar(EggLength, AvgProdAllAP(51).AvgProd, AvgProdAllAP(51).All95Conf, 'Color', 'k', 'LineWidth',3.5,'LineStyle',"--");
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 %title('Avg mRNA production per nucleus') 
 ylim([0 1400000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(7).AvgProd./MS2Conversion),(AvgProdAllAP(7).All95Conf./MS2Conversion),'Color',Colors(7).Color,'LineWidth',3.5,'LineStyle','-');
errorbar(EggLength, (AvgProdAllAP(37).AvgProd./MS2Conversion),(AvgProdAllAP(37).All95Conf./MS2Conversion),'Color',Colors(37).Color,'LineWidth',3.5,'LineStyle','--');
ylabel('total transcripts produced');
ylim([0 (1400000/MS2Conversion)]);
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
%saveas(gcf, [FigDirect filesep 'TotProdKrSE_Chr3_1xBcd','.pdf'],'pdf');
%% TFBS with Chr3 constructs
% SE w 1xBcd array
figure
 errorbar(EggLength,AvgProdAllAP(35).AvgProd,AvgProdAllAP(35).All95Conf,'Color',Colors(35).Color,'LineWidth',3.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(40).AvgProd, AvgProdAllAP(40).All95Conf,'Color',Colors(40).Color,'LineWidth',3.5,'LineStyle',':');
 errorbar(EggLength, AvgProdAllAP(51).AvgProd, AvgProdAllAP(51).All95Conf, 'Color', 'k', 'LineWidth',3.5,'LineStyle',"--");
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 %title('Avg mRNA production per nucleus') 
 ylim([0 1400000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(35).AvgProd./MS2Conversion),(AvgProdAllAP(35).All95Conf./MS2Conversion),'Color',Colors(35).Color,'LineWidth',3.5,'LineStyle','-');
errorbar(EggLength, (AvgProdAllAP(40).AvgProd./MS2Conversion),(AvgProdAllAP(40).All95Conf./MS2Conversion),'Color',Colors(40).Color,'LineWidth',3.5,'LineStyle',':');
errorbar(EggLength, (AvgProdAllAP(51).AvgProd./MS2Conversion),(AvgProdAllAP(51).All95Conf./MS2Conversion),'Color','k','LineWidth',3.5,'LineStyle','--');
ylabel('total transcripts produced');
ylim([0 (1400000/MS2Conversion)]);
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'TotProdKrSE_Chr3_1xBcd','.pdf'],'pdf');
%% Percent decrease in mRNA production homo vs hemi 
cc=input('Homozygous construct?');
cc2=input('Hemizygous construct?');
figure
mRNAChange=[AvgProdAllAP(cc2).AvgProd - AvgProdAllAP(cc).AvgProd];
mRNAChange=[mRNAChange./(AvgProdAllAP(cc2).AvgProd)];
plot(EggLength,(mRNAChange*100),'LineWidth',2.5,'Color',Colors(cc).Color);
hold on 
Homopeak=ExpressionBins{1,cc}; HemiPeak=ExpressionBins{1,cc2};
h=line([EggLength(Homopeak) EggLength(Homopeak)], [-20 20]);
h2=line([EggLength(HemiPeak) EggLength(HemiPeak)], [-20 20],'LineStyle','--'); 


figure
TopmRNAChange=[mRNAChange(ExpressionBins{1,cc}) mRNAChange(ExpressionBins{1,cc2})];
bar(TopmRNAChange.*100)
xticklabels({'Homo peak', 'Hemi peak'})
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('% decrease per allele expression');
title(ConstructList{cc});


%% Temp comparision hemizygotes
figure
 errorbar(EggLength,AvgProdAllAP(26).AvgProd,AvgProdAllAP(26).All95Conf,'Color',Colors(26).Color,'LineWidth',2.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(27).AvgProd, AvgProdAllAP(27).All95Conf,'Color',Colors(27).Color,'LineWidth',2.5,'LineStyle',':');
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 %title('Avg mRNA production per nucleus') 
 ylim([0 1400000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(26).AvgProd./MS2Conversion),(AvgProdAllAP(26).All95Conf./MS2Conversion),'Color',Colors(26).Color,'LineWidth',2.5,'LineStyle','-');
  errorbar(EggLength, (AvgProdAllAP(27).AvgProd./MS2Conversion),(AvgProdAllAP(27).All95Conf./MS2Conversion),'Color',Colors(27).Color,'LineWidth',2.5,'LineStyle',':');
ylabel('total transcripts produced');
ylim([0 (1400000/MS2Conversion)]);
set(gca, 'YColor','k');
%% 
  figure
 errorbar(EggLength,AvgProdAllAP(4).AvgProd,AvgProdAllAP(4).All95Conf,'Color',Colors(4).Color,'LineWidth',2.5,'LineStyle',':');
 hold on 
 errorbar(EggLength,AvgProdAllAP(9).AvgProd, AvgProdAllAP(9).All95Conf,'Color',Colors(9).Color,'LineWidth',2.5);
 legend('1x Distal','Both');
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 %title('Avg mRNA production per nucleus') 
 ylim([0 1100000]);
 xlim([0 100]);
 
 figure
 errorbar(EggLength,AvgProdAllAP(5).AvgProd,AvgProdAllAP(5).All95Conf,'Color',Colors(5).Color,'LineWidth',2.5,'LineStyle',':');
 hold on 
 errorbar(EggLength,AvgProdAllAP(9).AvgProd, AvgProdAllAP(9).All95Conf,'Color',Colors(9).Color,'LineWidth',2.5);
 legend('1x Proximal','Both');
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 %title('Avg mRNA production per nucleus') 
 ylim([0 1100000]);
 xlim([0 100]);
 
 
 
  
 
  figure
 errorbar(EggLength,AvgProdAllAP(5).AvgProd,AvgProdAllAP(5).All95Conf,'Color',Colors(5).Color,'LineWidth',2.5,'LineStyle',':');
 hold on 
 errorbar(EggLength,AvgProdAllAP(8).AvgProd,AvgProdAllAP(8).All95Conf,'Color',Colors(8).Color,'LineWidth',2.5,'LineStyle',':');
 legend('1x Proximal','Single 2xProximal');
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 ylim([0 1100000]);
 xlim([0 100]);
 
   figure
 errorbar(EggLength,AvgProdAllAP(4).AvgProd,AvgProdAllAP(4).All95Conf,'Color',Colors(4).Color,'LineWidth',2.5,'LineStyle',':');
 hold on 
 errorbar(EggLength,AvgProdAllAP(5).AvgProd,AvgProdAllAP(5).All95Conf,'Color',Colors(5).Color,'LineWidth',2.5,'LineStyle',':');
 errorbar(EggLength,AvgProdAllAP(10).AvgProd, AvgProdAllAP(10).All95Conf,'Color',Colors(10).Color,'LineWidth',2.5,'LineStyle',':');
 %legend('1x Distal','1x Proximal','1x Both');
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 ylim([0 1100000]);
 xlim([0 100]);

 %%
 %Compare SE and duplicates
 figure 
 hold on 
 errorbar(EggLength, AvgProdAllAP(7).AvgProd, AvgProdAllAP(7).All95Conf,'Color',Colors(7).Color,'LineWidth',2.5);
errorbar(EggLength, AvgProdAllAP(9).AvgProd,AvgProdAllAP(9).All95Conf,'Color',Colors(9).Color,'LineWidth',2.5);
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence')
 xlim([0 100]);
 ylim([0 900000]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(7).AvgProd./MS2Conversion),(AvgProdAllAP(7).All95Conf./MS2Conversion),'Color',Colors(7).Color,'LineWidth',2.5,'LineStyle','-');
  errorbar(EggLength, (AvgProdAllAP(9).AvgProd./MS2Conversion),(AvgProdAllAP(9).All95Conf./MS2Conversion),'Color',Colors(9).Color,'LineWidth',2.5,'LineStyle','-');
ylabel('transcripts produced');
ylim([0 (900000/MS2Conversion)]);
set(gca,'FontSize',fontsize,'FontName',fontname);
 print( [FigDirect filesep '2ProxBothTotalmRNA'],'-dsvg');
%legend('2x Proximal', 'Both','Location','best');
 %title('Avg mRNA production per nucleus') 
 
 
 
 %% Endogenous Distal RT vs 32C
 figure
 errorbar(EggLength,AvgProdAllAP(32).AvgProd,AvgProdAllAP(32).All95Conf,'Color',Colors(32).Color,'LineWidth',2.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(33).AvgProd, AvgProdAllAP(33).All95Conf,'Color',Colors(33).Color,'LineWidth',2.5,'LineStyle','-.');
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence')
 %title('Avg mRNA production per nucleus') 
 ylim([0 950000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(32).AvgProd./MS2Conversion),(AvgProdAllAP(32).All95Conf./MS2Conversion),'Color',Colors(32).Color,'LineWidth',2.5,'LineStyle','-');
  errorbar(EggLength, (AvgProdAllAP(33).AvgProd./MS2Conversion),(AvgProdAllAP(33).All95Conf./MS2Conversion),'Color',Colors(33).Color,'LineWidth',2.5,'LineStyle','-.');
ylabel('transcripts produced');
ylim([0 (950000/MS2Conversion)]);
%legend('Distal', 'Distal 32C','Location','best');

%17C
figure
 errorbar(EggLength,AvgProdAllAP(32).AvgProd,AvgProdAllAP(32).All95Conf,'Color',Colors(32).Color,'LineWidth',2.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(34).AvgProd, AvgProdAllAP(34).All95Conf,'Color',Colors(34).Color,'LineWidth',2.5,'LineStyle','-.');
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence')
 %title('Avg mRNA production per nucleus') 
 ylim([0 970000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(32).AvgProd./MS2Conversion),(AvgProdAllAP(32).All95Conf./MS2Conversion),'Color',Colors(32).Color,'LineWidth',2.5,'LineStyle','-');
  errorbar(EggLength, (AvgProdAllAP(34).AvgProd./MS2Conversion),(AvgProdAllAP(34).All95Conf./MS2Conversion),'Color',Colors(34).Color,'LineWidth',2.5,'LineStyle','-.');
ylabel('transcripts produced');
ylim([0 (970000/MS2Conversion)]);
%legend('Distal', 'Distal 32C','Location','best');
 %%
 %32C vs RT
 figure
 errorbar(EggLength,AvgProdAllAP(1).AvgProd,AvgProdAllAP(1).All95Conf,'Color',Colors(1).Color,'LineWidth',2.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(11).AvgProd, AvgProdAllAP(11).All95Conf,'Color',Colors(11).Color,'LineWidth',2.5,'LineStyle','-.');
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence')
 %title('Avg mRNA production per nucleus') 
 ylim([0 950000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(1).AvgProd./MS2Conversion),(AvgProdAllAP(1).All95Conf./MS2Conversion),'Color',Colors(1).Color,'LineWidth',2.5,'LineStyle','-');
  errorbar(EggLength, (AvgProdAllAP(11).AvgProd./MS2Conversion),(AvgProdAllAP(11).All95Conf./MS2Conversion),'Color',Colors(11).Color,'LineWidth',2.5,'LineStyle','-.');
ylabel('transcripts produced');
ylim([0 (950000/MS2Conversion)]);
%legend('Distal', 'Distal 32C','Location','best');
print( [FigDirect filesep 'DistTempCompTotalmRNA'],'-dsvg');


figure
 errorbar(EggLength,AvgProdAllAP(2).AvgProd,AvgProdAllAP(2).All95Conf,'Color',Colors(2).Color,'LineWidth',2.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(12).AvgProd, AvgProdAllAP(12).All95Conf,'Color',Colors(12).Color,'LineWidth',2.5,'LineStyle','-.');
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence')
 %title('Avg mRNA production per nucleus') 
 ylim([0 950000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(2).AvgProd./MS2Conversion),(AvgProdAllAP(2).All95Conf./MS2Conversion),'Color',Colors(2).Color,'LineWidth',2.5,'LineStyle','-');
  errorbar(EggLength, (AvgProdAllAP(12).AvgProd./MS2Conversion),(AvgProdAllAP(12).All95Conf./MS2Conversion),'Color',Colors(12).Color,'LineWidth',2.5,'LineStyle','--');
ylabel('transcripts produced');
ylim([0 (950000/MS2Conversion)]);
%legend('Proximal', 'Proximal 32C','Location','best');
 print( [FigDirect filesep 'ProxTempTotalmRNA'],'-dsvg');


 
 figure
 errorbar(EggLength,AvgProdAllAP(3).AvgProd,AvgProdAllAP(3).All95Conf,'Color',Colors(3).Color,'LineWidth',2.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(13).AvgProd, AvgProdAllAP(12).All95Conf,'Color',Colors(13).Color,'LineWidth',2.5,'LineStyle','-.');
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence')
 ylim([0 950000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(3).AvgProd./MS2Conversion),(AvgProdAllAP(3).All95Conf./MS2Conversion),'Color',Colors(3).Color,'LineWidth',2.5,'LineStyle','-');
  errorbar(EggLength, (AvgProdAllAP(13).AvgProd./MS2Conversion),(AvgProdAllAP(12).All95Conf./MS2Conversion),'Color',Colors(13).Color,'LineWidth',2.5,'LineStyle','-.');
ylabel('transcripts produced');
ylim([0 (950000/MS2Conversion)]);
%legend('Both Sep', 'Both Sep 32C','Location','best');
 print( [FigDirect filesep 'SepTempTotalmRNA'],'-dsvg');

 
figure
 errorbar(EggLength,AvgProdAllAP(9).AvgProd,AvgProdAllAP(9).All95Conf,'Color',Colors(9).Color,'LineWidth',2.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(14).AvgProd, AvgProdAllAP(13).All95Conf,'Color',Colors(14).Color,'LineWidth',2.5,'LineStyle','-.');
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 %title('Avg mRNA production per nucleus') 
 ylim([0 950000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(9).AvgProd./MS2Conversion),(AvgProdAllAP(9).All95Conf./MS2Conversion),'Color',Colors(9).Color,'LineWidth',2.5,'LineStyle','-');
  errorbar(EggLength, (AvgProdAllAP(14).AvgProd./MS2Conversion),(AvgProdAllAP(14).All95Conf./MS2Conversion),'Color',Colors(14).Color,'LineWidth',2.5,'LineStyle','-.');
ylabel('transcripts produced');
ylim([0 (950000/MS2Conversion)]);
%legend('Both', 'Both 32C','Location','best');
print( [FigDirect filesep 'BothTempTotalmRNA'],'-dsvg');

 
 figure
 errorbar(EggLength,AvgProdAllAP(7).AvgProd,AvgProdAllAP(7).All95Conf,'Color',Colors(7).Color,'LineWidth',2.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(15).AvgProd, AvgProdAllAP(15).All95Conf,'Color',Colors(15).Color,'LineWidth',2.5,'LineStyle','-.');
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 %title('Avg mRNA production per nucleus') 
 ylim([0 950000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(7).AvgProd./MS2Conversion),(AvgProdAllAP(7).All95Conf./MS2Conversion),'Color',Colors(7).Color,'LineWidth',2.5,'LineStyle','-');
  errorbar(EggLength, (AvgProdAllAP(15).AvgProd./MS2Conversion),(AvgProdAllAP(15).All95Conf./MS2Conversion),'Color',Colors(15).Color,'LineWidth',2.5,'LineStyle','-.');
ylabel('transcripts produced');
ylim([0 (950000/MS2Conversion)]);
%legend('2x Prox', '2x Prox 32C','Location','best');
 print( [FigDirect filesep '2xProxTempTotalmRNA'],'-dsvg');
 
 figure
 errorbar(EggLength,AvgProdAllAP(6).AvgProd,AvgProdAllAP(6).All95Conf,'Color',Colors(6).Color,'LineWidth',2.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(19).AvgProd, AvgProdAllAP(19).All95Conf,'Color',Colors(19).Color,'LineWidth',2.5,'LineStyle','-.');
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence')
 %title('Avg mRNA production per nucleus') 
 ylim([0 1150000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(6).AvgProd./MS2Conversion),(AvgProdAllAP(6).All95Conf./MS2Conversion),'Color',Colors(6).Color,'LineWidth',2.5,'LineStyle','-');
  errorbar(EggLength, (AvgProdAllAP(19).AvgProd./MS2Conversion),(AvgProdAllAP(19).All95Conf./MS2Conversion),'Color',Colors(19).Color,'LineWidth',2.5,'LineStyle','-.');
ylabel('transcripts produced');
ylim([0 (1150000/MS2Conversion)]);
%legend('Distal', 'Distal 32C','Location','best');
print( [FigDirect filesep '2xDistTempCompTotalmRNA'],'-dsvg');
 
 %% 17C vs RT
 figure
 errorbar(EggLength,AvgProdAllAP(1).AvgProd,AvgProdAllAP(1).All95Conf,'Color',Colors(1).Color,'LineWidth',2.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(18).AvgProd, AvgProdAllAP(18).All95Conf,'Color',Colors(18).Color,'LineWidth',2.5,'LineStyle','-.');
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence')
 %title('Avg mRNA production per nucleus') 
 ylim([0 950000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(1).AvgProd./MS2Conversion),(AvgProdAllAP(1).All95Conf./MS2Conversion),'Color',Colors(1).Color,'LineWidth',2.5,'LineStyle','-');
  errorbar(EggLength, (AvgProdAllAP(18).AvgProd./MS2Conversion),(AvgProdAllAP(18).All95Conf./MS2Conversion),'Color',Colors(18).Color,'LineWidth',2.5,'LineStyle','-.');
ylabel('transcripts produced','Color','k');
set(gca,'ycolor','k');
ylim([0 (950000/MS2Conversion)]);
%legend('Distal', 'Distal 32C','Location','best');
print( [FigDirect filesep 'DistC17CompTotalmRNA'],'-dsvg');

 figure
 errorbar(EggLength,AvgProdAllAP(2).AvgProd,AvgProdAllAP(2).All95Conf,'Color',Colors(2).Color,'LineWidth',2.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(22).AvgProd, AvgProdAllAP(22).All95Conf,'Color',Colors(22).Color,'LineWidth',2.5,'LineStyle','-.');
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence')
 %title('Avg mRNA production per nucleus') 
 ylim([0 950000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(2).AvgProd./MS2Conversion),(AvgProdAllAP(2).All95Conf./MS2Conversion),'Color',Colors(2).Color,'LineWidth',2.5,'LineStyle','-');
  errorbar(EggLength, (AvgProdAllAP(22).AvgProd./MS2Conversion),(AvgProdAllAP(22).All95Conf./MS2Conversion),'Color',Colors(22).Color,'LineWidth',2.5,'LineStyle','-.');
ylabel('transcripts produced','Color','k');
set(gca,'ycolor','k');
ylim([0 (950000/MS2Conversion)]);
%legend('Distal', 'Distal 32C','Location','best');
print( [FigDirect filesep 'ProxC17CompTotalmRNA'],'-dsvg');

figure
 errorbar(EggLength,AvgProdAllAP(9).AvgProd,AvgProdAllAP(9).All95Conf,'Color',Colors(9).Color,'LineWidth',2.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(20).AvgProd, AvgProdAllAP(20).All95Conf,'Color',Colors(20).Color,'LineWidth',2.5,'LineStyle','-.');
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence')
 %title('Avg mRNA production per nucleus') 
 ylim([0 1100000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(9).AvgProd./MS2Conversion),(AvgProdAllAP(9).All95Conf./MS2Conversion),'Color',Colors(9).Color,'LineWidth',2.5,'LineStyle','-');
  errorbar(EggLength, (AvgProdAllAP(20).AvgProd./MS2Conversion),(AvgProdAllAP(20).All95Conf./MS2Conversion),'Color',Colors(20).Color,'LineWidth',2.5,'LineStyle','-.');
ylabel('transcripts produced');
ylim([0 (1100000/MS2Conversion)]);
set(gca, 'ycolor','k')
%legend('Distal', 'Distal 32C','Location','best');
print( [FigDirect filesep 'BothC17CompTotalmRNA'],'-dsvg');

figure
 errorbar(EggLength,AvgProdAllAP(6).AvgProd,AvgProdAllAP(6).All95Conf,'Color',Colors(6).Color,'LineWidth',2.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(21).AvgProd, AvgProdAllAP(21).All95Conf,'Color',Colors(21).Color,'LineWidth',2.5,'LineStyle','-.');
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence')
 %title('Avg mRNA production per nucleus') 
 ylim([0 950000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(6).AvgProd./MS2Conversion),(AvgProdAllAP(6).All95Conf./MS2Conversion),'Color',Colors(6).Color,'LineWidth',2.5,'LineStyle','-');
  errorbar(EggLength, (AvgProdAllAP(21).AvgProd./MS2Conversion),(AvgProdAllAP(21).All95Conf./MS2Conversion),'Color',Colors(21).Color,'LineWidth',2.5,'LineStyle','-.');
ylabel('transcripts produced','Color','k');
set(gca,'ycolor','k');
ylim([0 (950000/MS2Conversion)]);
%legend('Distal', 'Distal 32C','Location','best');
print( [FigDirect filesep '2xDistC17CompTotalmRNA'],'-dsvg');

figure
 errorbar(EggLength,AvgProdAllAP(7).AvgProd,AvgProdAllAP(7).All95Conf,'Color',Colors(7).Color,'LineWidth',2.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(25).AvgProd, AvgProdAllAP(25).All95Conf,'Color',Colors(25).Color,'LineWidth',2.5,'LineStyle','-.');
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence')
 %title('Avg mRNA production per nucleus') 
 ylim([0 950000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(7).AvgProd./MS2Conversion),(AvgProdAllAP(7).All95Conf./MS2Conversion),'Color',Colors(7).Color,'LineWidth',2.5,'LineStyle','-');
  errorbar(EggLength, (AvgProdAllAP(25).AvgProd./MS2Conversion),(AvgProdAllAP(25).All95Conf./MS2Conversion),'Color',Colors(25).Color,'LineWidth',2.5,'LineStyle','-.');
ylabel('transcripts produced','Color','k');
set(gca,'ycolor','k');
ylim([0 (950000/MS2Conversion)]);
%legend('Distal', 'Distal 32C','Location','best');
print( [FigDirect filesep '2xProxC17CompTotalmRNA'],'-dsvg');

%% TFBS Deletion comparisons
figure
 errorbar(EggLength,AvgProdAllAP(6).AvgProd,AvgProdAllAP(6).All95Conf,'Color',Colors(6).Color,'LineWidth',2.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(23).AvgProd, AvgProdAllAP(23).All95Conf,'Color',Colors(23).Color,'LineWidth',2.5,'LineStyle',':');
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence')
 %title('Avg mRNA production per nucleus') 
 ylim([0 1500000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(6).AvgProd./MS2Conversion),(AvgProdAllAP(6).All95Conf./MS2Conversion),'Color',Colors(6).Color,'LineWidth',2.5,'LineStyle','-');
  errorbar(EggLength, (AvgProdAllAP(23).AvgProd./MS2Conversion),(AvgProdAllAP(23).All95Conf./MS2Conversion),'Color',Colors(23).Color,'LineWidth',2.5,'LineStyle',':');
ylabel('transcripts produced','Color','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
ylim([0 (1500000/MS2Conversion)]);
saveas(gcf, [FigDirect filesep 'DoubDist_DelBcdTotProdComp','.pdf'],'pdf');

figure
 errorbar(EggLength,AvgProdAllAP(7).AvgProd,AvgProdAllAP(7).All95Conf,'Color',Colors(7).Color,'LineWidth',2.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(24).AvgProd, AvgProdAllAP(24).All95Conf,'Color',Colors(24).Color,'LineWidth',2.5,'LineStyle',':');
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence')
 %title('Avg mRNA production per nucleus') 
 ylim([0 800000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(7).AvgProd./MS2Conversion),(AvgProdAllAP(7).All95Conf./MS2Conversion),'Color',Colors(7).Color,'LineWidth',2.5,'LineStyle','-');
  errorbar(EggLength, (AvgProdAllAP(24).AvgProd./MS2Conversion),(AvgProdAllAP(24).All95Conf./MS2Conversion),'Color',Colors(24).Color,'LineWidth',2.5,'LineStyle',':');
ylabel('transcripts produced','Color','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
ylim([0 (800000/MS2Conversion)]);
saveas(gcf, [FigDirect filesep 'DoubProx_DelHbTotalmRNAComp','.pdf'],'pdf');

figure
 errorbar(EggLength,AvgProdAllAP(9).AvgProd,AvgProdAllAP(9).All95Conf,'Color',Colors(9).Color,'LineWidth',2.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(28).AvgProd, AvgProdAllAP(28).All95Conf,'Color',Colors(28).Color,'LineWidth',2.5,'LineStyle',':');
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence')
 %title('Avg mRNA production per nucleus') 
 ylim([0 1500000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(9).AvgProd./MS2Conversion),(AvgProdAllAP(9).All95Conf./MS2Conversion),'Color',Colors(9).Color,'LineWidth',2.5,'LineStyle','-');
  errorbar(EggLength, (AvgProdAllAP(28).AvgProd./MS2Conversion),(AvgProdAllAP(28).All95Conf./MS2Conversion),'Color',Colors(28).Color,'LineWidth',2.5,'LineStyle',':');
ylabel('transcripts produced','Color','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
ylim([0 (1500000/MS2Conversion)]);
saveas(gcf, [FigDirect filesep 'SE_DelBcdTotmRNAComp','.pdf'],'pdf');

figure
 errorbar(EggLength,AvgProdAllAP(9).AvgProd,AvgProdAllAP(9).All95Conf,'Color',Colors(9).Color,'LineWidth',2.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(31).AvgProd, AvgProdAllAP(31).All95Conf,'Color',Colors(31).Color,'LineWidth',2.5,'LineStyle',':');
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence')
 title('SEdHb') 
 ylim([0 1500000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(9).AvgProd./MS2Conversion),(AvgProdAllAP(9).All95Conf./MS2Conversion),'Color',Colors(9).Color,'LineWidth',2.5,'LineStyle','-');
  errorbar(EggLength, (AvgProdAllAP(31).AvgProd./MS2Conversion),(AvgProdAllAP(31).All95Conf./MS2Conversion),'Color',Colors(31).Color,'LineWidth',2.5,'LineStyle',':');
ylabel('transcripts produced','Color','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
ylim([0 (1500000/MS2Conversion)]);
saveas(gcf, [FigDirect filesep 'SE_DelHbTotmRNAComp','.pdf'],'pdf');


%% Compare SE vs inverted orientation 
figure 
figure
 errorbar(EggLength,AvgProdAllAP(9).AvgProd,AvgProdAllAP(9).All95Conf,'Color',Colors(9).Color,'LineWidth',2.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(29).AvgProd, AvgProdAllAP(29).All95Conf,'Color',Colors(29).Color,'LineWidth',2.5,'LineStyle','-.');
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence')
 %title('Avg mRNA production per nucleus') 
 ylim([0 1500000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(9).AvgProd./MS2Conversion),(AvgProdAllAP(9).All95Conf./MS2Conversion),'Color',Colors(9).Color,'LineWidth',2.5,'LineStyle','-');
  errorbar(EggLength, (AvgProdAllAP(29).AvgProd./MS2Conversion),(AvgProdAllAP(29).All95Conf./MS2Conversion),'Color',Colors(29).Color,'LineWidth',2.5,'LineStyle','-.');
ylabel('transcripts produced','Color','k');
set(gca,'ycolor','k');
ylim([0 (1500000/MS2Conversion)]);
%legend('Distal', 'Distal 32C','Location','best');

 %% Compare embryos to ID outlyers 
 cc=input('Which construct to evaluate?');
 EmbEval=[]; EmbGroup=[];
 for ee=1:length(AvgProdAllAP(cc).EmbryosProd)
     EmbEval=[[EmbEval];[AvgProdAllAP(cc).EmbryosProd(ee).AllProds(:,ExpressionBins{1,cc})]]
     for bb=1:length([AvgProdAllAP(cc).EmbryosProd(ee).AllProds(:,ExpressionBins{1,cc})])
     EmbGroup=[[EmbGroup];ee];
     end
 end
 [p,tbl,stats]=anova1(EmbEval,EmbGroup);
 [c,m]=multcompare(stats);
KruskalTest=input('Do K-W test of medians?','s');
if KruskalTest=='y'
[p,tbl,stats]=kruskalwallis(EmbEval,EmbGroup);
[c,m]=multcompare(stats,'CType','bonferroni');
end

 
%% Save the mRNA production info in the Constructs folder
DistTotalProd=[AvgProdAllAP(1).AvgProd]; ProxTotalProd=[AvgProdAllAP(2).AvgProd];
DoubDistTotalProd=[AvgProdAllAP(6).AvgProd]; DoubProxTotalProd=[AvgProdAllAP(7).AvgProd];
DoubDistOldTotalProd=[AvgProdAllAP(30).AvgProd];
EndogDistTotalProd=[AvgProdAllAP(32).AvgProd]; EndogDistmRNACounts=[AvgProdAllAP(32).AvgmRNACounts];
BothTotalProd=[AvgProdAllAP(9).AvgProd]; BothmRNACounts=[AvgProdAllAP(9).AvgmRNACounts];
DistmRNACounts=[AvgProdAllAP(1).AvgmRNACounts]; ProxmRNACounts=[AvgProdAllAP(2).AvgmRNACounts];
DoubDistmRNACounts=[AvgProdAllAP(6).AvgmRNACounts]; DoubProxmRNACounts=[AvgProdAllAP(7).AvgmRNACounts];
Dist32CTotalProd=[AvgProdAllAP(11).AvgProd]; Prox32CTotalProd=[AvgProdAllAP(12).AvgProd];
Dist32CmRNACounts=[AvgProdAllAP(11).AvgmRNACounts]; Prox32CmRNACounts=[AvgProdAllAP(12).AvgmRNACounts];
Dist17CTotalProd=[AvgProdAllAP(18).AvgProd]; Prox17CTotalProd=[AvgProdAllAP(22).AvgProd];
Dist17CmRNACounts=[AvgProdAllAP(18).AvgmRNACounts]; Prox17CmRNACounts=[AvgProdAllAP(22).AvgmRNACounts];
save([DropboxFolder filesep 'Constructs' filesep 'AllTotalmRNAProd'],'AvgProdAllAP');
save([DropboxFolder filesep 'Constructs' filesep 'EndogDistTotalProd'], 'EndogDistTotalProd');
save([DropboxFolder filesep 'Constructs' filesep 'EndogDistmRNACounts'], 'EndogDistmRNACounts');
save([DropboxFolder filesep 'Constructs' filesep 'DistmRNACounts'],'DistmRNACounts');
save([DropboxFolder filesep 'Constructs' filesep 'ProxmRNACounts'],'ProxmRNACounts');
save([DropboxFolder filesep 'Constructs' filesep 'BothmRNACounts'],'BothmRNACounts');
save([DropboxFolder filesep 'Constructs' filesep 'BothTotalProd'],'BothTotalProd');
save([DropboxFolder filesep 'Constructs' filesep 'Duplicated enhancers' filesep 'DoubDistTotalProd'],'DoubDistTotalProd');
save([DropboxFolder filesep 'Constructs' filesep 'Duplicated enhancers' filesep 'DoubDistOldTotalProd'],'DoubDistOldTotalProd');
save([DropboxFolder filesep 'Constructs' filesep 'Duplicated enhancers' filesep 'DoubProxTotalProd'],'DoubProxTotalProd');
save([DropboxFolder filesep 'Constructs' filesep 'Duplicated enhancers' filesep 'DoubDistmRNACounts'],'DoubDistmRNACounts');
save([DropboxFolder filesep 'Constructs' filesep 'Duplicated enhancers' filesep 'DoubProxmRNACounts'],'DoubProxmRNACounts');
save([DropboxFolder filesep 'Constructs' filesep 'Altered temperature data' filesep 'Dist32CTotalProd'],'Dist32CTotalProd');
save([DropboxFolder filesep 'Constructs' filesep 'Altered temperature data' filesep 'Dist32CmRNACounts'],'Dist32CmRNACounts');
save([DropboxFolder filesep 'Constructs' filesep 'Altered temperature data' filesep 'Prox32CTotalProd'],'Prox32CTotalProd');
save([DropboxFolder filesep 'Constructs' filesep 'Altered temperature data' filesep 'Prox32CmRNACounts'],'Prox32CmRNACounts');
save([DropboxFolder filesep 'Constructs' filesep 'Altered temperature data' filesep 'Dist17CTotalProd'],'Dist17CTotalProd');
save([DropboxFolder filesep 'Constructs' filesep 'Altered temperature data' filesep 'Dist17CmRNACounts'],'Dist17CmRNACounts');
save([DropboxFolder filesep 'Constructs' filesep 'Altered temperature data' filesep 'Prox17CTotalProd'],'Prox17CTotalProd');
save([DropboxFolder filesep 'Constructs' filesep 'Altered temperature data' filesep 'Prox17CmRNACounts'],'Prox17CmRNACounts');

%% 2way anova of Duration vs AP position vs genotype 
%load constructs
ConstructList= {'KrDist';'KrProx';'KrBothSep';'KrDistEmpty';'KrDistDuplicN';'KrProxDuplic';'KrBoth'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location
%Ask if only want nc14 info
%ncUse=input('Want to only use nc14?','s');

%Count for each construct

TotalProdManovaVect=[];
APManovaVect=[];
ConManovaVect=[];
for cc=1:length(ConstructList)
    Data= LoadMS2SetsCS(ConstructList{cc});
    Datalength(cc)=length(Data);
    NEmbryos = length(Data);
    APbinID=[Data(1).APbinID];
    Label = ConstructList(cc);
   
    for ee=1:NEmbryos
        PrefixName=Data(ee).Prefix;
        if ncUse=='y'
        filename=[DropboxFolder filesep PrefixName filesep 'BurstPropertiesnc14.mat'];
        else
            filename=[DropboxFolder filesep PrefixName filesep 'BurstProperties.mat'];
        end 
        load(filename);
        NumberBursts(ee,cc)=length([BurstProperties.Duration]);
        for nn=1:length(BurstProperties)
            if isempty(BurstProperties(nn).TotalmRNA)
                TotalProdManovaVect=[TotalProdManovaVect, nan];
            else
            TotalProdManovaVect=[TotalProdManovaVect, BurstProperties(nn).TotalmRNA];
            end
            APManovaVect=[APManovaVect, BurstProperties(nn).APBin];
            ConManovaVect=[ConManovaVect, cc];
        end
    end
end

[p,tbl,stats]=anovan(TotalProdManovaVect,{APManovaVect, ConManovaVect},'sstype',1,'varnames',{'AP bin','Construct'})
%multcompare(stats,'Dimension',2)

%% Run small test sample 

% Load your constructs - listed as they are in DataStatus.xlsx
ConstructList= {'KrSEChrom3','KrSEChr3_Empty','Kr3_Chr3_1xBcd'}%{'KrDist','KrProx','KrBothSep','KrDistEmpty','KrProxEmpty','KrDistDuplicN','KrProxDuplic','Kr2xProxEmpty','KrBoth','KrBothEmpty','KrDist32C','KrProx32C','KrBothSep32C','KrBoth32C','Kr2xProx32C','HbEmpty','KrDistDuplicN','KrDist17C','Kr2xDist32C','KrBoth17C','Kr2xDist17C','KrProx17C','Kr2xDistdBcd','Kr2xProxdHb','Kr2xProx17C','Kr2xDistEmpty','Kr2xDistEmpty32C','KrSEdBcd','KrInvSE','Kr2xDistLessMS2','KrSEdHb','KrEndogDist','KrEndogDist32C','KrEndogDist17C','KrSEChrom3','Kr2xDistChrom3','Kr2xProxChrom3','Kr2xDistChrom3_Empty','Kr4_1xBcd','KrSEChr3_Empty','Kr3_1xBcd','Kr4_1xHb','Kr5_1xHb','Kr3_1xHb','Kr4_1xZld', 'Kr5_1xStat92E','Kr3_1xZld','Kr3_1xStat92E','Kr4_Chr2_Chr3','Kr4_6xBcd','Kr3_Chr3_1xBcd'} %{'KrDist','KrProx','KrBothSep', 'KrDistDuplicN', 'KrProxDuplic', 'KrBoth'};% %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location

%Ask if only want nc14 info - respond 'y' 
ncUse=input('Want to only use nc14? y/n','s');
SlopeUse=input('Want to use Slope calculations? y/n','s');

%Organize all integrated fluorescence recordings for each construct
AvgmRNAProd=[];
for cc=1:length(ConstructList)
    Data= LoadMS2SetsCS(ConstructList{cc});
    Datalength(cc)=length(Data);
    NEmbryos = length(Data);
    APbinID=[Data(1).APbinID];
    Label = ConstructList(cc);
    AvgmRNAProdCon=[];
    firsttime=1;
    Timez=[];
    ConmRNAProdAllAP=[];
    ConmRNAProdSE=[];
    ConProdSD=[];
    EmbsArray=[];
     mRNAProdAllAP=[];
     mRNAProdErrorAllAP=[];
% Load the data for each embryo
    for ee=1:NEmbryos
        PrefixName=Data(ee).Prefix;
        if SlopeUse=='y'
            filename=[DropboxFolder filesep PrefixName filesep 'BurstPropertiesSlope.mat'];
        elseif ncUse=='y'
        filename=[DropboxFolder filesep PrefixName filesep 'BurstPropertiesnc14.mat'];
        else
            filename=[DropboxFolder filesep PrefixName filesep 'BurstProperties.mat'];
        end 
        load(filename);
        CompPars=[DropboxFolder filesep PrefixName filesep 'CompiledParticles.mat'];
        load(CompPars);
        
        NumberBursts(ee,cc)=length([BurstProperties.Duration]);
        
        %seperate out by AP bin
        for aa=1:length(APbinID)
            ProdAP=[];
            ProdErrorAP=[];
            ProductionAP=find([BurstProperties.APBin]==APbinID(aa)); %find the nuclei with bursts 
            if isempty(ProductionAP)
                ProdAP=[ProdAP; nan];
            else
            for bb=1:length(ProductionAP)
                ProdAP=[ProdAP;BurstProperties(ProductionAP(bb)).TotalmRNA];  %put all mRNA outputs at a given AP value in a column going down
                ProdErrorAP=[ProdErrorAP;BurstProperties(ProductionAP(bb)).TotalmRNAError];           
            end
            end
            %make array with AP bins across as columns each page an embryo
            for bb=1:length(ProdAP)
                mRNAProdAllAP(bb,aa,ee)=ProdAP(bb);
            end
            for bb=1:length(ProdErrorAP)
                mRNAProdErrorAllAP(bb,aa,ee)=ProdErrorAP(bb);
            end
            mRNAProdAllAP(mRNAProdAllAP==0)=nan;  %Remove 0's recorded in BurstProperties to indicate nucleus exist without that spot of transcription active
            mRNAProdErrorAllAP(mRNAProdErrorAllAP==0)=nan;
            
            ProductionSD(ee,aa,cc)=nanstd(mRNAProdAllAP(:,aa,ee));
        ProductionSE(ee,aa,cc)=ProductionSD(ee,aa,cc)/sqrt(sum(~isnan(ProdAP))); %n=# of nuclei that produce mRNA
            clear ProductionAP
        end
        %Compile all raw duration data for a construct in one long column
        %to put in structure
        AvgProdAllAP(cc).nc14Time(ee)=ElapsedTime(end)-ElapsedTime(nc14);
        
         AvgProdAllAP(cc).EmbryosProd(ee).MeanProd=nanmean(mRNAProdAllAP(:,:,ee));
        AvgProdAllAP(cc).EmbryosProd(ee).SE=ProductionSE(ee,:,cc);
         AvgProdAllAP(cc).EmbryosProd(ee).SD=ProductionSD(ee,:,cc);
         AvgProdAllAP(cc).EmbryosProd(ee).AllProds=[mRNAProdAllAP(:,:,ee)];
         Timez=[Timez,(ElapsedTime(end)-ElapsedTime(nc14))];
    end
    %Combine the data from all embryos of a construct
    for bb=1:size(mRNAProdAllAP,3)
            ConmRNAProdAllAP=[ConmRNAProdAllAP; mRNAProdAllAP(:,:,bb)];
        end
        for bb=1:size(ProductionSD,3)
            ConProdSD=[ConProdSD; ProductionSD(:,:,bb)];
        end
        for bb=1:size(ProductionSE,3)
            ConmRNAProdSE=[ConmRNAProdSE;ProductionSE(:,:,bb)];
        end

        AvgProdAllAP(cc).nc14Time=mean(Timez);
        AvgProdAllAP(cc).AvgProd=nanmean(ConmRNAProdAllAP,1);  %Avg mRNA production of all embryos of a construct by AP position
        AvgProdAllAP(cc).ConSD=nanmean(ConProdSD,1);
        AvgProdAllAP(cc).ConSE=nanmean(ConmRNAProdSE,1);
        AvgProdAllAP(cc).AllProds=[ConmRNAProdAllAP];
        AvgProdAllAP(cc).AllSD=nanstd([AvgProdAllAP(cc).AllProds]);
        for aa=1:length(APbinID)
            %Perform SE calc by AP bin as there are different # of data
            %points in different AP bins
            AvgProdAllAP(cc).AllSE(aa)=(([AvgProdAllAP(cc).AllSD(aa)])./sqrt(sum(~isnan(AvgProdAllAP(cc).AllProds(:,aa)))));
        end
        AvgProdAllAP(cc).All95Conf=(AvgProdAllAP(cc).AllSE).*1.95;
        AvgProdAllAP(cc).EmbryoMeans=[EmbsArray];
        clear mRNAProdAllAP mRNAProdErrorAllAP;
end

% Get rid of single values for an AP bin
for cc=1:length(ConstructList)
    for aa=1:length(APbinID)
        if sum(~isnan(AvgProdAllAP(cc).AllProds(:,aa)))==1
            AvgProdAllAP(cc).AllSD(aa)=nan;
            AvgProdAllAP(cc).AvgProd(aa)=nan;
            AvgProdAllAP(cc).AllSE(aa)=nan;
            AvgProdAllAP(cc).All95Conf(aa)=nan;
        end
    end
end

%% Plot small test set
x = 1; y = 2; z = 3;
figure
 errorbar(EggLength,AvgProdAllAP(x).AvgProd,AvgProdAllAP(x).All95Conf,'Color',Colors(35).Color,'LineWidth',3.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(y).AvgProd, AvgProdAllAP(y).All95Conf,'Color',Colors(40).Color,'LineWidth',3.5,'LineStyle',':');
 errorbar(EggLength, AvgProdAllAP(z).AvgProd, AvgProdAllAP(z).All95Conf, 'Color', 'k', 'LineWidth',3.5,'LineStyle',"--");
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 %title('Avg mRNA production per nucleus') 
 ylim([0 1400000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(x).AvgProd./MS2Conversion),(AvgProdAllAP(x).All95Conf./MS2Conversion),'Color',Colors(35).Color,'LineWidth',3.5,'LineStyle','-');
errorbar(EggLength, (AvgProdAllAP(y).AvgProd./MS2Conversion),(AvgProdAllAP(y).All95Conf./MS2Conversion),'Color',Colors(40).Color,'LineWidth',3.5,'LineStyle',':');
errorbar(EggLength, (AvgProdAllAP(z).AvgProd./MS2Conversion),(AvgProdAllAP(z).All95Conf./MS2Conversion),'Color','k','LineWidth',3.5,'LineStyle','--');
ylabel('total transcripts produced');
ylim([0 (1400000/MS2Conversion)]);
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);

%% Fx for calculating error of exp boundaries
function [AntBound, PostBound] = CalcExpBounds(x,c)
% calculate the AP bins of c-th percentile peak expression from x data 
% edited to calc AP bins >= c, which is 50% of peak homo exp -- 20210401 RW
    
    ExpLim = c;
    MeanExp = nanmean(x); %calc the avg exp in each AP bin (columns)
    [o,PeakExp] = nanmax(MeanExp); %find AP bin of peak expression
    
    Ant_Exp = MeanExp(1:PeakExp-1);
    Post_Exp = MeanExp(PeakExp+1:end);
    
    %gives the AP bin of left most that is at least % expression
    %AntBound = min(find(Ant_Exp >= prctile(MeanExp, ExpLim)))
    %Use above threshold as cutoff
    AntBound = min(find(Ant_Exp >= c));
    %Deal w if there's nothing to the left that meets the % cutoff
     if isempty(AntBound)
         AntBound = PeakExp;
     end
    % AP bin right most that is at least % exp
    %PostBound = max(find(Post_Exp >= prctile(MeanExp, ExpLim)))+ PeakExp
    %Use above threshold as cutoff
    PostBound = max(find(Post_Exp >= c))+ PeakExp
     if isempty(PostBound)
         PostBound = PeakExp;
     end
    
    
end

    




