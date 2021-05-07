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
% Record the AP bin of max expression for each construct 
ExpressionBins=cell(2,(length(ConstructList)));
for cc=1:length(ConstructList)
    AvgProdAP=[nanmean(AvgProdAllAP(cc).AllProds)];
    APBintoUse=find((nanmean(AvgProdAllAP(cc).AllProds)==max(nanmean(AvgProdAllAP(cc).AllProds))));
    ExpressionBins{1,cc}=APBintoUse;
    ExpressionBins{2,cc}=ConstructList{cc};
end
save([DropboxFolder filesep 'Constructs' filesep 'ConstructExpressionBins'],'ExpressionBins');

%% Homozygous expression levels as fx of hemizygous levels - Figure 1
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

%% KrDistal in presence of enhancer only (no prom/MS2)
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
% 2 sample t-test
[h,p] = ttest2([AvgProdAllAP(63).AllProds(:,ExpressionBins{1,63})],[AvgProdAllAP(4).AllProds(:,ExpressionBins{1,4})]);
% 2 sample KS test
[h,p] = kstest2([AvgProdAllAP(63).AllProds(:,ExpressionBins{1,63})],[AvgProdAllAP(4).AllProds(:,ExpressionBins{1,4})]);

%% Reporter expression in presence of TFBS competitor arrays
% Graph of shadow pair with TFBS arrays

figure
SE_TFs = [41,47,44,48]; %constructs w the TFBS arrays on homologous chromosome
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

%% Dosage dependent effect of Bcd binding site arrays 
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

%% Compare reporter competition between insertion sites
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

% duplicated distal
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

ylabel('homozygous average total mRNA per allele');
xlabel('hemizygous average total mRNA per allele');
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'DupDist_HemiHomo_Chr2_3','.pdf'],'pdf');
