%% Compare AP expression boundaries of reporter constructs 
% Load your constructs - listed as they are in DataStatus.xlsx
ConstructList= {'KrDist','KrProx','KrBothSep','KrDistEmpty','KrProxEmpty','KrDistDuplicN','KrProxDuplic','Kr2xProxEmpty','KrBoth','KrBothEmpty','KrDist32C','KrProx32C','KrBothSep32C','KrBoth32C','Kr2xProx32C','HbEmpty','KrDistDuplicN','KrDist17C','Kr2xDist32C','KrBoth17C','Kr2xDist17C','KrProx17C','Kr2xDistdBcd','Kr2xProxdHb','Kr2xProx17C','Kr2xDistEmpty','Kr2xDistEmpty32C','KrSEdBcd','KrInvSE','Kr2xDistLessMS2','KrSEdHb','KrEndogDist','KrEndogDist32C','KrEndogDist17C','KrSEChrom3','Kr2xDistChrom3','Kr2xProxChrom3','Kr2xDistChrom3_Empty','Kr4_1xBcd','KrSEChr3_Empty','Kr3_1xBcd','Kr4_1xHb','Kr5_1xHb','Kr3_1xHb','Kr4_1xZld', 'Kr5_1xStat92E','Kr3_1xZld','Kr3_1xStat92E','Kr4_Chr2_Chr3','Kr4_6xBcd','Kr3_Chr3_1xBcd','Kr5_1xBcd','KrBoth_BcdGFP','KrDist_BcdGFP','Kr3_Chr2_Chr3','Kr3_6xBcd','HbP2_HbP2','Kr4_1xTwi','HbP2_6xBcd','Kr1_6xBcd','Kr4_3xBcd','4F_Kr4_Kr4','Kr2_Kr2NoProm','Kr3_3xBcd','4F_Kr4_0', 'HbP2_1xBcd'} %{'KrDist','KrProx','KrBothSep', 'KrDistDuplicN', 'KrProxDuplic', 'KrBoth'};% %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, ~, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location

% load data
Filename = [DropboxFolder filesep 'Constructs' filesep 'AllTotalmRNAProd'];
load(Filename)
Filename3 = [DropboxFolder filesep 'Constructs' filesep 'ConstructExpressionBins'];
load(Filename3);
PrefixName = '2018-10-24-Kr4_Kr4' %change this to load random file for compiledparticles
Filename2 = [DropboxFolder filesep PrefixName filesep 'CompiledParticles.mat'];
load(Filename2);
% set plotting info 
EggLength=APbinID.*100;

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

ArrayColor = [179 16 182] ./255;
Array3Color =[155 16 152] ./255;
Array6Color = [115 16 107]./255;


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
%save([DropboxFolder filesep 'Constructs' filesep 'AllTotalmRNAProd'],'AvgProdAllAP');

%Set which constructs (and corresponding peak expression AP bins) are interested in
SE_Bcd = [10,9,41,64,56];
APExpBins_SE = [ExpressionBins{1,[SE_Bcd]}];
DD_Bcd = [26,6,39,61,50];
APExpBins_DD = [ExpressionBins{1,[DD_Bcd]}];
%% Bootstrapping error bars for boundary shift
%use 1/2 max of homozygous construct
CutoffLim = nanmax([AvgProdAllAP(6).AvgProd]) *0.5;%input('What percentile to use?');
Ntimes = 1000; % set the number of times to try bootstrapping
counter=0;
for cc = [26,6,39,61,50]; % set constructs to do this for
    counter = counter+1;
    ConData = [AvgProdAllAP(cc).AllProds]; % sample from this
    [TrueAntBound, TruePostBound] = CalcExpBounds(ConData,CutoffLim); %calc from orig data
    % record for later use
    AntBound_DD(counter) = TrueAntBound;
    PostBound_DD(counter) = TruePostBound;

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
     % emperical bootstrap CI's
     Diff_All_Ant = BootBound_All(:,1) - TrueAntBound; %diff btwn orig data and resampled
     Diff_All_Post = BootBound_All(:,2) - TruePostBound;
     Diff_All_Ant = sort(Diff_All_Ant); % sort smallest to largest
     D_2_Ant(counter) = prctile(Diff_All_Ant,2.5); %bottom 2.5 percentile 
     D_98_Ant(counter) = prctile(Diff_All_Ant, 97.5); %top 2.5 percentile 
     Diff_All_Post = sort(Diff_All_Post);
     D_2_Post(counter) = prctile(Diff_All_Post,2.5); %bottom 2.5 percentile 
     D_98_Post(counter) = prctile(Diff_All_Post, 97.5); %top 2.5 percentile 
     EmpBoot(counter,1) = D_2_Ant(counter); %lower ant boundary
     EmpBoot(counter, 2) = D_98_Ant(counter); %higher ant boundary
     EmpBoot(counter,3) = D_2_Post(counter);
     EmpBoot(counter,4) = D_98_Post(counter);
end

% plot the data 
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
    errorbar(EggLength(AntBound_DD(bb)), bb, EmpBoot(bb,1), EmpBoot(bb,2), 'horizontal','Color','k','LineWidth',2.5);
    errorbar(EggLength(PostBound_DD(bb)), bb, EmpBoot(bb,3), EmpBoot(bb,4), 'horizontal','Color','k','LineWidth',2.5); 
    % Set transparency (undocumented)
    %set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h.Line.ColorData(1:3); 255*alpha])
    
    %h2= errorbar(EggLength(PostBound_DD(bb)),bb,BootBoundStd(bb,2),'horizontal','Color','k','LineWidth', 2.5);
    %set([h2.Bar, h2.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h2.Line.ColorData(1:3); 255*alpha])

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
saveas(gcf, [FigDirect filesep '2xDist_AntShifts_ErrorBars_EmpCI','.pdf'],'pdf');

%%
% Shadow pair
CutoffLim = nanmax([AvgProdAllAP(9).AvgProd]) *0.5;
counter=0;
for cc = [10,9,41,64,56]; % set constructs to do this for
    counter = counter+1;
    ConData = [AvgProdAllAP(cc).AllProds]; % sample from this
    [TrueAntBound, TruePostBound] = CalcExpBounds(ConData,CutoffLim); %calc from orig data
    % record for plotting use
     AntBound_SE(counter) = TrueAntBound; 
     PostBound_SE(counter) = TruePostBound;

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
     
     % emperical bootstrap CI's
     Diff_All_Ant = BootBound_All(:,1) - TrueAntBound; %diff btwn orig data and resampled
     Diff_All_Post = BootBound_All(:,2) - TruePostBound;
     Diff_All_Ant = sort(Diff_All_Ant); % sort smallest to largest
     D_2_Ant(counter) = prctile(Diff_All_Ant,2.5); %bottom 2.5 percentile 
     D_98_Ant(counter) = prctile(Diff_All_Ant, 97.5); %top 2.5 percentile 
     Diff_All_Post = sort(Diff_All_Post);
     D_2_Post(counter) = prctile(Diff_All_Post,2.5); %bottom 2.5 percentile 
     D_98_Post(counter) = prctile(Diff_All_Post, 97.5); %top 2.5 percentile 
     EmpBoot(counter,1) = D_2_Ant(counter); %lower ant boundary
     EmpBoot(counter, 2) = D_98_Ant(counter); %higher ant boundary
     EmpBoot(counter,3) = D_2_Post(counter);
     EmpBoot(counter,4) = D_98_Post(counter);
end

%plot the data
x_sizeNarrow =3.7;
y_sizeNarrow =5.52; 
figure
Grey =[0.5 0.5 0.5];
Array3Color =[155 16 152] ./255;
% Set transparency level (0:1) for errorbars
alpha = 0.65;

figure
Grey =[0.5 0.5 0.5];
Array3Color =[155 16 152] ./255;
SE_Colors = [Colors(9).Color; Colors(9).Color; ArrayColor; Array3Color; Array6Color];
for bb = 1:length(SE_Bcd)
    plot([EggLength(AntBound_SE(bb)) EggLength(PostBound_SE(bb))],[bb bb],'LineWidth',5,'Color',SE_Colors(bb,:));
    hold on
    %errorbar(EggLength(AntBound_SE(bb)),bb,BootBoundStd(bb,1),'horizontal','Color', 'k','LineWidth',2.5);
    errorbar(EggLength(AntBound_SE(bb)), bb, EmpBoot(bb,1), EmpBoot(bb,2), 'horizontal','Color','k','LineWidth',2.5);
    errorbar(EggLength(PostBound_SE(bb)), bb, EmpBoot(bb,3), EmpBoot(bb,4), 'horizontal','Color','k','LineWidth',2.5);

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
saveas(gcf, [FigDirect filesep 'SE_AntShifts_ErrorBars_EmpCI','.pdf'],'pdf');

%% Plot the data 

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

