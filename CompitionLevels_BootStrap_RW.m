%% Compare competition levels between reporters or configurations
% load data

% plotting info 
Filename = [DropboxFolder filesep 'Constructs' filesep 'AllTotalmRNAProd'];
load(Filename)

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
%% Bootstrap to get p value of competition, comparing all chr2 enhancer reporters
Ntimes = 1000; % set the number of times to try bootstrapping
counter=0;
for cc = [1,6,2,7,9]; % set constructs to do this for
    counter = counter+1;
    ConData = [AvgProdAllAP(cc).AllProds(:,ExpressionBins{1,cc})]; % sample from this
    if cc == 1
        HemiData = [AvgProdAllAP(4).AllProds(:,ExpressionBins{1,4})]; %where to sample hemi data
        dd =4;
    elseif cc == 2
        HemiData = [AvgProdAllAP(5).AllProds(:,ExpressionBins{1,5})];
        dd=5;
    elseif cc == 6
        HemiData = [AvgProdAllAP(26).AllProds(:,ExpressionBins{1,26})];
        dd=26;
    elseif cc == 7
        HemiData = [AvgProdAllAP(8).AllProds(:,ExpressionBins{1,8})];
        dd=8;
    elseif cc ==9
        HemiData = [AvgProdAllAP(10).AllProds(:,ExpressionBins{1,10})];
        dd=10;
    end
    %calc real fraction (hemi / homo)
    TrueFract(counter) = AvgProdAllAP(dd).AvgProd(ExpressionBins{1,dd}) / AvgProdAllAP(cc).AvgProd(ExpressionBins{1,cc});  %calc from orig data
    %Real test statistic
    T_True(counter) = TrueFract(counter) - 1; %Bc null hypothesis is that hemi = homo expression
    Boot_All=[];
    T_AllTest=[];
    Diff_All=[];
    %bootstrap error in exp bound shifts
    for tt = 1:Ntimes
        %get homozygous subsampled
        y = datasample(ConData,size(ConData,1),1); %Randomly sample rows, 1 long column 
        %hemizygous
        x = datasample(HemiData, size(HemiData,1),1); 
        %perform calculation, ie divide hemizygous by homozygous
        BootIndiv = nanmean(x) / nanmean(y);
        Boot_All = [Boot_All; BootIndiv];
        %emperical CI 
        Diff_Indiv = BootIndiv - TrueFract(counter);
        Diff_All = [Diff_All; Diff_Indiv];
    end
     %emperical bootstrap CI
     Diff_All = sort(Diff_All); % sort smallest to largest
     D_2(counter) = prctile(Diff_All,2.5); %bottom 2.5 percentile 
     D_98(counter) = prctile(Diff_All, 97.5); %top 2.5 percentile 
     EmpBoot(counter,1) = D_2(counter);
     EmpBoot(counter,2) = D_98(counter);
     
end

%% Make bar graph of fraction hemi expression using bootstrapped errorbars
figure
counter = 0;
for cc =[1,6,2,7,9]; %order so distals by distal, prox by prox
    counter =counter+1;
    %graph as % hemi higher than homo
    bar(counter, (TrueFract(counter)-1));
    hold on 
    errorbar(counter, (TrueFract(counter)-1),EmpBoot(counter,1), EmpBoot(counter,2),'LineWidth',2.5,'Color','k');
end
ylabel('% higher exp in hemizygotes');
ylim([0 1]);
yticks([0.2, 0.4, 0.6, 0.8, 1]);
yticklabels({'20','40','60','80','100'});
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 5]);
saveas(gcf, [FigDirect filesep 'HemiHomoPercComp_ErrorBars','.pdf'],'pdf');

%% Bootstrap for the error of % hemizygous expression comparing attP sites
 %Need to calc fraction of hemizygous expression after random sampling w
 %replacement 
 
 Ntimes = 1000; % set the number of times to try bootstrapping
counter=0;
 
 for cc = [6,36,9,35]
     counter = counter+1;
    ConData = [AvgProdAllAP(cc).AllProds(:,ExpressionBins{1,cc})]; % sample from this
    if cc == 6
        HemiData = [AvgProdAllAP(26).AllProds(:,ExpressionBins{1,26})];
        dd=26;
    elseif cc == 36
        HemiData = [AvgProdAllAP(38).AllProds(:,ExpressionBins{1,38})];
        dd=38;
    elseif cc == 9
        HemiData = [AvgProdAllAP(10).AllProds(:,ExpressionBins{1,10})];
        dd = 10;
    elseif cc == 35
        HemiData = [AvgProdAllAP(40).AllProds(:,ExpressionBins{1,40})];
        dd =40;
    end
     %calc real fraction (hemi / homo)
    TrueFract(counter) = AvgProdAllAP(dd).AvgProd(ExpressionBins{1,dd}) / AvgProdAllAP(cc).AvgProd(ExpressionBins{1,cc});  %calc from orig data
    TrueFract_PercHemi(counter) = 1/TrueFract(counter); % Calc % of hemi exp
    Boot_All=[];
    Boot_All_PH =[];
    Diff_All=[];
    Diff_All_PH=[];
    
    %bootstrap error in exp bound shifts
    for tt = 1:Ntimes
        %get homozygous subsampled
        y = datasample(ConData,size(ConData,1),1); %Randomly sample rows, 1 long column 
        %hemizygous
        x = datasample(HemiData, size(HemiData,1),1); 
        %perform calculation, ie divide hemizygous by homozygous
        BootIndiv = nanmean(x) / nanmean(y);
        Boot_All = [Boot_All; BootIndiv];
        %Also do inverse 
        BootIndiv_PH = nanmean(y) / nanmean(x);
        Boot_All_PH = [Boot_All_PH; BootIndiv_PH];
        %emperical CI 
        Diff_Indiv = BootIndiv - TrueFract(counter);
        Diff_All = [Diff_All; Diff_Indiv];
        %For perc hemi as well
        Diff_Indiv_PH = BootIndiv_PH - TrueFract_PercHemi(counter);
        Diff_All_PH = [Diff_All_PH; Diff_Indiv_PH];
    end
     %emperical bootstrap CI
     Diff_All = sort(Diff_All); % sort smallest to largest
     D_2(counter) = prctile(Diff_All,2.5); %bottom 2.5 percentile 
     D_98(counter) = prctile(Diff_All, 97.5); %top 2.5 percentile 
     EmpBoot(counter,1) = D_2(counter);
     EmpBoot(counter,2) = D_98(counter);
     % Also for % hemi calc
      %emperical bootstrap CI
     Diff_All_PH = sort(Diff_All_PH); % sort smallest to largest
     D_2_PH(counter) = prctile(Diff_All_PH,2.5); %bottom 2.5 percentile 
     D_98_PH(counter) = prctile(Diff_All_PH, 97.5); %top 2.5 percentile 
     EmpBoot_PH(counter,1) = D_2_PH(counter);
     EmpBoot_PH(counter,2) = D_98_PH(counter);
 end
    
%% Figure to compare competition on chr2 vs Ch3

% % higher expression in hemizygotes vs homozygotes
figure
counter = 0;
for cc =[6,36,9,35];
    counter =counter+1;
    %graph as % hemi higher than homo
    bar(counter, (TrueFract(counter)-1));
    hold on 
    errorbar(counter, (TrueFract(counter)-1),EmpBoot(counter,1), EmpBoot(counter,2),'LineWidth',2.5,'Color','k');
end
ylabel('% higher exp in hemizygotes');
ylim([0 1]);
yticks([0.2, 0.4, 0.6, 0.8, 1]);
yticklabels({'20','40','60','80','100'});
xticks([1 2 3 4]);
xticklabels({'DDCh2','DDCh3', 'SECh2','SECh3'});
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 5]);
saveas(gcf, [FigDirect filesep 'HemiHomoCompCompare_Chr2v3','.pdf'],'pdf');

% Perc hemi expression 
figure
counter = 0;
for cc =[6,36,9,35];
    counter =counter+1;
    %graph as % hemi higher than homo
    bar(counter, TrueFract_PercHemi(counter));
    hold on 
    errorbar(counter, TrueFract_PercHemi(counter),EmpBoot_PH(counter,1), EmpBoot_PH(counter,2),'LineWidth',2.5,'Color','k');
end
ylabel('% hemizygous expression');
ylim([0 1]);
yticks([0.2, 0.4, 0.6, 0.8, 1]);
yticklabels({'20','40','60','80','100'});
xticks([1 2 3 4]);
xticklabels({'DDCh2','DDCh3', 'SECh2','SECh3'});
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 5]);
saveas(gcf, [FigDirect filesep 'Ch2v3PercHemiComp','.pdf'],'pdf');

%Here make the inset figures for each construct 
% duplicated distal
figure
counter = 0; 
for cc = [6, 36];
        counter =counter+1;
    %graph as % hemi higher than homo
    bar(counter, TrueFract_PercHemi(counter),'FaceColor',Colors(6).Color);
    hold on 
    errorbar(counter, TrueFract_PercHemi(counter),EmpBoot_PH(counter,1), EmpBoot_PH(counter,2),'LineWidth',2.5,'Color','k');
end
ylabel('% hemizygous expression');
ylim([0 1]);
yticks([0.2, 0.4, 0.6, 0.8, 1]);
yticklabels({'20','40','60','80','100'});
xticks([1,2]);
xticklabels({'Chr2', 'Chr3'});
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'DupDist_RelHomoLvls_Chr2v3','.svg'],'svg');

% inset figure for shadow pair
figure
counter = 2; 
for cc = [9, 35];
        counter =counter+1;
    %graph as % hemi higher than homo
    bar(counter-2, TrueFract_PercHemi(counter),'FaceColor',Colors(9).Color);
    hold on 
    errorbar(counter-2, TrueFract_PercHemi(counter),EmpBoot_PH(counter,1), EmpBoot_PH(counter,2),'LineWidth',2.5,'Color','k');
end
ylabel('% hemizygous expression');
ylim([0 1]);
yticks([0.2, 0.4, 0.6, 0.8, 1]);
yticklabels({'20','40','60','80','100'});
xticks([1,2]);
xticklabels({'Chr2', 'Chr3'});
set(gca, 'YColor','k');
set(gca, 'Box','On','FontSize', fontsize, 'FontName', fontname,'LineWidth',0.5,'YColor','k');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[xLeft yTop 7 6]);
saveas(gcf, [FigDirect filesep 'SE_RelHomoLvls_Chr2v3','.svg'],'svg');
