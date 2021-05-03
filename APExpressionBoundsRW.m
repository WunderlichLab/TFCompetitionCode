%% Compare AP expression boundaries of reporter constructs 
% load data
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

