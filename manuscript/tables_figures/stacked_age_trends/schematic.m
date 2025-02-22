% schematic figure for ERL review
% Kristina Anderson Teixeira
% January 2021

clf; clear; close all;
set(0,'DefaultLegendAutoUpdate','off');

%plot settings
    %subplot positions
    plot_width_biome=0.27; 
    plot_width_age_f=0.57; 
    plot_width_age_s=0.507;
    plot_height=0.37;
    left_marg=.09;
    lower_marg=.1;
    plot_space_vertical=.1;
    plot_space_horizontal=.04;
    pos_biome_f=[left_marg lower_marg+1*(plot_height+plot_space_vertical) plot_width_biome plot_height];
    pos_biome_s=[left_marg lower_marg+0*(plot_height+plot_space_vertical) plot_width_biome plot_height];
    pos_age_f=[left_marg+plot_width_biome+plot_space_horizontal lower_marg+1*(plot_height+plot_space_vertical) plot_width_age_f plot_height];
    pos_age_s=[left_marg+plot_width_biome+plot_space_horizontal lower_marg+0*(plot_height+plot_space_vertical) plot_width_age_s plot_height];

%set facecolors:
    facecolor_in_fluxes= [0.2 0 .5;...  %' 'R_{auto}'
                        0 .7 1;...  %'BNPP'
                        0 .527 .2;... %'ANPP_{foliage}',
                        150/255 1 0 ;...%'ANPP_{woody.turnover}*'
                        255/256 255/256 0]; % NEP
    facecolor_stocks= [%0.1 0 .4;...  %'B_{root-coarse}*', ' 
                        0 .7 1;... % 'B_{root-fine}',
                        0 .527 .27 ;... % B_{foliage}',
                        170/255 1 0;... % 'B_{ag-wood}*', 
                        %1 .8 0;...  % 'DW_{standing}*', 
                        1 99/255 99/255; ... % 'DW_{down}',
                        0.7 0 .7 ];... %'OL'
                     
%% ~~~~~~~CREATE DATA FOR PLOTTING~~~~~~~~~~
age = linspace (1,100,100);
biome_names={'tropical', 'temperate', 'boreal'};

%FLUXES
% 1. by age
%for stacked plot
NEP=2.4*log10(age)-.04*(age-1); 
GPP= NEP + (1.99+1.44)* log10(age);
ANPP_foliage = 1.66* log10(age) +.18 - .01*(age-1);
ANPP_woody = 2.31* log10(age) -.76 - .02*(age-1);
ANPP_woody_turnover=1.7*log10(age)-1; %max(0,ANPP_woody-NEP);
BNPP_fine=1.46*log10(age) - .01*(age-1);
NPP=ANPP_foliage+ANPP_woody+BNPP_fine; %missing BNPP_coarse, but equations for BNPP_fine and BNPP are very similar
CUE=0.679-0.153*log10(age); %from DeLucia et al. 2007 (note that equation given mistakenly leaves of the log10(age)). Collati et al. 2020 gives an update of this.
R_auto= max(0,NPP.*(1./CUE-1)); %Rauto=NPP*(1/CUE-1), and CUE equation is above
R_het=ANPP_foliage+ANPP_woody_turnover+BNPP_fine;


%categorize by "in" (GPP components, including NEP) and "out" (Reco components)
in_fluxes = [R_auto; BNPP_fine; ANPP_foliage;  ANPP_woody_turnover; NEP]; %matrix with all fluxes for stacked plot
in_flux_names = {'R_{auto}', 'BNPP_{fine}', 'ANPP_{foliage}', 'woody turnover*',  'NEP (incl. NPP_{woody})*'};
out_fluxes = [R_auto; R_het];
out_flux_names = {'R_{auto}', 'R_{het}'};

%sum in and out fluxes
in_sum = sum(in_fluxes,1);
out_sum=sum(out_fluxes,1);

% 2. mature fluxes
Flux_ratio_temp_trop=0.5;
Flux_ratio_boreal_trop=0.25;

%for grouped stacked plot
mature_fluxes=0*ones(3,2,size(in_fluxes,1)); %empty 3D matrix (Group=biome, Stack=bar(GPP vs Reco), StackElement=component)
mature_fluxes(:,1,1:size(in_fluxes,1))=[in_fluxes(:,end)'; Flux_ratio_temp_trop*in_fluxes(:,end)' ; Flux_ratio_boreal_trop*in_fluxes(:,end)']; %components of GPP
mature_fluxes(:,2,1:size(out_fluxes,1))=[out_fluxes(:,end)'; Flux_ratio_temp_trop*out_fluxes(:,end)' ; Flux_ratio_boreal_trop*out_fluxes(:,end)'];

%others

%STOCKS
% 1. by age
%for stacked plot
B_woody=54* log10(age); 
DW= 5.27*log10(age);
OL= 4.5*log10(age);
B_foliage=max(0,1.2* log10(age)+1);
B_fineroot=max(0,1.2* log10(age)+1);

stocks = [B_fineroot; B_foliage; B_woody; DW; OL]; %matrix with all stocks for stacked plot
stock_names = {'B_{root-fine}', 'B_{foliage}', 'B_{woody}', 'DW_{tot}', 'OL'};




% 2. mature forests
Stock_ratio_by_biomes=[1 0.9 0.7; 1 0.9 0.7; 1 0.9 0.7; 1 1 .8; 1 1.5 2];

mature_stocks = [stocks(:,end)'; stocks(:,end)'; stocks(:,end)'].*Stock_ratio_by_biomes';


stocks_std = sum(mature_stocks, 2)'.*[.4 1 .4];

%% ~~~~~~~PLOTTING~~~~~~~~~~



figure (1)

subplot ('Position', pos_biome_f)
%stack as components of GPP and Reco, outline those.
plotBarStackGroups(mature_fluxes, biome_names, facecolor_in_fluxes); 
t = title('biome differences');
ylabel ('C fluxes (Mg C ha^{-1} yr^{-1})');
set(gca, 'YTick', []); %ticks off
set(gca, 'XTickLabel', {'Tropical' 'Temperate' 'Boreal'})
xtickangle(30)
text(.05,.9,'A.', 'Units','Normalized' ) 
box on

subplot ('Position', pos_age_f)
%flux plot:
h2=area (90:100, out_fluxes(2, 90:100)'); hold on; %will be covered just for legend.
    h2(1).FaceColor= [.7 0 .7];

h=area (age, in_fluxes'); hold on;
    for n=1:size(in_fluxes,1)
        h(n).FaceColor= facecolor_in_fluxes(n,:);
    end

h2=area (95:100, out_fluxes(:, 95:100)' ,'HandleVisibility','off'); hold on;
    h2(1).FaceColor= facecolor_in_fluxes(1,:);
    h2(1).EdgeColor= facecolor_in_fluxes(1,:);
    h2(2).FaceColor= [.7 0 .7];
    
plot(age, in_sum, '-b', 'LineWidth', 3.4); hold on;
plot(age, out_sum, '--r', 'LineWidth', 3.4);
plot (age(1:20), 2.4*R_het(end).*exp(-.27*age(1:20))+sum(out_sum(:,1:20),1), ':r', 'LineWidth', 2)

t = title('age trends'); 
xlabel ('stand age');
set(gca, 'XTick', []); %ticks off
set(gca, 'YTick', []); %ticks off
ylim([0, max(in_sum)+1]);
legend ([{'R_{het}'}, in_flux_names, {'GPP' 'R_{eco}' 'legacy R_{eco}'}], 'Location', 'BestOutside');
legend('Boxoff')
text(.05,.9,'B.', 'Units','Normalized' ) 

subplot ('Position', pos_biome_s)
h=bar(mature_stocks, 'stacked'); hold on;
    for n=1:size(mature_stocks,2)
        h(n).FaceColor= facecolor_stocks(n,:);
    end
er=errorbar([1:3],sum(mature_stocks, 2), 50-sum(mature_stocks, 2), stocks_std);
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
er.LineWidth= 1.7;
ylabel ('C  stocks (Mg C ha^{-1})')
set(gca, 'YTick', []); %ticks off
set(gca, 'XTickLabel', {'Tropical' 'Temperate' 'Boreal'})
xtickangle(30)
text(.05,.9,'C.', 'Units','Normalized' ) 

subplot ('Position', pos_age_s) %stocks by age:

h=area (age, stocks'); hold on;
plot(age, B_woody+B_foliage+B_fineroot, '-k', 'LineWidth', 3.4);
plot (age(1:20), B_woody(end).*exp(-.27*age(1:20))+sum(stocks(:,1:20),1), ':', 'Color',facecolor_stocks(4,:), 'LineWidth', 2)
xlabel ('stand age');
set(gca, 'XTick', []); %ticks off
set(gca, 'YTick', []); %ticks off


%set facecolor:
    for n=1:size(stocks,1)
        h(n).FaceColor= facecolor_stocks(n,:);
    end

legend ([stock_names,{'B_{tot}', 'legacy DW_{tot}'}], 'Location', 'BestOutside');
legend('Boxoff')
text(.05,.9,'D.', 'Units','Normalized' ) 

%% ~~~~~~~ SAVE FIGURE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print('schematic', '-dpng', '-r600')

    
