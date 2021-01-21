% schematic figure for ERL review
% Kristina Anderson Teixeira
% January 2021

clf; clear; close all;
set(0,'DefaultLegendAutoUpdate','off');


%plot settings

    %subplot positions
    plot_width=0.32;
    plot_height=0.15;
    plot_space_vertical=(1-plot_height*5)/5-.01;
    plot_space_horizontal=(1-plot_width*2.5)/2-.05;

    %color
    facecolor_cat=[.1 .4 .9; 0.8 0 0;  .9 .527 .527 ; 1 .9 .9 ]; %blue/ red

    
    %font
    fs= 10; %font size



%% ~~~~~~~CREATE DATA FOR PLOTTING~~~~~~~~~~
age = linspace (1,100,100);
biome_names={'tropical', 'temperate', 'boreal'};

%FLUXES
% 1. by age
%for stacked plot
NEP=2.4*log10(age)-.04*(age-1); %placeholder (peaks around 2 at 23 yrs age)
GPP= NEP + (1.99+1.44)* log10(age);
ANPP_foliage = 1.66* log10(age) +.18 - .01*(age-1);
ANPP_woody = 2.31* log10(age) -.76 - .02*(age-1);
ANPP_woody_turnover=1.7*log10(age)-1; %max(0,ANPP_woody-NEP);
BNPP=0.97*log10(age)+.67 - .01*(age-1);
NPP=ANPP_foliage+ANPP_woody+BNPP;
CUE=0.679-0.153*log10(age); %from DeLucia et al. 2007 (note that equation given mistakenly leaves of the log10(age)). Collati et al. 2020 gives an update of this.
R_auto= max(0,NPP.*(1./CUE-1)); %Rauto=NPP*(1/CUE-1), and CUE equation is above
R_het=ANPP_foliage+ANPP_woody_turnover+BNPP;


%categorize by "in" (GPP components, including NEP) and "out" (Reco components)
in_fluxes = [R_auto; BNPP; ANPP_foliage;  ANPP_woody_turnover; NEP]; %matrix with all fluxes for stacked plot
in_flux_names = {'R_{auto}', 'BNPP', 'ANPP_{foliage}', 'ANPP_{woody.turnover}',  'NEP'};
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
biomass=54* log10(age); 
DW= 2.6*log10(age);

stocks = [biomass; DW]; %matrix with all stocks for stacked plot
stock_names = {'B_{tot}', 'DW'};




% 2. mature forests
Stock_ratio_temp_trop=0.9;
Stock_ratio_boreal_trop=0.7;

mature_stocks = [stocks(:,end)'; Stock_ratio_temp_trop*stocks(:,end)' ; Stock_ratio_boreal_trop*stocks(:,end)'];


%% ~~~~~~~PLOTTING~~~~~~~~~~

%set facecolors:
facecolor_in_fluxes= [0.2 0 .7;...  %' 'R_{auto}'
                    0 .7 1;...  %'BNPP'
                    0 .527 .2;... %'ANPP_{foliage}',
                    150/255 1 0 ;...%'ANPP_{woody.turnover}*'
                    255/256 255/256 0]; % NEP
facecolor_stocks= [%0.1 0 .4;...  %'B_{root-coarse}*', ' 
                    %0 .7 1;... % 'B_{root-fine}',
                    %0 .527 .27 ;... % B_{foliage}',
                    170/255 1 0;... % 'B_{ag-wood}*', 
                    %1 .8 0;...  % 'DW_{standing}*', 
                    1 99/255 99/255 ... % 'DW_{down}',
                    %; 0.7 0 .7 ;... %'OL'
                    ];

figure (1)

subplot (2,3,1)
%stack as components of GPP and Reco, outline those.
plotBarStackGroups(mature_fluxes, biome_names, facecolor_in_fluxes); 
t = title('biome differences');
ylabel ('C fluxes (Mg C ha^{-1} yr^{-1})');
set(gca, 'XTick', []); %ticks off
set(gca, 'YTick', []); %ticks off

subplot (2,3,2:3)
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
    
plot(age, in_sum, '-b', 'LineWidth', 4); hold on;
plot(age, out_sum, '--r', 'LineWidth', 4);

t = title('age trends'); 
set(gca, 'XTick', []); %ticks off
set(gca, 'YTick', []); %ticks off
ylim([0, max(in_sum)+1]);
legend ([{'R_{het}'}, in_flux_names, {'GPP' 'Reco'}], 'Location', 'BestOutside');

subplot (2,3,4)
h=bar(mature_stocks, 'stacked');
    for n=1:size(mature_stocks,2)
        h(n).FaceColor= facecolor_stocks(n,:);
    end
ylabel ('C  stocks (Mg C ha^{-1})')
set(gca, 'XTickLabel', {'Tropical' 'Temperate' 'Boreal'})
xtickangle(30)
set(gca, 'YTick', []); %ticks off

subplot (2,3,5:6) %stocks by age:

h=area (age, stocks'); %, 'LineStyle','-'); 
hold on;
plot(age, biomass, '-k', 'LineWidth', 4);
xlabel ('stand age');
set(gca, 'XTick', []); %ticks off
set(gca, 'YTick', []); %ticks off


%set facecolor:
    for n=1:size(stocks,1)
        h(n).FaceColor= facecolor_stocks(n,:);
    end

legend ([stock_names,{'B_{tot}'}], 'Location', 'BestOutside');

%% ~~~~~~~ SAVE FIGURE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print('schematic', '-dpng', '-r600')

    
