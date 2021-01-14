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
R_auto = 2* log10(age);
NPP = R_auto;

fluxes = [R_auto; NPP]; %matrix with all fluxes for stacked plot
flux_names = {'R_{auto}', 'NPP'};

%others
GPP = sum(fluxes,1);

% 2. mature fluxes
Flux_ratio_temp_trop=0.5;
Flux_ratio_boreaal_trop=0.25;

%for stacked plot
mature_fluxes = [fluxes(:,end)'; Flux_ratio_temp_trop*fluxes(:,end)' ; Flux_ratio_boreaal_trop*fluxes(:,end)'];

%others

%STOCKS
% 1. by age
%for stacked plot
biomass=54* log10(age); 
DW= 2.6*log10(age);

stocks = [biomass; DW]; %matrix with all stocks for stacked plot
stock_names = {'B_{tot}', 'DW'};

%others



% 2. mature forests
Stock_ratio_temp_trop=0.9;
Stock_ratio_boreaal_trop=0.7;

mature_stocks = [stocks(:,end)'; Stock_ratio_temp_trop*stocks(:,end)' ; Stock_ratio_boreaal_trop*stocks(:,end)'];


%% ~~~~~~~PLOTTING~~~~~~~~~~

figure (1)
subplot (2,3,1:2)
%flux plot:
h=area (age, fluxes'); %, 'LineStyle','-'); 
hold on;
plot(age, GPP, '-g', 'LineWidth', 3);
t = title('age trends'); 
ylabel ('fluxes');
set(gca, 'XTick', []); %ticks off
set(gca, 'YTick', []); %ticks off


subplot (2,3,3)
%bar(mature_fluxes, 'stacked');
plotBarStackGroups(mature_fluxes, biome_names);
t = title('biome differences');
legend (flux_names);
set(gca, 'XTick', []); %ticks off
set(gca, 'YTick', []); %ticks off

subplot (2,3,4:5)
%flux plot:
h=area (age, stocks'); %, 'LineStyle','-'); 
hold on;
plot(age, biomass, '-k', 'LineWidth', 3);
ylabel ('stocks');
xlabel ('stand age');
set(gca, 'XTick', []); %ticks off
set(gca, 'YTick', []); %ticks off

subplot (2,3,6)
bar(mature_stocks, 'stacked');
legend (stock_names);
xlabel ('biome');
set(gca, 'XTick', []); %ticks off
set(gca, 'YTick', []); %ticks off


%% ~~~~~~~ SAVE FIGURE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print('schematic', '-dpng', '-r600')

    
