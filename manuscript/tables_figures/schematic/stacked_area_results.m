% age trend summaries for ERL review
% Kristina Anderson Teixeira
% January 2021

clf; clear; close all;
set(0,'DefaultLegendAutoUpdate','off');

%% ~~~~~~~READ IN DATA~~~~~~~~~~

age = linspace (1,100,100);
biome_names={'tropical', 'temperate', 'boreal'};

%FLUXES
% 1. by age
%for stacked plot
NEP=4*log10(age)-.08*(age-1); %placeholder (peaks around 2 at 23 yrs age)
GPP= 11.15 + (1.99+1.44)* log10(age);
ANPP_foliage = 1.66* log10(age) +.18;
ANPP_woody = 2.31* log10(age) -.76;
ANPP_woody_turnover=max(0,ANPP_woody-NEP);
BNPP=0.97*log10(age)+.67;
NPP=ANPP_foliage+ANPP_woody+BNPP;
R_auto = NPP.*(1+.001*age); %insufficient data
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

biomes= {'Tropical broadleaf','Temperate broadleaf', 'Temperate conifer', 'Boreal conifer'};
for b =3:3
%% ~~~~~~~PLOTTING~~~~~~~~~~

figure (b)
 
subplot (2,4,1:3)
%flux plot:
h=area (age, in_fluxes'); %, 'LineStyle','-'); 
hold on;
h=area (age, -1*out_fluxes'); %, 'LineStyle','-'); 
hold on;
plot(age, in_sum, '-b', 'LineWidth', 3); hold on;
plot(age, GPP, '--b', 'LineWidth', 3); hold on;
plot(age, NEP, '-w', 'LineWidth', 3);
t = title(biomes(b));
ylabel ('C  stocks (Mg C ha^{-1})')

subplot (2,4,4) %mature C fluxes
bar(mature_stocks(b,:), 'stacked');
legend (stock_names);


subplot (2,4,5:7)
%flux plot:
h=area (age, stocks'); %, 'LineStyle','-'); 
hold on;
plot(age, biomass, '-k', 'LineWidth', 3);
xlabel ('stand age');
ylabel ('C  stocks (Mg C ha^{-1})')

subplot (2,4,8) %mature C stocks
bar(mature_stocks, 'stacked');
legend (stock_names);
xlabel ('mature stands');


end