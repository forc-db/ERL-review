% age trend summaries for ERL review
% Kristina Anderson Teixeira
% January 2021

clf; clear; close all;
set(0,'DefaultLegendAutoUpdate','off');

%% SETTINGS
figures_dir='/Users/kteixeira/Dropbox (Smithsonian)/GitHub/ForC-db/ERL-review/manuscript/tables_figures';
working_dir='/Users/kteixeira/Dropbox (Smithsonian)/GitHub/ForC-db/ERL-review/manuscript/tables_figures/schematic';

%% READ IN & PREPARE DATA
%read in data:
cd(figures_dir)
model_summaries=readtable('SI_age_trend_model_summaries.csv');
%probably need to remove specieal formatting from model_summaries.Variable

%prepare some variables
variables = unique(model_summaries.Variable);
biomes= {'Tropical broadleaf','Temperate broadleaf', 'Temperate conifer', 'Boreal conifer'};
%{'BiomeTropical broadleaf','BiomeTemperate broadleaf','BiomeTemperate conifer','BiomeBoreal conifer'}
params_matrix=cell2table(num2cell(NaN*ones(length(variables), 10)),'VariableNames',{'variable', 'beta' 'betaTrB' 'betaTeB' 'betaTeN' 'betaBoN' 'intTrB' 'intTeB' 'intTeN' 'intBoN'});
params_matrix.variable=variables;

for v=1:length(variables) %cycle through variables and fill in params_matrix
    
end


for b =3:3
%% GENERATE EQUATIONS FOR PLOTTING
age = linspace (1,100,100);


% 1 FLUXES
% 1.1 define age trends
% 1.1.1 equations (currently entered manually; to be automated)
NEP=4*log10(age)-.08*(age-1); %placeholder (peaks around 2 at 23 yrs age)
GPP= 11.15 + (1.99+1.44)* log10(age);
ANPP_foliage = 1.66* log10(age) +.18;
ANPP_woody = 2.31* log10(age) -.76;
BNPP=0.97*log10(age)+.67;
NPP=(-.27+1.44)*log10(age)+11.15;
R_eco=.91*log10(age)+11.08;
R_het_soil=0.18*log10(age)+3.88;

% 1.1.2 calculated fluxes
NEP_calc=GPP-R_eco;
NPP_calc=ANPP_foliage+ANPP_woody+BNPP;
R_auto_calc1 = GPP.*(0.5+.001*age); %insufficient data
R_auto_calc2 = NPP_calc.*(1+.001*age); %insufficient data
R_auto_calc3 = R_eco-R_het_soil; %insufficient data

% 1.2 group for plotting
% "in" (GPP components) 
in_flux_names = {'R_{auto}', 'BNPP', 'ANPP_{foliage}', 'ANPP_{woody}'};
in_fluxes = [R_auto_calc1; BNPP; ANPP_foliage;  ANPP_woody]; %matrix with all fluxes for stacked plot
in_sum = sum(in_fluxes,1);

% "out" (Reco components)
out_flux_names = {'R_{auto}', 'R_{het}'};
out_fluxes = [R_auto_calc1; R_het_soil]; %matrix with all fluxes for stacked plot
out_sum = sum(out_fluxes,1);

% all fluxes
flux_names= [in_flux_names, out_flux_names];
fluxes = [in_fluxes; -1*out_fluxes];


% 1.2. define mature fluxes 

%THIS IS CURRENTLY JUST THE 100-YEAR VALUE FOR FLUXES. NEEDS TO BE REPLACED WITH MATURE BIOME MEANS!

%for grouped stacked plot

mature_fluxes = [fluxes(:,end)'; fluxes(:,end)' ; fluxes(:,end)'; fluxes(:,end)'];



% 2 STOCKS
% 2.1 define age trends
% 2.1.1 equations (currently entered manually; to be automated)
biomass=54* log10(age); 
DW= 2.6*log10(age);

% 2.2 group for plotting
stock_names = {'B_{tot}', 'DW'};
stocks = [biomass; DW]; %matrix with all stocks for stacked plot


% 2. mature forests

mature_stocks = [stocks(:,end)'; stocks(:,end)' ; stocks(:,end)'; stocks(:,end)'];

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
plot(age, NEP, '-w', 'LineWidth', 3);hold on;
plot(age, NEP_calc, '--w', 'LineWidth', 3);hold on;
plot(age, -R_eco, '--r', 'LineWidth', 3);
t = title(biomes(b));
ylabel ('C  stocks (Mg C ha^{-1})')

subplot (2,4,4) %mature C fluxes
bar([mature_fluxes(b,:); 0*ones(1, size(mature_fluxes,2))], 'stacked');
legend (flux_names, 'Location', 'BestOutside');
xlim([0.5 1.5])
set(gca, 'XTick', []); %ticks off


subplot (2,4,5:7)
%flux plot:
h=area (age, stocks'); %, 'LineStyle','-'); 
hold on;
plot(age, biomass, '-k', 'LineWidth', 3);
xlabel ('stand age');
ylabel ('C  stocks (Mg C ha^{-1})')

subplot (2,4,8) %mature C stocks
bar([mature_stocks(b,:); 0*ones(1, size(mature_stocks,2))], 'stacked');
legend (stock_names, 'Location', 'BestOutside');
xlabel ('mature stands');
xlim([0.5 1.5])
set(gca, 'XTick', []); %ticks off


end