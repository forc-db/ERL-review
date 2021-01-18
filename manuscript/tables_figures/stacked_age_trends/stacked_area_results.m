% age trend summaries for ERL review
% Kristina Anderson Teixeira
% January 2021

clf; clear; close all;
set(0,'DefaultLegendAutoUpdate','off');

%% SETTINGS
figures_dir='/Users/kteixeira/Dropbox (Smithsonian)/GitHub/ForC-db/ERL-review/manuscript/tables_figures';
working_dir='/Users/kteixeira/Dropbox (Smithsonian)/GitHub/ForC-db/ERL-review/manuscript/tables_figures/stacked_age_trends';

%% READ IN & PREPARE DATA
%read in data:
cd(figures_dir)
model_summaries=readtable('SI_age_trend_model_summaries.csv');
%remove specieal formatting from model_summaries.Variable:
for n=1:size(model_summaries,1)
    s=cell2mat(model_summaries.Variable(n));
    s(regexp(s,'[$,{,}]'))=[];
    model_summaries.Variable(n)=cellstr(s);
end

%prepare some variables
variables = unique(model_summaries.Variable);
biomes= {'Tropical broadleaf','Temperate broadleaf', 'Temperate conifer', 'Boreal conifer'};
%{'BiomeTropical broadleaf','BiomeTemperate broadleaf','BiomeTemperate conifer','BiomeBoreal conifer'}
params_matrix=cell2table(num2cell(NaN*ones(length(variables), 10)),'VariableNames',{'variable', 'beta' 'betaTrB' 'betaTeB' 'betaTeN' 'betaBoN' 'intTrB' 'intTeB' 'intTeN' 'intBoN'});
params_matrix.variable=variables;

for v=1:length(variables) %cycle through variables and fill in params_matrix
    var=variables(v);
    %beta
    if max(strcmp(model_summaries.Variable, var) + strcmp(model_summaries.Parameter, 'log10(stand.age)'))==2;
        params_matrix.beta(v)=table2array(model_summaries(strcmp(model_summaries.Variable, var) & strcmp(model_summaries.Parameter, 'log10(stand.age)'),3));
    else
        params_matrix.beta(v)=NaN;
    end
    %betaTrB
    if max(strcmp(model_summaries.Variable, var) + strcmp(model_summaries.Parameter, 'log10(stand.age):BiomeTropical broadleaf'))==2;
        params_matrix.betaTrB(v)=table2array(model_summaries(strcmp(model_summaries.Variable, var) & strcmp(model_summaries.Parameter, 'log10(stand.age):BiomeTropical broadleaf'),3));
    else
        params_matrix.betaTrB(v)=0;
    end
    %betaTeB
    if max(strcmp(model_summaries.Variable, var) + strcmp(model_summaries.Parameter, 'log10(stand.age):BiomeTemperate broadleaf'))==2;
        params_matrix.betaTeB(v)=table2array(model_summaries(strcmp(model_summaries.Variable, var) & strcmp(model_summaries.Parameter, 'log10(stand.age):BiomeTemperate broadleaf'),3));
    else
        params_matrix.betaTeB(v)=0;
    end
    %betaTeN
    if max(strcmp(model_summaries.Variable, var) + strcmp(model_summaries.Parameter, 'log10(stand.age):BiomeTemperate conifer'))==2;
        params_matrix.betaTeN(v)=table2array(model_summaries(strcmp(model_summaries.Variable, var) & strcmp(model_summaries.Parameter, 'log10(stand.age):BiomeTemperate conifer'),3));
    else
        params_matrix.betaTeN(v)=0;
    end
    %betaBoN
    if max(strcmp(model_summaries.Variable, var) + strcmp(model_summaries.Parameter, 'log10(stand.age):BiomeBoreal conifer'))==2;
        params_matrix.betaBoN(v)=table2array(model_summaries(strcmp(model_summaries.Variable, var) & strcmp(model_summaries.Parameter, 'log10(stand.age):BiomeBoreal conifer'),3));
    else
        params_matrix.betaBoN(v)=0;
    end
    %intTrB
    if max(strcmp(model_summaries.Variable, var) + strcmp(model_summaries.Parameter, 'BiomeTropical broadleaf'))==2;
        params_matrix.intTrB(v)=table2array(model_summaries(strcmp(model_summaries.Variable, var) & strcmp(model_summaries.Parameter, 'BiomeTropical broadleaf'),3));
    else
        params_matrix.intTrB(v)=0;
    end
    %intTeB
    if max(strcmp(model_summaries.Variable, var) + strcmp(model_summaries.Parameter, 'BiomeTemperate broadleaf'))==2;
        params_matrix.intTeB(v)=table2array(model_summaries(strcmp(model_summaries.Variable, var) & strcmp(model_summaries.Parameter, 'BiomeTemperate broadleaf'),3));
    else
        params_matrix.intTeB(v)=0;
    end
    %intTeN
    if max(strcmp(model_summaries.Variable, var) + strcmp(model_summaries.Parameter, 'BiomeTemperate conifer'))==2;
        params_matrix.intTeN(v)=table2array(model_summaries(strcmp(model_summaries.Variable, var) & strcmp(model_summaries.Parameter, 'BiomeTemperate conifer'),3));
    else
        params_matrix.intTeN(v)=0;
    end
    %intBoN
    if max(strcmp(model_summaries.Variable, var) + strcmp(model_summaries.Parameter, 'BiomeBoreal conifer'))==2;
        params_matrix.intBoN(v)=table2array(model_summaries(strcmp(model_summaries.Variable, var) & strcmp(model_summaries.Parameter, 'BiomeBoreal conifer'),3));
    else
        params_matrix.intBoN(v)=0;
    end
end

int_matrix=[params_matrix.intTrB, params_matrix.intTeB, params_matrix.intTeN, params_matrix.intBoN];
beta_matrix=[params_matrix.betaTrB, params_matrix.betaTeB, params_matrix.betaTeN, params_matrix.betaBoN];

for b =1:4
%% GENERATE EQUATIONS FOR PLOTTING
age = linspace (1,100,100);


% 1 FLUXES
% 1.1 define age trends
% 1.1.1 equations 
NEP_index=find(strcmp(params_matrix.variable,'NEP'));
NEP=(params_matrix.beta(NEP_index)+beta_matrix(NEP_index,b))* log10(age)+ int_matrix(NEP_index,b);

GPP_index=find(strcmp(params_matrix.variable,'GPP')); %not available for tropical forests
GPP=max(0,(params_matrix.beta(GPP_index)+beta_matrix(GPP_index,b))* log10(age)+ int_matrix(GPP_index,b));

ANPP_index=find(strcmp(params_matrix.variable,'ANPP'));
ANPP=max(0,(params_matrix.beta(ANPP_index)+beta_matrix(ANPP_index,b))* log10(age)+ int_matrix(ANPP_index,b));

ANPP_foliage_index=find(strcmp(params_matrix.variable,'ANPP_foliage'));
ANPP_foliage=max(0,(params_matrix.beta(ANPP_foliage_index)+beta_matrix(ANPP_foliage_index,b))* log10(age)+ int_matrix(ANPP_foliage_index,b));

ANPP_stem_index=find(strcmp(params_matrix.variable,'ANPP_stem'));
ANPP_stem=max(0,(params_matrix.beta(ANPP_stem_index)+beta_matrix(ANPP_stem_index,b))* log10(age)+ int_matrix(ANPP_stem_index,b));

ANPP_woody_index=find(strcmp(params_matrix.variable,'ANPP_woody'));
ANPP_woody=max(0,(params_matrix.beta(ANPP_woody_index)+beta_matrix(ANPP_woody_index,b))* log10(age)+ int_matrix(ANPP_woody_index,b));

BNPP_index=find(strcmp(params_matrix.variable,'BNPP'));
BNPP=max(0,(params_matrix.beta(BNPP_index)+beta_matrix(BNPP_index,b))* log10(age)+ int_matrix(BNPP_index,b));

BNPP_fine_index=find(strcmp(params_matrix.variable,'BNPP_fine'));
BNPP_fine=max(0,(params_matrix.beta(BNPP_fine_index)+beta_matrix(BNPP_fine_index,b))* log10(age)+ int_matrix(BNPP_fine_index,b));


NPP_index=find(strcmp(params_matrix.variable,'NPP'));
NPP=max(0,(params_matrix.beta(NPP_index)+beta_matrix(NPP_index,b))* log10(age)+ int_matrix(NPP_index,b));

R_eco_index=find(strcmp(params_matrix.variable,'R_eco'));
R_eco=max(0,(params_matrix.beta(R_eco_index)+beta_matrix(R_eco_index,b))* log10(age)+ int_matrix(R_eco_index,b));


R_soil_index=find(strcmp(params_matrix.variable,'R_soil'));
R_soil=max(0,(params_matrix.beta(R_soil_index)+beta_matrix(R_soil_index,b))* log10(age)+ int_matrix(R_soil_index,b));

R_het_soil_index=find(strcmp(params_matrix.variable,'R_het-soil'));
R_het_soil=max(0,(params_matrix.beta(R_het_soil_index)+beta_matrix(R_het_soil_index,b))* log10(age)+ int_matrix(R_het_soil_index,b));

R_root_index=find(strcmp(params_matrix.variable,'R_root'));
R_root=max(0,(params_matrix.beta(R_root_index)+beta_matrix(R_root_index,b))* log10(age)+ int_matrix(R_root_index,b));


% 1.1.2 calculated fluxes
NEP_calc=GPP-R_eco;


%options for CUE. lat one in list is the one that will be used.
CUE=0.679-0.153*log10(age); %from DeLucia et al. 2007 (note that equation given mistakenly leaves of the log10(age))
CUE=0.5;

if b>1  %extratropical forests
    ANPP_woody_calc=min(ANPP_stem, ANPP_woody);
    NPP_calc=ANPP_foliage+ANPP_woody_calc+BNPP;
    R_auto_ag_calc=R_eco-R_soil;
    R_auto_calc=R_root+R_auto_ag_calc;
elseif b==1  %tropical forests
    ANPP_woody_calc=max(0,ANPP-ANPP_foliage);
    NPP_calc=ANPP_foliage+ANPP_woody_calc+BNPP;
    %options for R_auto_calc (insufficient data for regression or to calculate as sum/difference). 
    %Last one in list is the one that will be used.
    R_auto_calc = NPP_calc.*(1./CUE-1); %Rauto=NPP*(1/CUE-1), and CUE equation is above
    R_auto_ag_calc=max(0,R_auto_calc-R_root);
end





BNPP_coarse_calc=max(0,BNPP-BNPP_fine);

% 1.2 group for plotting
% "in" (GPP components) 
in_flux_names = {'R_{root}', 'R_{auto-ag}*', 'BNPP', 'ANPP_{foliage}', 'ANPP_{woody}*'};
in_fluxes = [R_root; R_auto_ag_calc; BNPP; ANPP_foliage; ANPP_woody_calc ]; %matrix with all fluxes for stacked plot
in_sum = sum(in_fluxes,1);

% "out" (Reco components)
out_flux_names = { 'R_{het-soil}', 'R_{root}', 'R_{auto-ag}*'};
out_fluxes = [R_het_soil; R_root; R_auto_ag_calc ]; %matrix with all fluxes for stacked plot
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
% 2.1.1 equations

B_tot_index=find(strcmp(params_matrix.variable,'B_tot'));
B_tot=max(0,(params_matrix.beta(B_tot_index)+beta_matrix(B_tot_index,b))* log10(age)+ int_matrix(B_tot_index,b));

B_ag_index=find(strcmp(params_matrix.variable,'B_ag'));
B_ag=max(0,(params_matrix.beta(B_ag_index)+beta_matrix(B_ag_index,b))* log10(age)+ int_matrix(B_ag_index,b));

B_root_index=find(strcmp(params_matrix.variable,'B_root'));
B_root=max(0,(params_matrix.beta(B_root_index)+beta_matrix(B_root_index,b))* log10(age)+ int_matrix(B_root_index,b));

B_root_fine_index=find(strcmp(params_matrix.variable,'B_root-fine'));
B_root_fine=max(0,(params_matrix.beta(B_root_fine_index)+beta_matrix(B_root_fine_index,b))* log10(age)+ int_matrix(B_root_fine_index,b));

B_foliage_index=find(strcmp(params_matrix.variable,'B_foliage'));
B_foliage=max(0,(params_matrix.beta(B_foliage_index)+beta_matrix(B_foliage_index,b))* log10(age)+ int_matrix(B_foliage_index,b));

DW_tot_index=find(strcmp(params_matrix.variable,'DW_tot'));
DW_tot=max(0,(params_matrix.beta(DW_tot_index)+beta_matrix(DW_tot_index,b))* log10(age)+ int_matrix(DW_tot_index,b));

DW_down_index=find(strcmp(params_matrix.variable,'DW_down'));
DW_down=max(0,(params_matrix.beta(DW_down_index)+beta_matrix(DW_down_index,b))* log10(age)+ int_matrix(DW_down_index,b));

OL_index=find(strcmp(params_matrix.variable,'OL'));
OL=max(0,(params_matrix.beta(OL_index)+beta_matrix(OL_index,b))* log10(age)+ int_matrix(OL_index,b));


% 2.1.2 caculated stocks
DW_standing=max(0,DW_tot-DW_down);
B_root_coarse=max(0,B_root-B_root_fine);
B_tot_calc=B_ag+B_root_coarse+B_root_fine;
B_ag_wood_calc=max(0,B_ag-B_foliage);

% 2.2 group for plotting
stock_names = { 'B_{root-coarse}*', 'B_{root-fine}','B_{ag-wood}*', 'B_{foliage}','DW_{standing}*', 'DW_{down}', 'OL'};
stocks = [ B_root_coarse ; B_root_fine ; B_ag_wood_calc;  B_foliage; DW_standing; DW_down; OL]; %matrix with all stocks for stacked plot


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
%plot(age, in_sum, '--b', 'LineWidth', 3); hold on;
if b~= 1 
    plot(age, GPP, '-b', 'LineWidth', 3); hold on; % eddy flux: insufficient data for tropics
    plot(age, -R_eco, '-r', 'LineWidth', 3); % eddy flux: insufficient data for tropics
    %plot(age, -R_root-R_het_soil-R_auto_ag_calc, '--r', 'LineWidth', 3); % eddy flux: insufficient data for tropics
    plot(age, NEP, '-w', 'LineWidth', 3);hold on; % eddy flux: insufficient data for tropics
    plot(age, NEP_calc, '--w', 'LineWidth', 3);hold on;
end
plot(age, -R_soil, '-k', 'LineWidth', 3);hold on;
plot(age, -R_root-R_het_soil, '--k', 'LineWidth', 3);hold on;
t = title(biomes(b));
ylabel ('C  fluxes (Mg C ha^{-1})')

subplot (2,4,4) %mature C fluxes
bar([mature_fluxes(b,:); 0*ones(1, size(mature_fluxes,2))], 'stacked');
legend (flux_names, 'Location', 'BestOutside');
xlim([0.5 1.5])
set(gca, 'XTick', []); %ticks off


subplot (2,4,5:7)
%flux plot:
h=area (age, stocks'); hold on;
%compare biomass from different estimation methods:
%plot(age, B_tot, '-k', 'LineWidth', 3);
%plot(age, B_tot_calc, '--k', 'LineWidth', 3);
xlabel ('stand age (years)');
ylabel ('C  stocks (Mg C ha^{-1})')

subplot (2,4,8) %mature C stocks
bar([mature_stocks(b,:); 0*ones(1, size(mature_stocks,2))], 'stacked');
legend (stock_names, 'Location', 'BestOutside');
xlabel ('mature stands');
xlim([0.5 1.5])
set(gca, 'XTick', []); %ticks off

%figure (10+b)
%plot (age, ANPP, age, ANPP_foliage, age, ANPP_woody_calc, age, ANPP_foliage+ANPP_woody_calc)
%legend ('ANPP', 'ANPP_foliage', 'ANPP woody', 'foliage+woody')

%% ~~~~~~~ SAVE FIGURE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cd(working_dir)
savename=char(strcat(biomes(b),'_age_trends'));
print(savename, '-dpng', '-r300')

end