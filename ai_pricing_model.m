%Replication code for the theoretical model and analysis in The Rise of AI Pricing: Trends, Driving Forces, and Implications for Firm Performance (2026 JME)


clear all


%parameter values
alpha = .6;
gamma = .2;
betta = .75;
A=.18;
%This approach finds zbar given phi:
Phi = 1; kappa = 1; eta = .00001; zbar = eta*kappa + sqrt(4*eta*Phi);
rho = 1;
chi = .085;
xi = 1/gamma; mu_lb = .15; 


%functions *conditional on AI adoption*
%basic labor choice (eqn 12) (holds for non-adopters too)
Lb_fn = @(mu,q,w) (mu*Phi*betta/w)^(1/(1-betta));
%AI labor choice (eqn 23)
La_fn = @(mu,q,w) (mu*Phi*(alpha/w)^(1-gamma)*A^alpha *(gamma/q)^gamma)^(1/(1-alpha-gamma));
%Computing choice (eqn 21)
C_fn = @(mu,q,w) (w/q)*(gamma/alpha)*La_fn(mu,q,w);
%The labor ratio Lb/La:
%ratio_fn = @(mu,q,w) La_fn(mu,q,w)/Lb_fn(mu,q,w);
share_fn = @(mu,q,w) La_fn(mu,q,w)/(La_fn(mu,q,w)+Lb_fn(mu,q,w));
%number of factors N:
Nfactors_ai_fn = @(mu,q,w) Lb_fn(mu,q,w)^betta + A^alpha * (La_fn(mu,q,w))^alpha * (C_fn(mu,q,w))^gamma;

%functions *conditional on non-adoption*
Nfactors_noai_fn = @(mu,q,w) Lb_fn(mu,q,w)^betta;

%adoption decision functions
%mmarket size threshold:
minmu_fn = @(q,w) 1 / (Phi * rho * A^alpha) * (w/alpha)^alpha * (q/gamma)^gamma  * chi / (1 - (alpha+gamma))^(1-(alpha+gamma));
adopt_fn = @(q,w) min((mu_lb/minmu_fn(q,w))^xi,1);

%General functions (both cases):
%how many factors?
Nfactors_fn = @(mu,q,w) (mu>=minmu_fn(q,w))*Nfactors_ai_fn(mu,q,w) + (mu<minmu_fn(q,w))*Nfactors_noai_fn(mu,q,w);
%revenue
y_fn = @(mu,q,w) mu*(rho*Nfactors_fn(mu,q,w) + zbar^2 - eta^2*kappa^2)/(4*eta); 
%markup
mark_fn = @(mu,q,w) y_fn(mu,q,w)/(kappa*mu*(zbar-eta*kappa)/2)-1;


% Importing Min's data:

% Import the time series data
[~, ~, raw] = xlsread('Moments_for_Model.xlsx','Time-Series','A2:E15');
% Create output variable
data = reshape([raw{:}],size(raw));
% Allocate imported array to column variable names
Year = data(:,1);
AIPricingShare = data(:,2);
PricingShare = data(:,3);
AdoptionRate = data(:,4);
CostofAIUsage = data(:,5);

% Import the cross-sectional data
opts = spreadsheetImportOptions("NumVariables", 8);
% Specify sheet and range
opts.Sheet = "Cross-Section-2023";
opts.DataRange = "A1:H21";
% Specify column names and types
opts.VariableNames = ["BySalesGroup1", "logsales_avg", "g_share", "g_adoption", "gcount", "gadopt", "pcount", "apcount"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double"];
% Import the data
tbl = readtable("Moments_for_Model.xlsx", opts, "UseExcel", false);
% Convert to output type
BySalesGroup1 = tbl.BySalesGroup1;
logsales_avg = tbl.logsales_avg;
g_share = tbl.g_share;
g_adoption = tbl.g_adoption;
gcount = tbl.gcount;
gadopt = tbl.gadopt;
pcount = tbl.pcount;
apcount = tbl.apcount;
% Clear temporary variables
clear opts tbl

%Calculate q trend:
Xregmat = [Year ones(length(Year),1)];
reg_coeffs = (Xregmat'*Xregmat)^-1 *Xregmat'*log(CostofAIUsage);
q_trend = exp(reg_coeffs(2)+reg_coeffs(1)*Year);
%calculate adoption/q elasticity:
Xregmat = [log(q_trend) ones(length(Year),1)];
reg_coeffs = (Xregmat'*Xregmat)^-1 *Xregmat'*log(AdoptionRate);
adopt_elast= reg_coeffs(1); 

%Market size distribution:
data_msize = exp((0:.2:3)-2.5);

% Clear temporary variables
clearvars data raw;

%Solve model

share_qvec = 0*Year;
adopt_qvec = 0*Year;
share_msize_2010 = 0*data_msize;
revenue_msize_2010 = 0*data_msize;
markup_msize_2010 = 0*data_msize;
share_msize_2023 = 0*data_msize;
revenue_msize_2023 = 0*data_msize;
markup_msize_2023 = 0*data_msize;


for qq = 1:length(Year)
    share_qvec(qq) = share_fn(1,q_trend(qq),1); 
    adopt_qvec(qq) = adopt_fn(q_trend(qq),1);
end
for qq = 1:length(data_msize)
    share_msize_2010(qq) = share_fn(data_msize(qq),q_trend(1),1)*(data_msize(qq)>minmu_fn(q_trend(1),1));
    revenue_msize_2010(qq) = y_fn(data_msize(qq),q_trend(1),1);
    markup_msize_2010(qq) = mark_fn(data_msize(qq),q_trend(1),1);
    share_msize_2023(qq) = share_fn(data_msize(qq),q_trend(end),1)*(data_msize(qq)>minmu_fn(q_trend(end),1));
    revenue_msize_2023(qq) = y_fn(data_msize(qq),q_trend(end),1);
    markup_msize_2023(qq) = mark_fn(data_msize(qq),q_trend(end),1);    
end

share_msize_Delta = share_msize_2023 - share_msize_2010;
revenue_msize_Delta = revenue_msize_2023 - revenue_msize_2010;
markup_msize_Delta = markup_msize_2023 - markup_msize_2010;


%empirical predictions:
logrevenue_OLS = 1.2*share_msize_Delta;
logmarkup_OLS = .3*share_msize_Delta;
%choose intercept for plots:
logrevenue_OLS = logrevenue_OLS - mean(logrevenue_OLS) + mean(log(revenue_msize_Delta));
logmarkup_OLS = logmarkup_OLS - mean(logmarkup_OLS) + mean(log(markup_msize_Delta));

%% time plots

saveplots = 1;

close all
fig1 = figure(1);
hold on
plot(Year,q_trend,'LineWidth',2,'DisplayName','Trend')
plot(Year,CostofAIUsage,'LineWidth',2,'LineStyle','--','DisplayName','Data')
legend show
legend('Location','N')
hold off
xlabel('Year');
ylabel('Computing Price')
    set(gca, ...
      'Box'         , 'off'     , ...
      'TickDir'     , 'out'     , ...
      'XMinorTick'  , 'off'      , ...
      'LineWidth'   , 1         )
 fontsize(12,'points')
 

fig2 = figure(2);
hold on
plot(Year,share_qvec,'LineWidth',2,'DisplayName','Model')
plot(Year,AIPricingShare,'LineWidth',2,'LineStyle','--','DisplayName','Data')
legend show
legend('Location','N')
hold off
xlabel('Year');
ylabel('AI Share')
    set(gca, ...
      'Box'         , 'off'     , ...
      'TickDir'     , 'out'     , ...
      'XMinorTick'  , 'off'      , ...
      'LineWidth'   , 1         )
 fontsize(12,'points')

fig3 = figure(3);
hold on
plot(Year,adopt_qvec,'LineWidth',2,'DisplayName','Model')
plot(Year,AdoptionRate,'LineWidth',2,'LineStyle','--','DisplayName','Data')
legend show
legend('Location','N')
hold off
xlabel('Year');
ylabel('Adoption')
    set(gca, ...
      'Box'         , 'off'     , ...
      'TickDir'     , 'out'     , ...
      'XMinorTick'  , 'off'      , ...
      'LineWidth'   , 1         )
 fontsize(12,'points')



 if saveplots == 1
 saveas(figure(1),'figures/computing_price.png')
 saveas(figure(2),'figures/pricing_labor_share.png')
 saveas(figure(3),'figures/pricing_adoption.png')
 end

 
% cross-section plots

fig4 = figure(4);
hold on
plot(log(revenue_msize_2023(2:end)),share_msize_2023(2:end),'LineWidth',2,'DisplayName','Model')
plot(logsales_avg,g_share,'LineWidth',2,'LineStyle','--','DisplayName','Data')
legend show
legend('Location','N')
hold off
xlabel('Log Sales');
ylabel('AI Share');

    set(gca, ...
      'Box'         , 'off'     , ...
      'TickDir'     , 'out'     , ...
      'XMinorTick'  , 'off'      , ...
      'LineWidth'   , 1         )
 fontsize(12,'points')

  if saveplots == 1
 saveas(figure(4),'figures/share_vs_sales.png')
 end
