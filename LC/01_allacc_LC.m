% Quadratic programming to estimate fluxes in metabolic network
% MOMA model for flux distributions
% Model details see below
% Contact: tonghao0605@gmail.com

% module load /apps/Modules/devel/Matlab-2020a

addpath(genpath('/winmounts/tong/winhome/MATLAB/glpk-4.48'));
addpath(genpath('/winmounts/tong/winhome/MATLAB/glpkmex'));
addpath(genpath('/winmounts/tong/winhome/MATLAB/opencobra-cobratoolbox-7be8e9b'));
changeCobraSolver('glpk');

addpath('/winmounts/tong/winhome/mpidir/8.AlleleFlux/00_GxE/LC');
cd /winmounts/tong/winhome/mpidir/8.AlleleFlux/00_GxE/LC/cvx/
cvx_setup

cd /winmounts/tong/winhome/mpidir/8.AlleleFlux/00_GxE/LC


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% flux distribution in new environment using 
% condition specific biomass function and the biomass ratio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% non-zero flux only
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

aramodel = readCbModel('ArabidopsisCoreModel.xml');

fluxcall = load('AraCore_allacc_opt.mat','fluxall');
fluxcall = fluxcall.fluxall;

S = full(aramodel.S);
%[Srow Scol] = size(S);
c = aramodel.c.';
%n = Scol;

idnzero = csvread('nonzeroid.csv',0,0);

%fluxc = fluxc(idnzero);
[n m] = size(fluxcall);

S = S(:,idnzero);
[Srow Scol] = size(S);

c = c(idnzero);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% scale biomass value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

biom = csvread('araid67_final_new_env.csv');
biom1 = biom(:,2); %% high N % E1
biom2 = biom(:,4); %% low C % E2

biom_m = biom1; %% biomass in measurement E1
biom_p = fluxcall(find(c==1),:).'; %% biomass in model E1

biom_me = biom2; %% biomass in measurement E2

biomratio = biom_me.*biom_p./biom_m;
biomratio = round(biomratio,5); %% biomass in model E1

biomsall = biomratio;

%SEM = std(biomratio)/sqrt(length(biomratio));           
%ts = [-2.576,2.576]; %%% 99% bimass interval
%CI = mean(biomsall)+ts*SEM;                      
ee = [-0.00001,0.00001];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% quadratic programming by cvx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

biofunc = csvread('biofunallaccfinal_LC.csv',0,0);
fluxall = [];
optvalall = [];

for i=1:67,

A = 1./fluxcall(:,i).';
b = 1;

%A = eye(n);
%b = fluxc;

carb = 6;
oxy = 61;
starch = 20;
suc = 31;

ccou = repelem(0,n);
ccou(carb) = 1;
ccou(oxy) = -3.81445;

ccol = repelem(0,n);
ccol(carb) = 1;
ccol(oxy) = -0.93815;

cssu = repelem(0,n);
cssu(starch) = 1;
cssu(suc) = -3.3694;

cssl = repelem(0,n);
cssl(starch) = 1;
cssl(suc) = -0.7898;

lb = aramodel.lb(idnzero);
ub = aramodel.ub(idnzero);

%%%%%% biomass + error %%%%%%%%%%%%%%%%%%

SS = S;
SS(:,n) = biofunc(:,i);

x = [];

cvx_begin quiet

	variable x(n);
	minimize(norm(A*x-b));

	subject to

	SS*x == repelem(0,Srow).';	
	lb <= x <= ub;
	
	%c*x == 0.002;
	biomsall(i)+ee(1) <= c*x <= biomsall(i)+ee(2);
	%0.00001 <= c*x <= 0.004;

	ccou*x <= 0;
	ccol*x >= 0;
	cssu*x <= 0;
	cssl*x >= 0;

cvx_end

fluxall = [fluxall,x];
optvalall = [optvalall,cvx_optval];

end 

csvwrite('AraCore_allacc_LC.csv',fluxall);

save('AraCore_allacc_LC.mat','fluxall');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


