% Quadratic programming to estimate fluxes in metabolic network
% MOMA model for flux distributions
% Model details see below
% Contact: tonghao0605@gmail.com

%module load /apps/Modules/MATLAB/R2016b

addpath(genpath('/apps/MATLAB/tomlab'));
addpath(genpath('/apps/MATLAB/glpk-4.48'));
addpath(genpath('/apps/MATLAB/glpkmex'));
addpath(genpath('/apps/MATLAB/opencobra-cobratoolbox-7be8e9b'));
changeCobraSolver('glpk');

addpath('/home/tong/fluxGxE/opt');
cd /home/tong/fluxGxE/opt

cd /home/tong/fluxGxE/opt/cvx/
cvx_setup

cd /home/tong/fluxGxE/opt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % col0 flux distribution% % % % % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model = readCbModel('ArabidopsisCoreModel.xml');

carb = 6;
oxy = 85;
starch = 22;
suc = 43;

co_ratio = 10.7/3.72;
ss_ratio = 6.32/2.45;

model.S(end+1,[carb oxy]) = [1 -co_ratio];
model.b(end+1) = 0;
model.mets(end+1) = {'carb/oxy ratio'};
model.metNames(end+1) = {'carb/oxy ratio'};

model.S(end+1,[starch suc]) = [1 -ss_ratio];
model.b(end+1) = 0;
model.mets(end+1) = {'starch/suc ratio'};
model.metNames(end+1) = {'starch/suc ratio'};

sol = optimizeCbModel(model);
v=sol.x;

%%%%%% col0 flux %%%%%%
fluxc = v; 

%c = model.c.';
%zmax = fluxc(find(c==1),:); % max biomass in model with ratio (CO,SS)
%zmax = 0.003124746706309;
					%%% the max with col0 biomass function

idnzero = union(find(fluxc>1E-5),find(fluxc<-1E-5)); %% non-zero flux id
idzero = intersect(find(fluxc<1E-5),find(fluxc>-1E-5)); %% zero flux id

csvwrite('nonzeroid.csv',idnzero);
csvwrite('zeroid.csv',idzero);

save('fluxcol0.mat','fluxc')
csvwrite('fluxcol0.csv',fluxc);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % accession fluxes by constraining % % % % % 
% % % % % the biomass ratio % % % % % 
% % % % % biomass functions % % % % % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% zmax for all accession biomass function
%%%% ratio (CO,SS) free
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

biofunc = csvread('biofunallaccfinal.csv',0,0);

zmaxall=[];

for i=1:67,

model = readCbModel('ArabidopsisCoreModel.xml');

model.S(:,549) = biofunc(:,i);

S = full(model.S);
[Srow Scol] = size(S);
c = model.c.';
n = Scol;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = c;

c1 = S;

carb = 6;
oxy = 85;
starch = 22;
suc = 43;

ox = repelem(0,n);
ox(carb) = 1;
c3 = ox;

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

c41 = ccou;
c42 = ccol;

c51 = cssu;
c52 = cssl;

Aeq = vertcat(c1,c3,c41,c42,c51,c52);

r1 = repelem(0,Srow);
beq = [r1,0,0,0,0,0];

ctype = [repelem('S',Srow),'L','U','L','U','L'];
vartype = repelem('C',n);
sense = -1;

lb = [model.lb.'];
ub = [model.ub.'];

[xmin, fmin, status, extra] = glpk (f, Aeq, beq, lb, ub, ctype, vartype, sense);

 % max biomass in model with ratio (CO,SS)
						%%% the max with col biomass function
zmaxall = [zmaxall,fmin];
 
end

%aa = max(zmaxall);
%bb = mean(zmaxall);
%zmax = 0.003957206534872; %%% max overall accessions 
zmax = 0.003503620892060; %%% average overall accessions 
%zmax = zmax - 0.0001;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% nonzero flux only
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

aramodel = readCbModel('ArabidopsisCoreModel.xml');

fluxc = load('fluxcol0.mat','fluxc');
fluxc = fluxc.fluxc;

S = full(aramodel.S);
[Srow Scol] = size(S);
c = aramodel.c.';
n = Scol;


idnzero = csvread('nonzeroid.csv',0,0);

fluxc = fluxc(idnzero);
[n m] = size(fluxc);

S = S(:,idnzero);
[Srow Scol] = size(S);

c = c(idnzero);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% scale biomass value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%bioID = csvread('lineID68_ID.csv',0,2);
biom = csvread('araid67_final_new.csv',0,1);
biom = biom(:,1);

biomc_m = biom(15,:); %% Col0 in measurement
biomc_p = fluxc(find(c==1),:); %% Col0 in model

biomratio = biom*biomc_p/biomc_m;

%biomall = [biom;biomc_m];
%biommax = max(biomall); 
biommax = max(biomratio); % max biomass in measurement

delta = 0.00011;
biomsall = (zmax-delta)*biomratio/biommax; %% scale biomass by maximum

biomsall = round(biomsall,5);

SEM = std(biomsall)/sqrt(length(biomsall));           
ts = [-1.645,1.645]; %%% 90% bimass interval
%CI = mean(biomsall)+ts*SEM;                      
ee = ts*SEM;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% quadratic programming by cvx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = 1./fluxc.';
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

lb = [aramodel.lb(idnzero)];
ub = [aramodel.ub(idnzero)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

biofunc = csvread('biofunallaccfinal.csv',0,0);

fluxall = [];
optvalall = [];
SS = S;

for i=1:67,

%%%%%% biomass + error %%%%%%%%%%%%%%%%%%

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
	%0.001 <= c*x <= 0.002;

	ccou*x <= 0;
	ccol*x >= 0;
	cssu*x <= 0;
	cssl*x >= 0;

cvx_end


fluxall = [fluxall,x];
optvalall = [optvalall,cvx_optval];

end


csvwrite('AraCore_allacc_opt.csv',fluxall);

save('AraCore_allacc_opt.mat','fluxall');


%writeCbModel(model, 'fileName', 'AraAccessionModel.sbml','format','sbml')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


