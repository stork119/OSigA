function R = solvetrajMRE(name,Times,y0,obs,par_in)

[parn, parv, parnames] = textread([pwd, '/models/','/',name,'/'  ,name, '.par'], '%s %f %q');
addpath([pwd,'/models/','/',name,'/symbolic/']);
addpath([pwd, '/models','/',name]);


[parn, parv, parnames] = textread([pwd, '/models/','/',name,'/'  ,name, '.par'], '%s %f %q');
stoichiometry=load([pwd, '/models/','/',name,'/' ,name,'_stoich.txt']);
parv=par_in;
size_sto=size(stoichiometry); 
nvar=size_sto(1);
par=parv;
npar=length(parn);
all_equations=str2func([name,'_MRE']);
ext_nvar=2*nvar+(nvar-1)*nvar/2;
%x = linspace(init_T,init_T+freq*N,N);
tspan=[0,5+max(Times)];

mysolution=ode15s(all_equations,tspan,y0(1:nvar),[],parv);
traj=deval(mysolution,5+Times);

R=[Times;traj];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








