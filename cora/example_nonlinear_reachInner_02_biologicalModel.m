% 
% function res = example_nonlinear_reachInner_02_biologicalModel
% example_nonlinear_reachInner_02_biologicalModel - example for the
%    computation of an inner approximation of the reachable set for
%    nonlinear dynamics, using the approaches in [1] and [2]
%
% Syntax:
%    res = example_nonlinear_reachInner_02_biologicalModel
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 
%
% References:
%    [1] M. Wetzlinger, A. Kulmburg, and M. Althoff. "Inner approximations
%        of reachable sets for nonlinear systems using the Minkowski
%        difference". IEEE Control Systems Letters, 2024.

% Authors:       Mark Wetzlinger
% Written:       17-December-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% System Dynamics ---------------------------------------------------------

sys = nonlinearSys('biologicalModel',@biologicalModel);

%% 


% Parameters --------------------------------------------------------------

params.tFinal = 0.2;
R0 = interval(0.99*ones(7,1),1.01*ones(7,1));

%% 

% Reachability Settings ---------------------------------------------------

options_outer.alg = 'lin-adaptive';

%% 

% 1. Minkdiff algorithm
options_inner_Minkdiff.algInner = 'minkdiff';
options_inner_Minkdiff.timeStep = 0.01;
options_inner_Minkdiff.tensorOrder = 2;
options_inner_Minkdiff.compOutputSet = false;

% Reachability Analysis ---------------------------------------------------

%% 
params.R0 = polytope(R0);
start = tic();
[Rin_Minkdiff,Rout_Minkdiff] = reachInner(sys,params,options_inner_Minkdiff);
ellapsed = toc(start)
%% 

% Visualization -----------------------------------------------------------

figure;
projDims = {[1,2],[3,4],[5,6],[1,7],[2,3],[4,5],[6,7],[1,3],[2,4]};

for p=1:length(projDims)

    % labels
    subplot(3,3,p); hold on; box on;
    xlabel("x_" + projDims{p}(1)); ylabel("x_" + projDims{p}(2));
    useCORAcolors('CORA:contDynamics',3);
    
    % outer approximation 
    
    % inner approximation
    plot(Rin_Minkdiff,projDims{p},'DisplayName','Inner approximation (Minkowski difference)');
    plot(Rin_Minkdiff.R0,projDims{p},'DisplayName','Initial set');
    plot(Rin_Minkdiff.timePoint.set{end},projDims{p},'DisplayName','Final set');
    plot(Rout_Minkdiff.timePoint.set{end},projDims{p},'DisplayName','Final set');
    % initial set

end

% completed
res = true;

% ------------------------------ END OF CODE ------------------------------
