% ------------------------------ BEGIN CODE -------------------------------

% System Dynamics ---------------------------------------------------------

boat = @(x,u) [cos(x(3)); sin(x(3)); 0];
sys = nonlinearSys('boat',boat);


% Parameters --------------------------------------------------------------

params.tFinal = 20.0;
R0 = interval([-1;-1;-1+pi/4],[1;1;1+pi/4]);


% Reachability Settings ---------------------------------------------------

% options_outer.alg = 'lin-adaptive';

% options_inner_Minkdiff.algInner = 'minkdiff';
% options_inner_Minkdiff.timeStep = 0.01;
% options_inner_Minkdiff.compOutputSet = false;
% options_inner_Minkdiff.tensorOrder = 2;

options1.algInner = 'scale';
options1.splits = 2;
options1.iter = 2;
options1.orderInner = 5;
options1.scaleFac = 0.95;
options1.timeStep = 0.01;                           
options1.taylorTerms = 10;                            
options1.zonotopeOrder = 50;       
options1.intermediateOrder = 20;
options1.errorOrder = 10;



% Reachability Analysis ---------------------------------------------------

start = tic();

params.R0 = zonotope(R0);
Rout = reach(sys,params,options_outer);

params.R0 = polytope(R0);
Rin = reachInner(sys,params,options1);
P = Rin.timePoint.set{end};
isempty(P)
% Rout.timePoint.set{end}
ellapsed = toc(start)

% Visualization -----------------------------------------------------------

% figure; hold on; box on;
% xlabel('x_1'); ylabel('x_2');
% 
% useCORAcolors("CORA:manual");
% plot(Rout.timePoint.set{end},[1,2],'DisplayName','Outer approximation');
% plot(Rin.timePoint.set{end},[1,2],'DisplayName','Inner approximation');

% xlim([0.65 0.75]); ylim([0.485 0.585]);
% xlim([-10 30]); ylim([-10 30]);

% ------------------------------ END OF CODE ------------------------------