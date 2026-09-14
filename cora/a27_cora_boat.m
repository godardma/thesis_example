% System Dynamics ---------------------------------------------------------

boat = @(x,u) [cos(x(3)); sin(x(3)); 0];
sys = nonlinearSys('boat',boat);

% Parameters --------------------------------------------------------------

params.tFinal = 20.0;
R0 = interval([-1;-1;-1],[1;1;1]);


% Reachability Settings ---------------------------------------------------

% outer

options_outer.alg = 'lin-adaptive';

% inner

% minkdiff

options_inner_Minkdiff.algInner = 'minkdiff';
options_inner_Minkdiff.timeStep = 0.005;
options_inner_Minkdiff.compOutputSet = false;
options_inner_Minkdiff.tensorOrder = 2;

% scale

options_inner_scale.algInner = 'scale';
options_inner_scale.splits = 2;
options_inner_scale.iter = 2;
optoptions_inner_scaleions1.orderInner = 5;
options_inner_scale.scaleFac = 0.95;
options_inner_scale.timeStep = 0.05;                           
options_inner_scale.taylorTerms = 10;                            
options_inner_scale.zonotopeOrder = 50;       
options_inner_scale.intermediateOrder = 20;
options_inner_scale.errorOrder = 10;



% Reachability Analysis ---------------------------------------------------

start = tic();

% outer
params.R0 = zonotope(R0);
Rout = reach(sys,params,options_outer);

% inner
params.R0 = polytope(R0);

% choose method for inner
Rin = reachInner(sys,params,options_inner_Minkdiff);
% Rin = reachInner(sys,params,options_inner_scale);

P = Rin.timePoint.set{end};

if representsa(P, 'emptySet')
    disp('The inner reachable set is empty.');
else
    disp('The inner reachable set is not empty.');
end
ellapsed = toc(start);
fprintf('Elapsed time for reachability analysis: %.2f seconds\n', ellapsed);

% Visualization -----------------------------------------------------------

figure; hold on; box on;
xlabel('x_1'); ylabel('x_2');

useCORAcolors("CORA:manual");
% plot Rout timePoint every 5 second
times = cell2mat(Rout.timePoint.time);

for t = 0:5:20

    [~,idx] = min(abs(times - t));

    plot(Rout.timePoint.set{idx}, [1,2], ...
        'DisplayName', sprintf('t = %.1f', times(idx)));

end

xlim([-5 35]); ylim([-30 30]);

