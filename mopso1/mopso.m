function [REP]= mopso(c,iw,max_iter,lower_bound,upper_bound,swarm_size,rep_size,grid_size,alpha,beta,gamma,mu,problem)
%mopso is an implementation of multi objective particle swarm optimization
%technique for a minimization problem
%% initialize parameters
if nargin==0
    c = [0.1,0.2]; % [cognitive acceleration, social acceleration] coefficients
    iw = [0.5 0.001]; % [starting, ending] inertia weight
    max_iter = 100; % maximum iterations
    lower_bound = zeros(1,3); % lower bound of vars
    upper_bound = pi/2*ones(1,3); % upper bound of vars
    swarm_size=100; % swarm size
    rep_size=100; % Repository Size
    grid_size=7; % Number of Grids per Dimension
    alpha=0.1; % Inflation Rate
    beta=2; % Leader Selection Pressure
    gamma=2; % Deletion Selection Pressure
    mu=0.1; % Mutation Rate
    problem=@prob; % objective function
end
%% initialize particles
fprintf('Initializing swarm ...\n')
w = @(it) ((max_iter - it) - (iw(1) - iw(2)))/max_iter + iw(2);
pm = @(it) (1-(it-1)/(max_iter-1))^(1/mu);
swarm(1,swarm_size) = Particle();
for i = 1:swarm_size
    swarm(i)=Particle(lower_bound,upper_bound,problem);
    retry = 0;
    while swarm(i).infeasablity > 0 && retry < 100
        swarm(i)=Particle(lower_bound,upper_bound,problem);
        retry = retry + 1;
    end
end
REP = Repository(swarm,rep_size,grid_size,alpha,beta,gamma);
%% Loop
fprintf('Starting the optimization loop ...\n')
for it=1:max_iter
    leader = REP.SelectLeader();
    wc = w(it); %current inertia weight
    pc=pm(it); %current mutation rate
    for i =1:swarm_size %update particles
        swarm(i)=swarm(i).update(wc,c,pc,leader,problem);
    end
    REP = REP.update(swarm);
    Title = sprintf('Iteration %d, Number of Rep Members = %d',it,length(REP.swarm));
    PlotCosts(swarm,REP.swarm,Title)
    disp(Title);
end
end
function PlotCosts(swarm,rep,Title)
figure(1)
feasable_swarm = swarm([swarm.infeasablity]==0);
infeasable_swarm = swarm([swarm.infeasablity]>0);
LEG = {};
if ~isempty(feasable_swarm)
    swarm_costs=vertcat(feasable_swarm.cost);
    plot(swarm_costs(:,1),swarm_costs(:,2),'go')
    hold on
    LEG = {'Current feasable SWARM'};
    Title = sprintf([Title '\nfeasable swarm=%d'],length(feasable_swarm));
end
if ~isempty(infeasable_swarm)
    swarm_costs=vertcat(infeasable_swarm.cost);
    plot(swarm_costs(:,1),swarm_costs(:,2),'ro')
    hold on
    LEG = [LEG, 'Current infeasable SWARM'];
    if contains(Title,newline)
        Title = sprintf([Title ', infeasable swarm=%d'],length(infeasable_swarm));
    else
        Title = sprintf([Title '\ninfeasable swarm=%d'],length(infeasable_swarm));
    end
end
rep_costs=vertcat(rep.cost);
plot(rep_costs(:,1),rep_costs(:,2),'b*')
xlabel('1^{st} Objective')
ylabel('2^{nd} Objective')
grid on
hold off
title(Title)
legend([LEG ,'REPASITORY'],'location','best')
drawnow
end