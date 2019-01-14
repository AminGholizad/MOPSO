function [REP]= mopso(c,iw,max_iter,lower_bound,upper_bound,swarm_size,rep_size,grid_size,alpha,beta,gamma,objective,constraint)
%mopso is an implementation of multi objective particle swarm optimization
%technique for a minimization problem
%   when called mopso() it solves the provided example or can be edited
%   here!
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
    objective=@fitness; % objective function
    constraint=@constraints; % constraints function
end
%% initialize particles
w = @(it) ((max_iter - it) - (iw(1) - iw(2)))/max_iter + iw(2);
pm = @(it) (1-(it-1)/(max_iter-1))^(1/mu);
swarm(1,swarm_size) = Particle();
for i =1:swarm_size
    swarm(i)=Particle(lower_bound,upper_bound,objective,constraint);
    while ~all(swarm(i).isFeasable)
        swarm(i)=Particle(lower_bound,upper_bound,objective,constraint);
    end
end
for i=1:swarm_size
    swarm(i) = swarm(i).updateDomination(swarm,i);
end
REP = Repository(swarm,rep_size,grid_size,alpha,beta,gamma);
%% Loop
for it=1:max_iter
    leader = REP.SelectLeader();
    wc = w(it); %current inertia weight
    pc=pm(it); %current mutation rate
    for i =1:swarm_size %update particles
        swarm(i)=swarm(i).update(wc,c,pc,leader,objective,constraint);
    end
    REP = REP.update(swarm);
    figure(1)
    PlotCosts(swarm,REP.swarm)
    drawnow
    disp(['Iteration ' num2str(it) ': Number of Rep Members = ' num2str(length(REP.swarm))]);
end
end
function PlotCosts(swarm,rep)
i = 1;
while i<=length(swarm)
    if isnan(swarm(i).cost)
        swarm(i)=[];
    else
        i=i+1;
    end
end
LEG = {};
if ~isempty(swarm)
    pop_costs=cell2mat({swarm.cost}');
    plot(pop_costs(:,1),pop_costs(:,2),'ko')
    hold on
    LEG = {'Current SWARM'};
end
rep_costs=cell2mat({rep.cost}');
plot(rep_costs(:,1),rep_costs(:,2),'r*')
xlabel('1^{st} Objective')
ylabel('2^{nd} Objective')
grid on
hold off
legend([LEG ,'REPASITORY'])
end
