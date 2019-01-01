classdef Repository
    properties
        swarm
        rep_size
        Grid
        grid_size
        alpha
        beta
        gamma
    end
    
    methods
        function obj = Repository(swarm,rep_size,grid_size,alpha,beta,gamma)
            if nargin>0
                obj.rep_size = rep_size;
                obj.swarm = swarm(~[swarm.isDominated]);
                obj.grid_size=grid_size;
                obj.alpha=alpha;
                obj.beta = beta;
                obj.gamma = gamma;
                obj.Grid=obj.grid();
                for i = 1:length(obj.swarm)
                    obj.swarm(i) = obj.swarm(i).updateGridIndex(obj.Grid);
                end
            end
        end
        function Grid = grid(obj)
            C = cell2mat({obj.swarm.cost}');
            cmin = min(C,[],2);
            cmax = max(C,[],2);
            dc = cmax - cmin;
            cmin = cmin - obj.alpha * dc;
            cmax = cmax + obj.alpha * dc;
            nObj = size(C,2);
            empty_grid.LB = [];
            empty_grid.UB = [];
            Grid = repmat(empty_grid,nObj,1);
            for j = 1:nObj
                cj = linspace(cmin(j),cmax(j),obj.grid_size+1);
                Grid(j).LB = [-inf, cj];
                Grid(j).UB = [cj, +inf];
            end
        end
        function leader = SelectLeader(obj)
            GI = [obj.swarm.GridIndex];
            OC = unique(GI);
            N = zeros(size(OC));
            for k = 1:length(OC)
                N(k) = length(find(GI==OC(k)));
            end
            P = exp(-obj.beta*N);
            P = P/sum(P);
            sci = Repository.RouletteWheelSelection(P);
            sc = OC(sci);
            SCM = find(GI==sc);
            smi = randi([1 length(SCM)]);
            sm = SCM(smi);
            leader = obj.swarm(sm);
        end
        function obj = DeleteOneRepMemebr(obj)
            GI=[obj.swarm.GridIndex];
            OC=unique(GI);
            N=zeros(size(OC));
            for k=1:length(OC)
                N(k)=length(find(GI==OC(k)));
            end
            P=exp(obj.gamma*N);
            P=P/sum(P);
            sci=Repository.RouletteWheelSelection(P);
            sc=OC(sci);
            SCM=find(GI==sc);
            smi=randi([1 length(SCM)]);
            sm=SCM(smi);
            obj.swarm(sm)=[];
        end
        function obj = update(obj,swarm)
            obj.swarm = [obj.swarm,swarm(~[swarm.isDominated])];
            for i=1:length(obj.swarm)
                obj.swarm(i) = obj.swarm(i).updateDomination(obj.swarm,i);
            end
            obj.swarm = obj.swarm(~[obj.swarm.isDominated]);
            obj.Grid=obj.grid();
            for i = 1:length(obj.swarm)
                obj.swarm(i) = obj.swarm(i).updateGridIndex(obj.Grid);
            end
            Extra=length(obj.swarm)-obj.rep_size;
            for e=1:Extra
                obj=obj.DeleteOneRepMemebr();
            end
        end
    end
    methods (Static)
        function i = RouletteWheelSelection(P)
            i = find(rand<=cumsum(P),1,'first');
        end
    end
end

