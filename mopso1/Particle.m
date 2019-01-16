classdef Particle
    properties
        x
        l
        u
        v
        cost
        infeasablity
        pBest
        pBestCost
        pBestinfeasablity
        GridIndex
        isDominated
    end
    methods
        function obj = Particle(lower,upper,problem)
            if nargin > 0
                obj.GridIndex = 0;
                obj.isDominated = false;
                obj.x = unifrnd(lower,upper);
                obj.l = lower;
                obj.u = upper;
                obj.v = zeros(1,max(length(lower),length(upper)));
                [obj.cost, obj.infeasablity] = problem(obj.x);
                obj.pBest = obj.x;
                obj.pBestCost = obj.cost;
                obj.pBestinfeasablity = obj.infeasablity;
            end
        end
        function obj = update(obj,w,c,pm,gBest,problem)
            obj = obj.updateV(w,c,gBest);
            obj = obj.updateX();
            [obj.cost, obj.infeasablity] = problem(obj.x);
            obj = obj.applyMutatation(pm,problem);
            obj = obj.updatePbest();
        end
        function obj = updateV(obj,w,c,gBest)
            obj.v = w.*obj.v + c(1).*rand.*(obj.pBest-obj.x) + c(2).*rand.*(gBest.x-obj.x);
        end
        function obj = updateX(obj)
            i=find(or(obj.x+obj.v>obj.u,obj.x+obj.v<obj.l));
            obj.v(i) = -obj.v(i);
            obj.x = max(min(obj.x+obj.v,obj.u),obj.l);
        end
        function obj = updatePbest(obj)
            if obj.infeasablity == 0
                if obj.pBestinfeasablity > 0
                    obj.pBest = obj.x;
                    obj.pBestCost = obj.cost;
                    obj.pBestinfeasablity = obj.infeasablity;
                elseif all(obj.pBestCost >= obj.cost) && any(obj.pBestCost > obj.cost)
                    obj.pBest = obj.x;
                    obj.pBestCost = obj.cost;
                    obj.pBestinfeasablity = obj.infeasablity;
                end
            else
                if obj.pBestinfeasablity >= obj.infeasablity
                    obj.pBest = obj.x;
                    obj.pBestCost = obj.cost;
                    obj.pBestinfeasablity = obj.infeasablity;
                end
            end
        end
        function obj = applyMutatation(obj,pm,problem)
            if rand<pm
                X=obj.Mutate(pm);
                [X.cost,X.infeasablity]=problem(X.x);
                if X.dominates(obj)
                    obj=X;
                elseif ~obj.dominates(X)
                    if rand<0.5
                        obj=X;
                    end
                end
            end
        end
        function obj=Mutate(obj,pm)
            nVar=numel(obj.x);
            j=randi([1 nVar]);
            dx=pm*(obj.u(j)-obj.l(j));
            lb=max(obj.x(j)-dx,obj.l(j));
            ub=min(obj.x(j)+dx,obj.u(j));
            obj.x(j)=unifrnd(lb,ub);
        end
        function d = dominates(obj,obj1)
            if obj.infeasablity == 0
                if obj1.infeasablity == 0
                    if all(obj.cost <= obj1.cost) &&  any(obj.cost < obj1.cost)
                        d = true;
                    else
                        d = false;
                    end
                else
                    d = true;
                end
            elseif obj1.infeasablity == 0
                d = false;
            elseif obj.infeasablity < obj1.infeasablity
                d = true;
            else
                d = false;
            end
        end
        function obj=updateGridIndex(obj,Grid)
            nObj=length(obj.cost);
            nGrid=length(Grid(1).LB);
            GridSubIndex=zeros(1,nObj);
            for j=1:nObj
                GridSubIndex(j)=find(obj.cost(j)<Grid(j).UB,1,'first');
            end
            obj.GridIndex=GridSubIndex(1);
            for j=2:nObj
                obj.GridIndex=obj.GridIndex-1;
                obj.GridIndex=nGrid*obj.GridIndex;
                obj.GridIndex=obj.GridIndex+GridSubIndex(j);
            end
        end
    end
    methods (Static)
        function swarm = updateDomination(swarm)
            for index = 1:length(swarm)
            swarm(index).isDominated = false;
                for i = 1:length(swarm)
                    if i == index
                        continue
                    end
                    if swarm(i).dominates(swarm(index))
                        swarm(index).isDominated = true;
                        break
                    end
                end
            end
        end
    end
end

