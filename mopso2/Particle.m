classdef Particle
    properties
        x
        l
        u
        v
        cost
        isFeasable
        pBest
        pBestCost
        GridIndex
        isDominated
    end
    methods
        function obj = Particle(lower,upper,objective,constraint)
            if nargin > 0
                obj.GridIndex = 0;
                obj.isDominated = false;
                obj.x = unifrnd(lower,upper);
                obj.l = lower;
                obj.u = upper;
                obj.v = zeros(1,max(length(lower),length(upper)));
                obj.isFeasable = constraint(obj.x);
                if all(obj.isFeasable)
                    obj.cost = objective(obj.x);
                    obj.pBestCost = obj.cost;
                    obj.pBest = obj.x;
                else
                    obj.cost = NaN;
                    obj.pBestCost = NaN;
                    obj.pBest = NaN;
                end
            end
        end
        function obj = update(obj,w,c,pm,gBest,objective,constraint)
            obj = obj.updateV(w,c,gBest);
            obj = obj.updateX();
            obj.isFeasable = constraint(obj.x);
            if all(obj.isFeasable)
                obj.cost = objective(obj.x);
            else
                obj.cost = NaN;
            end
            obj = obj.applyMutatation(pm,objective);
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
            if all(obj.isFeasable)
                if isnan(obj.pBestCost)
                    obj.pBest = obj.x;
                    obj.pBestCost = obj.cost;
                elseif all(obj.pBestCost >= obj.cost) && any(obj.pBestCost > obj.cost)
                    obj.pBest = obj.x;
                    obj.pBestCost = obj.cost;
                end
            else
                if isnan(obj.pBestCost)
                    if rand<0.5
                        obj.pBest = obj.x;
                        obj.pBestCost = obj.cost;
                    end
                end
            end
        end
        function obj = applyMutatation(obj,pm,objective)
            if rand<pm
                X=obj.Mutate(pm);
                X.cost=objective(X.x);
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
            if all(obj.isFeasable)
                if all(obj1.isFeasable)
                    if all(obj.cost <= obj1.cost) &&  any(obj.cost < obj1.cost)
                        d = true;
                    else
                        d = false;
                    end
                else
                    d = true;
                end
            elseif all(obj1.isFeasable)
                d = false;
            elseif sum(obj.isFeasable)>=sum(obj1.isFeasable)
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
        function obj = updateDomination(obj,swarm,index)
            obj.isDominated = false;
            for i = 1:length(swarm)
                if i == index
                    continue
                end
                if swarm(i).dominates(obj)
                    obj.isDominated = true;
                    break
                end
            end
        end
    end
    methods (Static)
    end
end
