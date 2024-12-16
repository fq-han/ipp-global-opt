function [xk5,xk_hist,xkerror,nsamples,it] = my_PSO_optimization(example_idx,d, MaxIt,term_tol,xinit)
%% Problem Definition
[CostFunction,~,xex] = choose_example(1,1,d,example_idx);
nVar=d;            % Number of Decision Variables
VarSize=[1 nVar];   % Size of Decision Variables Matrix
VarMin=-10;         % Lower Bound of Variables
VarMax= 10;         % Upper Bound of Variables

xkerror = zeros(MaxIt,1); xk_hist = zeros(MaxIt,d);
%% PSO Parameters
nPop=d*40;        % Population Size (Swarm Size)
% PSO Parameters
w=1;            % Inertia Weight
wdamp=0.99;     % Inertia Weight Damping Ratio
c1=1.5;         % Personal Learning Coefficient
c2=2.0;         % Global Learning Coefficient
 
% Velocity Limits
VelMax=0.1*(VarMax-VarMin);
VelMin=-VelMax;
%% Initialization
empty_particle.Position=[];
empty_particle.Cost=[];
empty_particle.Velocity=[];
empty_particle.Best.Position=[];
empty_particle.Best.Cost=[];
particle=repmat(empty_particle,nPop,1);
GlobalBest.Cost=inf;
for i=1:nPop
    
    % Initialize Position
    particle(i).Position=unifrnd(VarMin,VarMax,VarSize);
    particle(i).Position=reshape(xinit,[1,nVar])+randn(VarSize);
    % Initialize Velocity
    particle(i).Velocity=zeros(VarSize);
    
    % Evaluation
    particle(i).Cost=CostFunction(particle(i).Position);
    
    % Update Personal Best
    particle(i).Best.Position=particle(i).Position;
    particle(i).Best.Cost=particle(i).Cost;
    
    % Update Global Best
    if particle(i).Best.Cost<GlobalBest.Cost
        
        GlobalBest=particle(i).Best;
        
    end
    
end
BestCost=zeros(MaxIt,1);
nsamples = nPop;
%% PSO Main Loop
for it=1:MaxIt
    % max_pos_change = 0;
    for i=1:nPop
        ori_pos = particle(i).Position;
        % Update Velocity
        particle(i).Velocity = w*particle(i).Velocity ...
            +c1*rand(VarSize).*(particle(i).Best.Position-particle(i).Position) ...
            +c2*rand(VarSize).*(GlobalBest.Position-particle(i).Position);
        
        % Apply Velocity Limits
        particle(i).Velocity = max(particle(i).Velocity,VelMin);
        particle(i).Velocity = min(particle(i).Velocity,VelMax);
        
        % Update Position
        particle(i).Position = particle(i).Position + particle(i).Velocity;
        
        % Velocity Mirror Effect
        IsOutside=(particle(i).Position<VarMin | particle(i).Position>VarMax);
        particle(i).Velocity(IsOutside)=-particle(i).Velocity(IsOutside);
        
        % Apply Position Limits
        particle(i).Position = max(particle(i).Position,VarMin);
        particle(i).Position = min(particle(i).Position,VarMax);
        
        % Evaluation
        particle(i).Cost = CostFunction(particle(i).Position);

        % Update Personal Best
        if particle(i).Cost<particle(i).Best.Cost
            
            particle(i).Best.Position=particle(i).Position;
            particle(i).Best.Cost=particle(i).Cost;
            
            % Update Global Best
            if particle(i).Best.Cost<GlobalBest.Cost
                GlobalBest=particle(i).Best;    
            end    
        end
        % max_pos_change = max(max_pos_change,norm(particle(i).Position-ori_pos,inf));
    end
     BestCost(it)=GlobalBest.Cost;
    nsamples = nsamples + nPop;
    % disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
    
    w=w*wdamp;
    xk_hist(it,:) = GlobalBest.Position;
    xkerror(it) = my_error_opt(GlobalBest.Position,xex);
    if xkerror(it) < term_tol
        xk5 = GlobalBest.Position;
        return
    end
    % if max_pos_change < term_tol*1e-2
    %     fprintf('PSO trapped');
    %     xk5 = GlobalBest.Position;
    %     return;
    % end
end
xk5 = GlobalBest.Position;
end 
