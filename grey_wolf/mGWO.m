
function [Alpha_score,Alpha_pos,Convergence_curve] = mGWO(SearchAgents_no,Max_iter,lb,ub,dim,fobj,K)

%初始化 alpha beta and delta_pos
Alpha_pos = zeros(1,dim);%生成1行dim列的alpha狼的位置
Alpha_score = inf;%change this to -inf for maximization problems

Beta_pos = zeros(1,dim);
Beta_score = inf;%change this to -inf for maximinzation problems

Delta_pos = zeros(1,dim);
Delta_score = inf;


BestPosition = 0;
%初始化代理位置

Positions = initialization(SearchAgents_no,dim,ub,lb);

Convergence_curve = zeros(1,Max_iter);

l=0; %循环计数器

%进入循环

while l<Max_iter
    
    %遍历所有狼的个数
    for i=1:size(Positions,1)
        
        %返回超过搜索界限的代理
        Flag4ub  = Positions(i,:) > ub;
        Flag4lb = Positions(i,:)<lb;
        
        %如果两个之一为true则结果为true ~(Flag4ub+Flag4lb)为0 前半部分则为0 不参与计算
        %只有后半部分参与计算 如果两者都超过了界限 则重新复制为界限值 的加和 成为了中间
        %如果只有其中一个界限超过了 那么只给一个界限 这个等式写的很好 应该记住
        Positions(i,:) = (Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;

        %计算适应度
        fitness = fobj(Positions(i,:));
        
       
      
        
        %更新alhpa beta delta
        
        %这里实现的稍微有一点问题 应该是大于
        if fitness < Alpha_score
            Alpha_score = fitness;%更新值
            Alpha_pos = Positions(i,:);  %更新位置   
            %重新设置当前最优位置
            BestPosition = i;
        end
        
        if fitness > Alpha_score && fitness < Beta_score
            Beta_score = fitness;
            Beta_pos = Positions(i,:);
        end
        
        
        if fitness > Alpha_score && fitness > Beta_score && fitness < Delta_score
            Delta_score = fitness;
            Delta_pos = Positions(i,:);
        end
    end         
    
        %%%选择 前K个最优狼的位置，计算和其他狼和最有狼的距离
        a = 2-l*((2)/Max_iter); % a decreases linearly from 2 to 0                
        
         %计算完每一个狼之后，对其中的K个狼进行位置更新
       
        
        
        %更新位置搜索代理 包括omega
        for i=1:size(Positions,1)
            
            for j=1:size(Positions,2)
                
                r1 = rand();
                r2 = rand();
                
                A1 = 2*a*r1-a;
                C1 = 2*r2;
                D_alpha = abs(C1*Alpha_pos(j)-Positions(i,j));
                X1 = Alpha_pos(j) - A1*D_alpha;
                
                r1 = rand();
                r2 = rand();
                
                A2 = 2*a*r1-a;
                C2 = 2*r2;
                
                D_beta = abs(C2*Beta_pos(j)-Positions(i,j));
                X2 = Beta_pos(j)-A2*D_beta;
                                              
                %omega的值
                A3=2*a*r1-a; % Equation (3.3)
                C3=2*r2; % Equation (3.4)

                D_delta=abs(C3*Delta_pos(j)-Positions(i,j)); % Equation (3.5)-part 3
                X3=Delta_pos(j)-A3*D_delta; % Equation (3.5)-part 3             
            
                Positions(i,j)=(X1+X2+X3)/3;% Equation (3.7)
                
            end           
            
        end        
        
      
        [ Positions , Alpha_pos , Alpha_score ] = disBandW(Positions,fobj,SearchAgents_no,dim,Alpha_pos,Alpha_score,BestPosition,10);
        
        l = l+1;
        Convergence_curve(l) = Alpha_score;

end

