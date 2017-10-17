%灰狼算法的劣势在于没有交流

%改进思路 随机一部分狼 跑的较快 复刻围捕合作行为  那些较优的狼步伐较大
function [Alpha_score,Alpha_pos,Convergence_curve] = cGWO(SearchAgents_no,Max_iter,lb,ub,dim,fobj)
%当前全局最优解

gBest = zeros(1,dim);
gBestScore = inf; 

BestPosition = 0;
pBestScore = zeros(SearchAgents_no);%默认是一列
pBest = zeros(SearchAgents_no,dim);

Vmax = 6;%最大速度
wMax = 0.9;
wMin = 0.2;

%初始化 alpha beta and delta_pos
Alpha_pos = zeros(1,dim);%生成1行dim列的alpha狼的位置
Alpha_score = inf;%change this to -inf for maximization problems

Beta_pos = zeros(1,dim);
Beta_score = inf;%change this to -inf for maximinzation problems

Delta_pos = zeros(1,dim);
Delta_score = inf;

%初始化代理位置
vel = zeros(SearchAgents_no,dim);

Positions = initialization(SearchAgents_no,dim,ub,lb);

Convergence_curve = zeros(1,Max_iter);

l=0; %循环计数器

%进入循环

twpw2 = zeros(1,dim);

for i=1:SearchAgents_no
    pBestScore(i) = inf;%默认是最大值
end

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
        
%          twpw = Positions(i,:);
        
%         %令每一浪执行p次随机游走
%         for p=1:3
%             
%             for j=1:size(Positions,2);
%                 
%                  %8是搜索方向 2是步长
%                  twpw2(1,j) = twpw(1,j) + sin(2*pi*p/3)*2;
%                  Positions(i,j) = twpw2(1,j);
%                  
%                  ft = fobj(Positions(i,:));
%                  if ft < fitness
%                      fitness = ft;                  
%                  else
%                      Positions(i,:) = twpw;
%                  end
%                  
%             end
%         end
%         
        %当计算每一个狼的时候
        %随机取一个狼的相反值
%         wolfs =  randi(SearchAgents_no,1,1);
%         第一
%         wolfa = wolfs(1,1);
% 
%         获取这两个狼的适应度
%         fia = fobj(Positions(wolfa,:));
%         
%         tempPosition =  (Positions(wolfa,:) + Positions(i,:))/2;
%         获取适应度值
%         fitnesstemp = fobj(tempPosition);        
%         
%         if fitnesstemp < fia && fitnesstemp < fitness          
%             fitness = fitnesstemp; 
%         end
        %更新alhpa beta delta
        
        %这里实现的稍微有一点问题 应该是大于
        if fitness < gBestScore            
            gBest = Positions(i,:);
            gBestScore = fitness;
            BestPosition = i;
        end
        
        if (pBestScore(i) > fitness)
            pBestScore(i) = fitness;
            pBest(i,:) = Positions(i,:);
        end
        
        if fitness < Alpha_score
            Alpha_score = fitness;%更新值
            Alpha_pos = Positions(i,:);  %更新位置                    
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
   %  x=0:0.1:600; 
%y=sigmf(x,[1/60 300]); 
    %改进1 对参数进行改进
    
    %利用sigmoid函数
       a =  -(1/(1+(exp(-(10/Max_iter)*(l-Max_iter/2)))))*2+2;
      % a = 2-2*(1/(exp(1)-1)*(exp(l/Max_iter)-1)); % a decreases linearly from 2 to 0
%        a = 2-l*((2)/Max_iter); % a decreases linearly from 2 to 0    
         w = wMax-l*((wMax-wMin)/Max_iter);
        %更新位置搜索代理 包括omega
        for i=1:size(Positions,1)
           
            for j=1:size(Positions,2)
                        
                r1 = rand();
                r2 = rand();
 
%                 A1 = 2*a*r1-a;
%                 C1 = 2*r2;
%                 D_alpha = C1*Alpha_pos(j)-Positions(i,j) + 2*r1*(2-a)*1*(pBest(i,j)-Positions(i,j));
%                 X1 =vel(i,j) +  Alpha_pos(j) - A1*D_alpha;
%                 
%                 
%                 vel1 = w*vel(i,j)+2*r1*(pBest(i,j)-Positions(i,j))+2*r2*(Alpha_pos(j)-Positions(i,j));
%                    
%                 vel2 = w*vel(i,j)+2*r1*(pBest(i,j)-Positions(i,j))+2*r2*(Beta_pos(j)-Positions(i,j));
%                 
%                 vel3 = w*vel(i,j)+2*r1*(pBest(i,j)-Positions(i,j))+2*r2*(Delta_pos(j)-Positions(i,j));
%                 
%                 if (vel(i,j)>Vmax)
%                     vel(i,j) = Vmax;
%                 end
% 
%                 if(vel(i,j)<-Vmax)
%                     vel(i,j) = -Vmax;
%                 end
%                 
%                 vel(i,j) = (vel1+vel2+vel3)/3;
%                 Positions(i,j)  =Positions(i,j)+vel(i,j);
                
                A1 = 2*a*r1-a;
                C1 = 2*r2;
                D_alpha = C1*Alpha_pos(j)-Positions(i,j) + 2*r1*(2-a)*1*(pBest(i,j)-Positions(i,j));
                X1 =vel(i,j) +  Alpha_pos(j) - A1*D_alpha;
                
                r1 = rand();
                r2 = rand();
                
                A2 = 2*a*r1-a;
                C2 = 2*r2;
                
                D_beta = C2*Beta_pos(j)-Positions(i,j) + 2*r1*(2-a)*1.3*(pBest(i,j)-Positions(i,j));
                X2 = vel(i,j) + Beta_pos(j)-A2*D_beta;
                
                r1 = rand();
                r2 = rand();
                
                %omega的值
                A3=2*a*r1-a; % Equation (3.3)
                C3=2*r2; % Equation (3.4)

                D_delta=C3*Delta_pos(j)-Positions(i,j) + 2*r1*(2-a)*1.5*(pBest(i,j)-Positions(i,j)); % Equation (3.5)-part 3
                X3= vel(i,j) + Delta_pos(j)-A3*D_delta; % Equation (3.5)-part 3             
             
                r1 = rand();
                r2 = rand();
                
                %全局最优解已知是Alpha_score，再加一个他以前的最优解，
                A4=2*a*r1-a; % Equation (3.3)
                C4=2*r2; % Equation (3.4)
      
                P_delta=C4*pBest(i,j)-Positions(i,j) + 2*r1*1.7*(pBest(i,j)-Positions(i,j)); % Equation (3.5)-part 3
                X4= vel(i,j) + pBest(i,j)-A4*P_delta; % Equation (3.5)-part 3 
                
                %再来一
                %全局最优解已知是Alpha_score，再加一个他以前的最优解，
                A5=2*a*r1-a; % Equation (3.3)
                C5=2*r2; % Equation (3.4)
              
                %不仅要考虑三头头浪的 还要考虑当前狼本身
                Positions(i,j)=(X1+X2+X3+X4)/4 ;% Equation (3.7)
                  
                 
%                 Positions(i,j) =  rand()*abs(Positions(i,j)/(vel(i,j) + Positions(i,j))) * Positions(i,j) + rand()*abs(vel(i,j)/(vel(i,j) + Positions(i,j))) * vel(i,j);
                         
            end           
        end       
        
        %对每一个个体狼引进行进速度的概念
       % [ Positions , Alpha_pos , Alpha_score ] = disBandW(Positions,fobj,SearchAgents_no,dim,Alpha_pos,Alpha_score,BestPosition,3);
        
         
        orderFitnessPosition = orderFitness(Positions,fobj);
        
        %对一半的狼进行信息更改
        for i=fix(SearchAgents_no/2):SearchAgents_no
            
            r = randperm(5,1);
            
            if r==1
            
                 %被排序的狼
                 %是加还是减
                 si = rand();
                 if si>=0.5
                     Positions(orderFitnessPosition(i,2),:) = Alpha_pos + (ub-lb*rand()+lb);     
                 else
                      Positions(orderFitnessPosition(i,2),:) = Alpha_pos - (ub-lb*rand()+lb);
                 end
                               
            elseif r==2
                
                 si = rand();
                 if si>=0.5
                     Positions(orderFitnessPosition(i,2),:) = Beta_pos + (ub-lb*rand()+lb);     
                 else
                      Positions(orderFitnessPosition(i,2),:) = Beta_pos - (ub-lb*rand()+lb);
                 end
            elseif r==3
                       
                 si = rand();
                 if si>=0.5
                     Positions(orderFitnessPosition(i,2),:) = Delta_pos + (ub-lb*rand()+lb);     
                 else
                      Positions(orderFitnessPosition(i,2),:) = Delta_pos - (ub-lb*rand()+lb);
                 end
                 
             elseif r==4
                       
                 si = rand();
                 if si>=0.5
                     Positions(orderFitnessPosition(i,2),:) = (ub-lb*rand()+lb);     
                 else
                      Positions(orderFitnessPosition(i,2),:) =  (ub-lb*rand()+lb);
                 end
             elseif r==5
                       
                 si = rand();
                 if si>=0.5
                     Positions(orderFitnessPosition(i,2),:) = pBest(i,:) + (ub-lb*rand()+lb);     
                 else
                      Positions(orderFitnessPosition(i,2),:) = pBest(i,:) - (ub-lb*rand()+lb);
                 end
                 
            end
            
            %等概率的随机产生1到四的整数
        end
        l = l+1;
        Convergence_curve(l) = gBestScore;
end

