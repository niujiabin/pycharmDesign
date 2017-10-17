function [ cg_curve ] = PSO( N,Max_iteration,lb,ub,dim,fobj )
%PSO 此处显示有关此函数的摘要
%   此处显示详细说明

%Particual Swarm Optimization
%粒子群算法的相关信息

%粒子在每一位空间上的每一维都有一个最大速度Vdmax，用来对粒子的速度进行限制
Vmax = 6;%最大速度
noP = N;%粒子的数量
wMax = 0.9;
wMin = 0.2;

%初始化 合适的c1 c2可以加快收敛速度 又不容易陷入局部最优

%学习因子 ，或称加速系数 

%分别调节pBest和Gbest方向飞行的最大步长
%决定粒子个体经验和群体经验对粒子运行轨迹的影响，反应离子群之间的信息交流

%如果c1=0 粒子只有群体经验 收敛速度快 但容易陷入局部最优
%c2=0 粒子没有群体信息 粒子各行其是 得到解得概率比较小 因此一般设置c1=c2
%这样就有了相同的影响力 较低的c1 c2是的粒子徘徊在远离目标的区域 较高的c1和c2
%值产生的运动或越过目标区域
c1=2;
c2=2;



iter = Max_iteration;
vel = zeros(noP,dim);
pBestScore = zeros(noP);%默认是一列
pBest = zeros(noP,dim);%一个例子本身所找到的最好解，即个体极值，另一个极值是整个粒子群中所有粒子在历代
%搜索中所达到的最优解纪委全局极值。
%找到两个最好解后，接下来是Pso中最重要的加速过程 每个粒子不断改变其在解空间中的速度，
%以尽可能地朝pBest和gBest所指向的区域飞去。
gBest = zeros(1,dim);%当前最好解
cg_curve = zeros(1,iter);

%随机初始化agents
pos = initialization(noP,dim,ub,lb);


for i=1:noP
    pBestScore(i) = inf;%默认是最大值
end

gBestScore = inf;


for l=1:iter
    
    %是否超界限的判断
    Flag4ub = pos(i,:)>ub;
    Flag4lb = pos(i,:)<lb;
    
    pos(i,:) = (pos(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
    
    for i=1:size(pos,1)
        %计算目标函数值 对于每一个粒子
        fitness = fobj(pos(i,:));
        
        %每一个得最优解
        if (pBestScore(i) > fitness)
            pBestScore(i) = fitness;
            pBest(i,:) = pos(i,:);
        end
        
        %判断全局最优解
        if (gBestScore > fitness)
            gBestScore = fitness;
            gBest = pos(i,:);
        end
    end
    
    %更新PSO的w权值 随着迭代次数的改变而改变
    w = wMax-l*((wMax-wMin)/iter);
    
    %更新速度和粒子的位置
    for i=1:size(pos,1)
        for j=1:size(pos,2)
            vel(i,j) = w*vel(i,j)+c1*rand()*(pBest(i,j)-pos(i,j))+rand()*c2*(gBest(j)-pos(i,j));
            
            %不能超过最大速度
            if (vel(i,j)>Vmax)
                vel(i,j) = Vmax;
            end
            
            if(vel(i,j)<-Vmax)
                vel(i,j) = -Vmax;
            end
            
            pos(i,j)  =pos(i,j)+vel(i,j);
        end
    end
    cg_curve(l) = gBestScore;
end




end

