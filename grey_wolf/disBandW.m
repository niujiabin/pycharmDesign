
%返回值 是最新的 Positions 狼的相对位置已经改变了·
function [ Positions,Alpha_pos,Alpha_score ] = disBandW(Positions,fobj,SearchAgents_no,dim,Alpha_pos,Alpha_score,BestPosition,K)

    %初始化 distance 用于 记录   最优距离
    distance = zeros(SearchAgents_no,2);

    %对所有狼和最优进行运算   最优狼不会处理
    for i = 1: SearchAgents_no
        
        %第i个狼的欧氏距离
        for j = 1:dim 
        
            distance(i,1) = distance(i,1) +  (Positions(i,j) - Alpha_pos(1,j)).^2;
        
        end
        
        %获取欧氏距离
        distance(i,1) = sqrt(distance(i,1));
        
        %第i个狼的序号
        distance(i,2) = i;
    end    
    %根据欧氏距离进行排序
    distance = sortrows(distance,1);   
    %获取前K个狼
    %先看最优狼    
    alfitpos  = Alpha_pos * unifrnd(-2,2);
    alfitness = fobj(alfitpos);
    
   % if alfitness < Alpha_score       
       %这里应该记录最优值的原位置 所以这里出错了  最优位置如果小于不发生
    
    Positions(BestPosition,:) =  alfitpos;
    Alpha_score = alfitness;
    
    
    % end
%     if als < Alpha_score
%         
%         Position
%     end
        
    for i=2:K+1       
        %获取狼的位置
        tempP = distance(i,2);
        %获取新的由随机数产生的新的个体狼的位置       
        fitpre =  fobj(Positions(tempP,:));
        
        temP2 = Positions(tempP,:) + distance(i,1) *  unifrnd(-2,2);        
        %获取新位置的个体值
        fitnow = fobj(temP2);
        
        %判断随机更改后的个体值 和 原来的做对比             
        if fitnow < fitpre           
            %如果比原来好，就更改
            Positions(i,:) = temP2;                     
        end
        %没有原来的好 就不做更改
    end
    
    
end

