
%%灰狼算法的并行算法   2017-03-21 实现方法
%%Aparallel grey wolf optimizer 
function [Alpha_score,Alpha_pos,Convergence_curve] = paGWO(SearchAgents_no,Max_iter,lb,ub,dim,fobj)

%初始化 alpha beta and delta_pos
%第一个小组
FirstGroupAlpha_pos = zeros(1,dim);%生成1行dim列的alpha狼的位置
FirstGroupAlpha_score = inf;%change this to -inf for maximization problems

FirstGroupBeta_pos = zeros(1,dim);
FirstGroupBeta_score = inf;%change this to -inf for maximinzation problems

FirstGroupDelta_pos = zeros(1,dim);
FirstGroupDelta_score = inf;

%第二个小组
SecondGroupAlpha_pos = zeros(1,dim);%生成1行dim列的alpha狼的位置
SecondGroupAlpha_score = inf;%change this to -inf for maximization problems

SecondGroupBeta_pos = zeros(1,dim);
SecondGroupBeta_score = inf;%change this to -inf for maximinzation problems

SecondGroupDelta_pos = zeros(1,dim);
SecondGroupDelta_score = inf;

%初始化代理位置

Positions = initialization(SearchAgents_no,dim,ub,lb);

Convergence_curve = zeros(1,Max_iter);



%Generate initial population and divide them into 2 subgroups
%第一个小组的狼数
SubGroupWolfNumber = SearchAgents_no/2;
%第一个小组 取前SubGroupWolfNumber行
SubGroupsFirstPosition = Positions(1:SubGroupWolfNumber,:);
%第二个小组 取SubGroupWolfNumber+1行 到 最后一行
SubGroupsSecondPosition = Positions(SubGroupWolfNumber+1:SearchAgents_no,:);
l=0; %循环计数器

%进入循环
while l<Max_iter
    
    %计算第一小组
    %遍历所有狼的个数
    for i=1:size(SubGroupsFirstPosition,1)
        
        %返回超过搜索界限的代理
        Flag4ub  = SubGroupsFirstPosition(i,:) > ub;
        Flag4lb = SubGroupsFirstPosition(i,:)<lb;
        
        %如果两个之一为true则结果为true ~(Flag4ub+Flag4lb)为0 前半部分则为0 不参与计算
        %只有后半部分参与计算 如果两者都超过了界限 则重新复制为界限值 的加和 成为了中间
        %如果只有其中一个界限超过了 那么只给一个界限 这个等式写的很好 应该记住
        SubGroupsFirstPosition(i,:) = (SubGroupsFirstPosition(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;

        %计算适应度
        fitness = fobj(SubGroupsFirstPosition(i,:));
        
        %更新alhpa beta delta
        
        %这里实现的稍微有一点问题 应该是大于
        if fitness < FirstGroupAlpha_score
            FirstGroupAlpha_score = fitness;%更新值
            FirstGroupAlpha_pos = SubGroupsFirstPosition(i,:);  %更新位置                    
        end
        
        if fitness > FirstGroupAlpha_score && fitness < FirstGroupBeta_score
            FirstGroupBeta_score = fitness;
            FirstGroupBeta_pos = SubGroupsFirstPosition(i,:);
        end
        
        
        if fitness > FirstGroupAlpha_score && fitness > FirstGroupBeta_score && fitness < FirstGroupDelta_score
            FirstGroupDelta_score = fitness;
            FirstGroupDelta_pos = SubGroupsFirstPosition(i,:);
        end
    end    
     
    
       %计算第二小组
        a = 2-l*((2)/Max_iter); % a decreases linearly from 2 to 0
        %%计算第二个小组
     for i=1:size(SubGroupsSecondPosition,1)
        
        %返回超过搜索界限的代理
        Flag4ub  = SubGroupsSecondPosition(i,:) > ub;
        Flag4lb = SubGroupsSecondPosition(i,:)<lb;
        
        %如果两个之一为true则结果为true ~(Flag4ub+Flag4lb)为0 前半部分则为0 不参与计算
        %只有后半部分参与计算 如果两者都超过了界限 则重新复制为界限值 的加和 成为了中间
        %如果只有其中一个界限超过了 那么只给一个界限 这个等式写的很好 应该记住
        SubGroupsSecondPosition(i,:) = (SubGroupsSecondPosition(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;

        %计算适应度
        fitness = fobj(SubGroupsSecondPosition(i,:));
        
        %更新alhpa beta delta
        
        %这里实现的稍微有一点问题 应该是大于
        if fitness < SecondGroupAlpha_score
            SecondGroupAlpha_score = fitness;%更新值
            SecondGroupAlpha_pos = SubGroupsSecondPosition(i,:);  %更新位置                    
        end
        
        if fitness > SecondGroupAlpha_score && fitness < SecondGroupBeta_score
            SecondGroupBeta_score = fitness;
            SecondGroupBeta_pos = SubGroupsSecondPosition(i,:);
        end
               
        if fitness > SecondGroupAlpha_score && fitness > SecondGroupBeta_score && fitness < SecondGroupDelta_score
            SecondGroupDelta_score = fitness;
            SecondGroupDelta_pos = SubGroupsSecondPosition(i,:);
        end
        
     end     
       
        %更新位置搜索代理 包括omega
        for i=1:size(SubGroupsFirstPosition,1)
            
            for j=1:size(SubGroupsFirstPosition,2)
                
                r1 = rand();
                r2 = rand();
                
                A1 = 2*a*r1-a;
                C1 = 2*r2;
                D_alpha = abs(C1*FirstGroupAlpha_pos(j)-SubGroupsFirstPosition(i,j));
                X1 = FirstGroupAlpha_pos(j) - A1*D_alpha;
                
                r1 = rand();
                r2 = rand();
                
                A2 = 2*a*r1-a;
                C2 = 2*r2;
                
                D_beta = abs(C2*FirstGroupBeta_pos(j)-SubGroupsFirstPosition(i,j));
                X2 = FirstGroupBeta_pos(j)-A2*D_beta;
                
                %omega的值
                A3=2*a*r1-a; % Equation (3.3)
                C3=2*r2; % Equation (3.4)

                D_delta=abs(C3*FirstGroupDelta_pos(j)-SubGroupsFirstPosition(i,j)); % Equation (3.5)-part 3
                X3=FirstGroupDelta_pos(j)-A3*D_delta; % Equation (3.5)-part 3             

                SubGroupsFirstPosition(i,j)=(X1+X2+X3)/3;% Equation (3.7)  
            end
            
        end
        %第二个小组更新所有的狼
         %更新位置搜索代理 包括omega
        for i=1:size(SubGroupsSecondPosition,1)
            
            for j=1:size(SubGroupsSecondPosition,2)
                
                r1 = rand();
                r2 = rand();
                
                A1 = 2*a*r1-a;
                C1 = 2*r2;
                D_alpha = abs(C1*SecondGroupAlpha_pos(j)-SubGroupsSecondPosition(i,j));
                X1 = SecondGroupAlpha_pos(j) - A1*D_alpha;
                
                r1 = rand();
                r2 = rand();
                
                A2 = 2*a*r1-a;
                C2 = 2*r2;
                
                D_beta = abs(C2*SecondGroupBeta_pos(j)-SubGroupsSecondPosition(i,j));
                X2 = SecondGroupBeta_pos(j)-A2*D_beta;
                
                %omega的值
                A3=2*a*r1-a; % Equation (3.3)
                C3=2*r2; % Equation (3.4)

                D_delta=abs(C3*SecondGroupDelta_pos(j)-SubGroupsSecondPosition(i,j)); % Equation (3.5)-part 3
                X3=SecondGroupDelta_pos(j)-A3*D_delta; % Equation (3.5)-part 3             

                SubGroupsSecondPosition(i,j)=(X1+X2+X3)/3;% Equation (3.7)  
            end
            
        end
        
        %%分别对第一个小组 和 第二个小组 进行排序 按照适应度大小从大到小
        %%返回值代表的是2 3 4 6 77 1 依次是第几个狼的位置
        
        %%排序完毕 
        OrderedFirstGroupByFitness = orderFitness(SubGroupsFirstPosition,fobj);
        OrderedSecondGroupByFitness = orderFitness(SubGroupsSecondPosition,fobj);
        
        %%两组狼都已经结束了
        %Update the positions of members by(5)-(11) finish
       
        %Perform oppisition-based learning on subgroups
        
        %遍历两个狼群进行 oppsition  ！！！！！！为什么文献没有详细的说明 这个策略到底是什么策略 比例是多少 是怎么变化的
        
        %用总迭代次数 除以 本批次狼  可以计算second剩余的大小
        
        %初始比例是0.5一般一般             
        %改变转换因子
        OppositionScale = 0.5 + l*(0.5/Max_iter);

        %最好的狼改变的维数
        currentWolfchangeno = 1;
        
        %改变维数的增加 每增加一个狼其改变的程度越     
        %第一个小组进行opposition操作
        for index=1:SubGroupWolfNumber
            
            wolfno = OrderedFirstGroupByFitness(index,2);
             
            %前半部分  firsthalf 的狼  前半部分的狼越来越多
            if SubGroupWolfNumber * OppositionScale>=index
                
                %从当前维数中随机取出几维 其中改变的维数随迭代次数的变化而变化
                %从OrderedFirstGroupByFitness 中取出一个狼原来的编号            
                %取随机数 1-dim中选择 currentWolfchangeno 个随机数                                            
                %将要被转变的维数 
                changedimdata = randi(dim,currentWolfchangeno,1);
                
                for k = 1:currentWolfchangeno                    
                    changedimno = changedimdata(k,1);
                    %改变这一列
                    SubGroupsFirstPosition(wolfno,changedimno) = -SubGroupsFirstPosition(wolfno,changedimno);                
                end               
                %所有列都变换完了 下一次需要转换的维数就要加一个了
                currentWolfchangeno = currentWolfchangeno+1;      
            
            %以下是需要全部转换的狼
            else
                %所有列都转换就可以了
                SubGroupsFirstPosition(wolfno,:) = -SubGroupsFirstPosition(wolfno,:);                         
            end
        end
        
        currentWolfchangeno = 1;
        for index=1:(SearchAgents_no-SubGroupWolfNumber)
            
            wolfno = OrderedSecondGroupByFitness(index,2);
             
            %前半部分  firsthalf 的狼  前半部分的狼越来越多
            if (SearchAgents_no-SubGroupWolfNumber) * OppositionScale>=index
                
                %从当前维数中随机取出几维 其中改变的维数随迭代次数的变化而变化
                %从OrderedFirstGroupByFitness 中取出一个狼原来的编号            
                %取随机数 1-dim中选择 currentWolfchangeno 个随机数                                            
                %将要被转变的维数 
                changedimdata = randi(dim,currentWolfchangeno,1);
                
                for k = 1:currentWolfchangeno                    
                    changedimno = changedimdata(k,1);
                    %改变这一列
                    SubGroupsSecondPosition(wolfno,changedimno) = -SubGroupsSecondPosition(wolfno,changedimno);                
                end               
                %所有列都变换完了 下一次需要转换的维数就要加一个了
                currentWolfchangeno = currentWolfchangeno+1;      
            
            %以下是需要全部转换的狼
            else
                %所有列都转换就可以了
                SubGroupsSecondPosition(wolfno,:) = -SubGroupsSecondPosition(wolfno,:);                         
            end
        end  

        
        %--------到此为止所有狼的转换工作已经结束可-----%
        %开始进行交叉工作
        %重新计算适应度
        OrderedFirstGroupByFitness = orderFitness(SubGroupsFirstPosition,fobj);
        OrderedSecondGroupByFitness = orderFitness(SubGroupsSecondPosition,fobj);
        
        %分别把找到每组前三 和 最后 三个
        
       SubGroupsSecondPosition(OrderedSecondGroupByFitness(size(SubGroupsSecondPosition,1),2),:)  =  SubGroupsFirstPosition(OrderedFirstGroupByFitness(1,2),:);
       SubGroupsSecondPosition(OrderedSecondGroupByFitness(size(SubGroupsSecondPosition,1)-1,2),:)  =  SubGroupsFirstPosition(OrderedFirstGroupByFitness(2,2),:);
       SubGroupsSecondPosition(OrderedSecondGroupByFitness(size(SubGroupsSecondPosition,1)-2,2),:)  =  SubGroupsFirstPosition(OrderedFirstGroupByFitness(2,2),:);
       
       SubGroupsFirstPosition(OrderedFirstGroupByFitness(size(SubGroupsFirstPosition,1),2),:)  =  SubGroupsSecondPosition(OrderedSecondGroupByFitness(1,2),:);
       SubGroupsFirstPosition(OrderedFirstGroupByFitness(size(SubGroupsFirstPosition,1)-1,2),:)  =  SubGroupsSecondPosition(OrderedSecondGroupByFitness(2,2),:);
       SubGroupsFirstPosition(OrderedFirstGroupByFitness(size(SubGroupsFirstPosition,1)-2,2),:)  =  SubGroupsSecondPosition(OrderedSecondGroupByFitness(2,2),:);
       
       %交换完毕
         
       l = l+1;
       %比较较好的值
       
       if OrderedSecondGroupByFitness(1,1) > OrderedFirstGroupByFitness(1,1) 
      
            Convergence_curve(l) =  OrderedFirstGroupByFitness(1,1) ;
            Alpha_score = OrderedFirstGroupByFitness(1,1) ;
            Alpha_pos = OrderedFirstGroupByFitness(1,2) ;
       else
            Convergence_curve(l) =  OrderedSecondGroupByFitness(1,1) ;
             Alpha_score = OrderedSecondGroupByFitness(1,1) ;
            Alpha_pos = OrderedSecondGroupByFitness(1,2) ;
       end
      
end

end