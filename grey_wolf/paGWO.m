
%%�����㷨�Ĳ����㷨   2017-03-21 ʵ�ַ���
%%Aparallel grey wolf optimizer 
function [Alpha_score,Alpha_pos,Convergence_curve] = paGWO(SearchAgents_no,Max_iter,lb,ub,dim,fobj)

%��ʼ�� alpha beta and delta_pos
%��һ��С��
FirstGroupAlpha_pos = zeros(1,dim);%����1��dim�е�alpha�ǵ�λ��
FirstGroupAlpha_score = inf;%change this to -inf for maximization problems

FirstGroupBeta_pos = zeros(1,dim);
FirstGroupBeta_score = inf;%change this to -inf for maximinzation problems

FirstGroupDelta_pos = zeros(1,dim);
FirstGroupDelta_score = inf;

%�ڶ���С��
SecondGroupAlpha_pos = zeros(1,dim);%����1��dim�е�alpha�ǵ�λ��
SecondGroupAlpha_score = inf;%change this to -inf for maximization problems

SecondGroupBeta_pos = zeros(1,dim);
SecondGroupBeta_score = inf;%change this to -inf for maximinzation problems

SecondGroupDelta_pos = zeros(1,dim);
SecondGroupDelta_score = inf;

%��ʼ������λ��

Positions = initialization(SearchAgents_no,dim,ub,lb);

Convergence_curve = zeros(1,Max_iter);



%Generate initial population and divide them into 2 subgroups
%��һ��С�������
SubGroupWolfNumber = SearchAgents_no/2;
%��һ��С�� ȡǰSubGroupWolfNumber��
SubGroupsFirstPosition = Positions(1:SubGroupWolfNumber,:);
%�ڶ���С�� ȡSubGroupWolfNumber+1�� �� ���һ��
SubGroupsSecondPosition = Positions(SubGroupWolfNumber+1:SearchAgents_no,:);
l=0; %ѭ��������

%����ѭ��
while l<Max_iter
    
    %�����һС��
    %���������ǵĸ���
    for i=1:size(SubGroupsFirstPosition,1)
        
        %���س����������޵Ĵ���
        Flag4ub  = SubGroupsFirstPosition(i,:) > ub;
        Flag4lb = SubGroupsFirstPosition(i,:)<lb;
        
        %�������֮һΪtrue����Ϊtrue ~(Flag4ub+Flag4lb)Ϊ0 ǰ�벿����Ϊ0 ���������
        %ֻ�к�벿�ֲ������ ������߶������˽��� �����¸���Ϊ����ֵ �ļӺ� ��Ϊ���м�
        %���ֻ������һ�����޳����� ��ôֻ��һ������ �����ʽд�ĺܺ� Ӧ�ü�ס
        SubGroupsFirstPosition(i,:) = (SubGroupsFirstPosition(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;

        %������Ӧ��
        fitness = fobj(SubGroupsFirstPosition(i,:));
        
        %����alhpa beta delta
        
        %����ʵ�ֵ���΢��һ������ Ӧ���Ǵ���
        if fitness < FirstGroupAlpha_score
            FirstGroupAlpha_score = fitness;%����ֵ
            FirstGroupAlpha_pos = SubGroupsFirstPosition(i,:);  %����λ��                    
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
     
    
       %����ڶ�С��
        a = 2-l*((2)/Max_iter); % a decreases linearly from 2 to 0
        %%����ڶ���С��
     for i=1:size(SubGroupsSecondPosition,1)
        
        %���س����������޵Ĵ���
        Flag4ub  = SubGroupsSecondPosition(i,:) > ub;
        Flag4lb = SubGroupsSecondPosition(i,:)<lb;
        
        %�������֮һΪtrue����Ϊtrue ~(Flag4ub+Flag4lb)Ϊ0 ǰ�벿����Ϊ0 ���������
        %ֻ�к�벿�ֲ������ ������߶������˽��� �����¸���Ϊ����ֵ �ļӺ� ��Ϊ���м�
        %���ֻ������һ�����޳����� ��ôֻ��һ������ �����ʽд�ĺܺ� Ӧ�ü�ס
        SubGroupsSecondPosition(i,:) = (SubGroupsSecondPosition(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;

        %������Ӧ��
        fitness = fobj(SubGroupsSecondPosition(i,:));
        
        %����alhpa beta delta
        
        %����ʵ�ֵ���΢��һ������ Ӧ���Ǵ���
        if fitness < SecondGroupAlpha_score
            SecondGroupAlpha_score = fitness;%����ֵ
            SecondGroupAlpha_pos = SubGroupsSecondPosition(i,:);  %����λ��                    
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
       
        %����λ���������� ����omega
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
                
                %omega��ֵ
                A3=2*a*r1-a; % Equation (3.3)
                C3=2*r2; % Equation (3.4)

                D_delta=abs(C3*FirstGroupDelta_pos(j)-SubGroupsFirstPosition(i,j)); % Equation (3.5)-part 3
                X3=FirstGroupDelta_pos(j)-A3*D_delta; % Equation (3.5)-part 3             

                SubGroupsFirstPosition(i,j)=(X1+X2+X3)/3;% Equation (3.7)  
            end
            
        end
        %�ڶ���С��������е���
         %����λ���������� ����omega
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
                
                %omega��ֵ
                A3=2*a*r1-a; % Equation (3.3)
                C3=2*r2; % Equation (3.4)

                D_delta=abs(C3*SecondGroupDelta_pos(j)-SubGroupsSecondPosition(i,j)); % Equation (3.5)-part 3
                X3=SecondGroupDelta_pos(j)-A3*D_delta; % Equation (3.5)-part 3             

                SubGroupsSecondPosition(i,j)=(X1+X2+X3)/3;% Equation (3.7)  
            end
            
        end
        
        %%�ֱ�Ե�һ��С�� �� �ڶ���С�� �������� ������Ӧ�ȴ�С�Ӵ�С
        %%����ֵ�������2 3 4 6 77 1 �����ǵڼ����ǵ�λ��
        
        %%������� 
        OrderedFirstGroupByFitness = orderFitness(SubGroupsFirstPosition,fobj);
        OrderedSecondGroupByFitness = orderFitness(SubGroupsSecondPosition,fobj);
        
        %%�����Ƕ��Ѿ�������
        %Update the positions of members by(5)-(11) finish
       
        %Perform oppisition-based learning on subgroups
        
        %����������Ⱥ���� oppsition  ������������Ϊʲô����û����ϸ��˵�� ������Ե�����ʲô���� �����Ƕ��� ����ô�仯��
        
        %���ܵ������� ���� ��������  ���Լ���secondʣ��Ĵ�С
        
        %��ʼ������0.5һ��һ��             
        %�ı�ת������
        OppositionScale = 0.5 + l*(0.5/Max_iter);

        %��õ��Ǹı��ά��
        currentWolfchangeno = 1;
        
        %�ı�ά�������� ÿ����һ������ı�ĳ̶�Խ     
        %��һ��С�����opposition����
        for index=1:SubGroupWolfNumber
            
            wolfno = OrderedFirstGroupByFitness(index,2);
             
            %ǰ�벿��  firsthalf ����  ǰ�벿�ֵ���Խ��Խ��
            if SubGroupWolfNumber * OppositionScale>=index
                
                %�ӵ�ǰά�������ȡ����ά ���иı��ά������������ı仯���仯
                %��OrderedFirstGroupByFitness ��ȡ��һ����ԭ���ı��            
                %ȡ����� 1-dim��ѡ�� currentWolfchangeno �������                                            
                %��Ҫ��ת���ά�� 
                changedimdata = randi(dim,currentWolfchangeno,1);
                
                for k = 1:currentWolfchangeno                    
                    changedimno = changedimdata(k,1);
                    %�ı���һ��
                    SubGroupsFirstPosition(wolfno,changedimno) = -SubGroupsFirstPosition(wolfno,changedimno);                
                end               
                %�����ж��任���� ��һ����Ҫת����ά����Ҫ��һ����
                currentWolfchangeno = currentWolfchangeno+1;      
            
            %��������Ҫȫ��ת������
            else
                %�����ж�ת���Ϳ�����
                SubGroupsFirstPosition(wolfno,:) = -SubGroupsFirstPosition(wolfno,:);                         
            end
        end
        
        currentWolfchangeno = 1;
        for index=1:(SearchAgents_no-SubGroupWolfNumber)
            
            wolfno = OrderedSecondGroupByFitness(index,2);
             
            %ǰ�벿��  firsthalf ����  ǰ�벿�ֵ���Խ��Խ��
            if (SearchAgents_no-SubGroupWolfNumber) * OppositionScale>=index
                
                %�ӵ�ǰά�������ȡ����ά ���иı��ά������������ı仯���仯
                %��OrderedFirstGroupByFitness ��ȡ��һ����ԭ���ı��            
                %ȡ����� 1-dim��ѡ�� currentWolfchangeno �������                                            
                %��Ҫ��ת���ά�� 
                changedimdata = randi(dim,currentWolfchangeno,1);
                
                for k = 1:currentWolfchangeno                    
                    changedimno = changedimdata(k,1);
                    %�ı���һ��
                    SubGroupsSecondPosition(wolfno,changedimno) = -SubGroupsSecondPosition(wolfno,changedimno);                
                end               
                %�����ж��任���� ��һ����Ҫת����ά����Ҫ��һ����
                currentWolfchangeno = currentWolfchangeno+1;      
            
            %��������Ҫȫ��ת������
            else
                %�����ж�ת���Ϳ�����
                SubGroupsSecondPosition(wolfno,:) = -SubGroupsSecondPosition(wolfno,:);                         
            end
        end  

        
        %--------����Ϊֹ�����ǵ�ת�������Ѿ�������-----%
        %��ʼ���н��湤��
        %���¼�����Ӧ��
        OrderedFirstGroupByFitness = orderFitness(SubGroupsFirstPosition,fobj);
        OrderedSecondGroupByFitness = orderFitness(SubGroupsSecondPosition,fobj);
        
        %�ֱ���ҵ�ÿ��ǰ�� �� ��� ����
        
       SubGroupsSecondPosition(OrderedSecondGroupByFitness(size(SubGroupsSecondPosition,1),2),:)  =  SubGroupsFirstPosition(OrderedFirstGroupByFitness(1,2),:);
       SubGroupsSecondPosition(OrderedSecondGroupByFitness(size(SubGroupsSecondPosition,1)-1,2),:)  =  SubGroupsFirstPosition(OrderedFirstGroupByFitness(2,2),:);
       SubGroupsSecondPosition(OrderedSecondGroupByFitness(size(SubGroupsSecondPosition,1)-2,2),:)  =  SubGroupsFirstPosition(OrderedFirstGroupByFitness(2,2),:);
       
       SubGroupsFirstPosition(OrderedFirstGroupByFitness(size(SubGroupsFirstPosition,1),2),:)  =  SubGroupsSecondPosition(OrderedSecondGroupByFitness(1,2),:);
       SubGroupsFirstPosition(OrderedFirstGroupByFitness(size(SubGroupsFirstPosition,1)-1,2),:)  =  SubGroupsSecondPosition(OrderedSecondGroupByFitness(2,2),:);
       SubGroupsFirstPosition(OrderedFirstGroupByFitness(size(SubGroupsFirstPosition,1)-2,2),:)  =  SubGroupsSecondPosition(OrderedSecondGroupByFitness(2,2),:);
       
       %�������
         
       l = l+1;
       %�ȽϽϺõ�ֵ
       
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