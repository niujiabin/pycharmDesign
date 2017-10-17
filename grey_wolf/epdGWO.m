
function [Alpha_score,Alpha_pos,Convergence_curve] = GWO(SearchAgents_no,Max_iter,lb,ub,dim,fobj)

%��ʼ�� alpha beta and delta_pos
Alpha_pos = zeros(1,dim);%����1��dim�е�alpha�ǵ�λ��
Alpha_score = inf;%change this to -inf for maximization problems

Beta_pos = zeros(1,dim);
Beta_score = inf;%change this to -inf for maximinzation problems

Delta_pos = zeros(1,dim);
Delta_score = inf;

%��ʼ������λ��

Positions = initialization(SearchAgents_no,dim,ub,lb);

Convergence_curve = zeros(1,Max_iter);

l=0; %ѭ��������

%����ѭ��

while l<Max_iter
    
    %���������ǵĸ���
    for i=1:size(Positions,1)
        
        %���س����������޵Ĵ���
        Flag4ub  = Positions(i,:) > ub;
        Flag4lb = Positions(i,:)<lb;
        
        %�������֮һΪtrue����Ϊtrue ~(Flag4ub+Flag4lb)Ϊ0 ǰ�벿����Ϊ0 ���������
        %ֻ�к�벿�ֲ������ ������߶������˽��� �����¸���Ϊ����ֵ �ļӺ� ��Ϊ���м�
        %���ֻ������һ�����޳����� ��ôֻ��һ������ �����ʽд�ĺܺ� Ӧ�ü�ס
        Positions(i,:) = (Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;

        %������Ӧ��
        fitness = fobj(Positions(i,:));
        
       
      
        
        %����alhpa beta delta
        
        %����ʵ�ֵ���΢��һ������ Ӧ���Ǵ���
        if fitness < Alpha_score
            Alpha_score = fitness;%����ֵ
            Alpha_pos = Positions(i,:);  %����λ��                    
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
        a = 2-l*((2)/Max_iter); % a decreases linearly from 2 to 0        
        %����λ���������� ����omega
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
                              
                %omega��ֵ
                A3=2*a*r1-a; % Equation (3.3)
                C3=2*r2; % Equation (3.4)

                D_delta=abs(C3*Delta_pos(j)-Positions(i,j)); % Equation (3.5)-part 3
                X3=Delta_pos(j)-A3*D_delta; % Equation (3.5)-part 3             
          
                Positions(i,j)=(X1+X2+X3)/3;% Equation (3.7) 
                
            end
            
            
        end
        
        l = l+1;
        
        %������Ӧ��ֵ ��������
        
        orderFitnessPosition = orderFitness(Positions,fobj);
        
        %�������Ⱥ��Ϣ���и���
        
        f = rand();
        
        if f >=0.8
       
             %��һ����ǽ�����Ϣ����
        for i=SearchAgents_no/2:SearchAgents_no
            
            r = randperm(4,1);
            
            if r==1
            
                 %���������
                 %�Ǽӻ��Ǽ�
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
                 
            end
            
            %�ȸ��ʵ��������1���ĵ�����
        end
        end
            
        
       
        
        Convergence_curve(l) = Alpha_score;
end

