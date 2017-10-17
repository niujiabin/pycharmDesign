%�����㷨����������û�н���

%�Ľ�˼· ���һ������ �ܵĽϿ� ����Χ��������Ϊ  ��Щ���ŵ��ǲ����ϴ�
function [Alpha_score,Alpha_pos,Convergence_curve] = cGWO(SearchAgents_no,Max_iter,lb,ub,dim,fobj)
%��ǰȫ�����Ž�

gBest = zeros(1,dim);
gBestScore = inf; 

BestPosition = 0;
pBestScore = zeros(SearchAgents_no);%Ĭ����һ��
pBest = zeros(SearchAgents_no,dim);

Vmax = 6;%����ٶ�
wMax = 0.9;
wMin = 0.2;

%��ʼ�� alpha beta and delta_pos
Alpha_pos = zeros(1,dim);%����1��dim�е�alpha�ǵ�λ��
Alpha_score = inf;%change this to -inf for maximization problems

Beta_pos = zeros(1,dim);
Beta_score = inf;%change this to -inf for maximinzation problems

Delta_pos = zeros(1,dim);
Delta_score = inf;

%��ʼ������λ��
vel = zeros(SearchAgents_no,dim);

Positions = initialization(SearchAgents_no,dim,ub,lb);

Convergence_curve = zeros(1,Max_iter);

l=0; %ѭ��������

%����ѭ��

twpw2 = zeros(1,dim);

for i=1:SearchAgents_no
    pBestScore(i) = inf;%Ĭ�������ֵ
end

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
        
%          twpw = Positions(i,:);
        
%         %��ÿһ��ִ��p���������
%         for p=1:3
%             
%             for j=1:size(Positions,2);
%                 
%                  %8���������� 2�ǲ���
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
        %������ÿһ���ǵ�ʱ��
        %���ȡһ���ǵ��෴ֵ
%         wolfs =  randi(SearchAgents_no,1,1);
%         ��һ
%         wolfa = wolfs(1,1);
% 
%         ��ȡ�������ǵ���Ӧ��
%         fia = fobj(Positions(wolfa,:));
%         
%         tempPosition =  (Positions(wolfa,:) + Positions(i,:))/2;
%         ��ȡ��Ӧ��ֵ
%         fitnesstemp = fobj(tempPosition);        
%         
%         if fitnesstemp < fia && fitnesstemp < fitness          
%             fitness = fitnesstemp; 
%         end
        %����alhpa beta delta
        
        %����ʵ�ֵ���΢��һ������ Ӧ���Ǵ���
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
   %  x=0:0.1:600; 
%y=sigmf(x,[1/60 300]); 
    %�Ľ�1 �Բ������иĽ�
    
    %����sigmoid����
       a =  -(1/(1+(exp(-(10/Max_iter)*(l-Max_iter/2)))))*2+2;
      % a = 2-2*(1/(exp(1)-1)*(exp(l/Max_iter)-1)); % a decreases linearly from 2 to 0
%        a = 2-l*((2)/Max_iter); % a decreases linearly from 2 to 0    
         w = wMax-l*((wMax-wMin)/Max_iter);
        %����λ���������� ����omega
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
                
                %omega��ֵ
                A3=2*a*r1-a; % Equation (3.3)
                C3=2*r2; % Equation (3.4)

                D_delta=C3*Delta_pos(j)-Positions(i,j) + 2*r1*(2-a)*1.5*(pBest(i,j)-Positions(i,j)); % Equation (3.5)-part 3
                X3= vel(i,j) + Delta_pos(j)-A3*D_delta; % Equation (3.5)-part 3             
             
                r1 = rand();
                r2 = rand();
                
                %ȫ�����Ž���֪��Alpha_score���ټ�һ������ǰ�����Ž⣬
                A4=2*a*r1-a; % Equation (3.3)
                C4=2*r2; % Equation (3.4)
      
                P_delta=C4*pBest(i,j)-Positions(i,j) + 2*r1*1.7*(pBest(i,j)-Positions(i,j)); % Equation (3.5)-part 3
                X4= vel(i,j) + pBest(i,j)-A4*P_delta; % Equation (3.5)-part 3 
                
                %����һ
                %ȫ�����Ž���֪��Alpha_score���ټ�һ������ǰ�����Ž⣬
                A5=2*a*r1-a; % Equation (3.3)
                C5=2*r2; % Equation (3.4)
              
                %����Ҫ������ͷͷ�˵� ��Ҫ���ǵ�ǰ�Ǳ���
                Positions(i,j)=(X1+X2+X3+X4)/4 ;% Equation (3.7)
                  
                 
%                 Positions(i,j) =  rand()*abs(Positions(i,j)/(vel(i,j) + Positions(i,j))) * Positions(i,j) + rand()*abs(vel(i,j)/(vel(i,j) + Positions(i,j))) * vel(i,j);
                         
            end           
        end       
        
        %��ÿһ�������������н��ٶȵĸ���
       % [ Positions , Alpha_pos , Alpha_score ] = disBandW(Positions,fobj,SearchAgents_no,dim,Alpha_pos,Alpha_score,BestPosition,3);
        
         
        orderFitnessPosition = orderFitness(Positions,fobj);
        
        %��һ����ǽ�����Ϣ����
        for i=fix(SearchAgents_no/2):SearchAgents_no
            
            r = randperm(5,1);
            
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
             elseif r==5
                       
                 si = rand();
                 if si>=0.5
                     Positions(orderFitnessPosition(i,2),:) = pBest(i,:) + (ub-lb*rand()+lb);     
                 else
                      Positions(orderFitnessPosition(i,2),:) = pBest(i,:) - (ub-lb*rand()+lb);
                 end
                 
            end
            
            %�ȸ��ʵ��������1���ĵ�����
        end
        l = l+1;
        Convergence_curve(l) = gBestScore;
end

