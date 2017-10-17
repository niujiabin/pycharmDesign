
%����ֵ �����µ� Positions �ǵ����λ���Ѿ��ı��ˡ�
function [ Positions,Alpha_pos,Alpha_score ] = disBandW(Positions,fobj,SearchAgents_no,dim,Alpha_pos,Alpha_score,BestPosition,K)

    %��ʼ�� distance ���� ��¼   ���ž���
    distance = zeros(SearchAgents_no,2);

    %�������Ǻ����Ž�������   �����ǲ��ᴦ��
    for i = 1: SearchAgents_no
        
        %��i���ǵ�ŷ�Ͼ���
        for j = 1:dim 
        
            distance(i,1) = distance(i,1) +  (Positions(i,j) - Alpha_pos(1,j)).^2;
        
        end
        
        %��ȡŷ�Ͼ���
        distance(i,1) = sqrt(distance(i,1));
        
        %��i���ǵ����
        distance(i,2) = i;
    end    
    %����ŷ�Ͼ����������
    distance = sortrows(distance,1);   
    %��ȡǰK����
    %�ȿ�������    
    alfitpos  = Alpha_pos * unifrnd(-2,2);
    alfitness = fobj(alfitpos);
    
   % if alfitness < Alpha_score       
       %����Ӧ�ü�¼����ֵ��ԭλ�� �������������  ����λ�����С�ڲ�����
    
    Positions(BestPosition,:) =  alfitpos;
    Alpha_score = alfitness;
    
    
    % end
%     if als < Alpha_score
%         
%         Position
%     end
        
    for i=2:K+1       
        %��ȡ�ǵ�λ��
        tempP = distance(i,2);
        %��ȡ�µ���������������µĸ����ǵ�λ��       
        fitpre =  fobj(Positions(tempP,:));
        
        temP2 = Positions(tempP,:) + distance(i,1) *  unifrnd(-2,2);        
        %��ȡ��λ�õĸ���ֵ
        fitnow = fobj(temP2);
        
        %�ж�������ĺ�ĸ���ֵ �� ԭ�������Ա�             
        if fitnow < fitpre           
            %�����ԭ���ã��͸���
            Positions(i,:) = temP2;                     
        end
        %û��ԭ���ĺ� �Ͳ�������
    end
    
    
end

