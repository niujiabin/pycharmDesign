function [ cg_curve ] = PSO( N,Max_iteration,lb,ub,dim,fobj )
%PSO �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

%Particual Swarm Optimization
%����Ⱥ�㷨�������Ϣ

%������ÿһλ�ռ��ϵ�ÿһά����һ������ٶ�Vdmax�����������ӵ��ٶȽ�������
Vmax = 6;%����ٶ�
noP = N;%���ӵ�����
wMax = 0.9;
wMin = 0.2;

%��ʼ�� ���ʵ�c1 c2���Լӿ������ٶ� �ֲ���������ֲ�����

%ѧϰ���� ����Ƽ���ϵ�� 

%�ֱ����pBest��Gbest������е���󲽳�
%�������Ӹ��徭���Ⱥ�徭����������й켣��Ӱ�죬��Ӧ����Ⱥ֮�����Ϣ����

%���c1=0 ����ֻ��Ⱥ�徭�� �����ٶȿ� ����������ֲ�����
%c2=0 ����û��Ⱥ����Ϣ ���Ӹ������� �õ���ø��ʱȽ�С ���һ������c1=c2
%������������ͬ��Ӱ���� �ϵ͵�c1 c2�ǵ������ǻ���Զ��Ŀ������� �ϸߵ�c1��c2
%ֵ�������˶���Խ��Ŀ������
c1=2;
c2=2;



iter = Max_iteration;
vel = zeros(noP,dim);
pBestScore = zeros(noP);%Ĭ����һ��
pBest = zeros(noP,dim);%һ�����ӱ������ҵ�����ý⣬�����弫ֵ����һ����ֵ����������Ⱥ����������������
%���������ﵽ�����Ž��ίȫ�ּ�ֵ��
%�ҵ�������ý�󣬽�������Pso������Ҫ�ļ��ٹ��� ÿ�����Ӳ��ϸı����ڽ�ռ��е��ٶȣ�
%�Ծ����ܵس�pBest��gBest��ָ��������ȥ��
gBest = zeros(1,dim);%��ǰ��ý�
cg_curve = zeros(1,iter);

%�����ʼ��agents
pos = initialization(noP,dim,ub,lb);


for i=1:noP
    pBestScore(i) = inf;%Ĭ�������ֵ
end

gBestScore = inf;


for l=1:iter
    
    %�Ƿ񳬽��޵��ж�
    Flag4ub = pos(i,:)>ub;
    Flag4lb = pos(i,:)<lb;
    
    pos(i,:) = (pos(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
    
    for i=1:size(pos,1)
        %����Ŀ�꺯��ֵ ����ÿһ������
        fitness = fobj(pos(i,:));
        
        %ÿһ�������Ž�
        if (pBestScore(i) > fitness)
            pBestScore(i) = fitness;
            pBest(i,:) = pos(i,:);
        end
        
        %�ж�ȫ�����Ž�
        if (gBestScore > fitness)
            gBestScore = fitness;
            gBest = pos(i,:);
        end
    end
    
    %����PSO��wȨֵ ���ŵ��������ĸı���ı�
    w = wMax-l*((wMax-wMin)/iter);
    
    %�����ٶȺ����ӵ�λ��
    for i=1:size(pos,1)
        for j=1:size(pos,2)
            vel(i,j) = w*vel(i,j)+c1*rand()*(pBest(i,j)-pos(i,j))+rand()*c2*(gBest(j)-pos(i,j));
            
            %���ܳ�������ٶ�
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

