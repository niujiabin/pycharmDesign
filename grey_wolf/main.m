function [GWO_cg_curve,GWO_cg_curve1,GWO_cg_curve2] = main(SearchAgents_no,Function_name,Max_iteration)
%%����C�����ϰ� ����A��ȡֵ�����Ƿ񿿽����� ����Զ������
%SearchAgents_no = 50;
%14 ������ò�����ý�  F13����Ⱥ�㷨��ʱ��ò�����ý� F12 PSO�����׵õ���ý� GWO��̫��
%F10 F11 ������������ٶ� ��PSOҪ�� F9���Կ�����������ٶ� ��PSOЧ��Ҫ��
%F6��Ҫ�Ľ�
%Max_iteration = 10;
K = 10;
[lb,ub,dim,fobj] =  Get_Function_details(Function_name);
[Best_score,Best_pos,GWO_cg_curve] = GWO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
%%���ڸ�ά���Ӽ�ֵ ����Ⱥ�㷨���Ǻܺõ�
PSO_cg_curve = PSO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
%��λ��ʱ����Щ���ȵ���
[Best_score1,Best_pos1,GWO_cg_curve1] =GWO1(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
%F12 ���������û����ô���� ��ֵ��Ƚ϶� PSO�㷨���������ֳ�����
%�����������̫���� �ﵽ�ֲ���ֵ����  �����������ָ����  ��̬�ļ���ָ���� ��һ���ǵ������� ���ݵ�ǰ����� ÿ�λ�������50����ֹͣ��
%�����ٶ�̫�� �е㲻̫��
[Best_score2,Best_pos2,GWO_cg_curve2] =wdGWO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
%���������ĳЩ������Ǻõ� ��������Ծ ���Ǻܺ�ѽ ���������Ծ��Ծ ���ҵ��Ǹ��Ϻõ�ֵ ��ʱ��ȴ�Ҳ���
%[Best_score3,Best_pos3,GWO_cg_curve3] =paGWO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
%������Ⱥ�㷨��˼�� 
[Best_score4,Best_pos4,GWO_cg_curve4] =cGWO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
[Best_score5,Best_pos5,GWO_cg_curve5] =mGWO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj,K);
[Best_score6,Best_pos6,GWO_cg_curve6] =epdGWO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
%ָ�����ڴ�С
figure('Position',[500 500 660 290]);
%��ͼ

subplot(1,2,1);
func_plot(Function_name);

title('Parameter space');
xlabel('x_1')
ylabel('x_2');
zlabel([Function_name,'(x_1,x_2)']);


%Draw objective space
subplot(1,2,2);
semilogy(GWO_cg_curve,'Color','r')
hold on
semilogy(PSO_cg_curve,'Color','b')
hold on
semilogy(GWO_cg_curve1,'Color','g')
hold on
semilogy(GWO_cg_curve2,'Color','y')%wd���������Ϲд�� �ܲ�ܲ� ʵ����һ���ӾͿ������� �����������ȫ��0 ����������
%hold on
%semilogy(GWO_cg_curve3,'Color','m')%wd���������Ϲд�� �ܲ�ܲ� ʵ����һ���ӾͿ������� �����������ȫ��0 ����������
hold on
semilogy(GWO_cg_curve4,'Color','c')
hold on
semilogy(GWO_cg_curve5,'Color','m')%wd���������Ϲд�� �ܲ�ܲ� ʵ����һ���ӾͿ������� �����������ȫ��0 ����������

hold on
semilogy(GWO_cg_curve6,'Color','k')%wd���������Ϲд�� �ܲ�ܲ� ʵ����һ���ӾͿ������� �����������ȫ��0 ����������

title('Objective space')
xlabel('Iteration');
ylabel('Best score obtained so far');

axis tight
grid on
box on
legend('GWO','PSO','GWO1','wdGWO','cGWO','mGWO','epdGWO')

display(['The best solution obtained by GWO is : ', num2str(Best_pos)]);
display(['The best optimal value of the objective funciton found by GWO is : ', num2str(Best_score)]);

        



