function [GWO_cg_curve,GWO_cg_curve1,GWO_cg_curve2] = main(SearchAgents_no,Function_name,Max_iteration)
%%向量C代表障碍 向量A的取值代表是否靠近礼物 还是远离猎物
%SearchAgents_no = 50;
%14 大多数得不到最好解  F13用狼群算法有时候得不到最好解 F12 PSO更容易得到最好解 GWO不太好
%F10 F11 可以提高收敛速度 比PSO要好 F9可以考虑提高收敛速度 比PSO效果要好
%F6需要改进
%Max_iteration = 10;
K = 10;
[lb,ub,dim,fobj] =  Get_Function_details(Function_name);
[Best_score,Best_pos,GWO_cg_curve] = GWO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
%%对于高维复杂极值 粒子群算法还是很好的
PSO_cg_curve = PSO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
%低位的时候有些过度迭代
[Best_score1,Best_pos1,GWO_cg_curve1] =GWO1(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
%F12 这个函数并没有那么优秀 极值点比较多 PSO算法的优势体现出来了
%这个方法收敛太快了 达到局部极值点了  可以引进多个指导狼  动态的加入指导狼 不一定非得是三个 根据当前的情况 每次基本不到50代就停止了
%收敛速度太快 有点不太好
[Best_score2,Best_pos2,GWO_cg_curve2] =wdGWO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
%这个方法在某些情况下是好的 总容易跳跃 不是很好呀 但是最后跳跃跳跃 能找到那个较好的值 有时候却找不到
%[Best_score3,Best_pos3,GWO_cg_curve3] =paGWO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
%加入狼群算法的思想 
[Best_score4,Best_pos4,GWO_cg_curve4] =cGWO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
[Best_score5,Best_pos5,GWO_cg_curve5] =mGWO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj,K);
[Best_score6,Best_pos6,GWO_cg_curve6] =epdGWO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
%指定窗口大小
figure('Position',[500 500 660 290]);
%画图

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
semilogy(GWO_cg_curve2,'Color','y')%wd的这个就是瞎写的 很差很差 实验结果一下子就看出来了 他里面的论文全是0 根本不可能
%hold on
%semilogy(GWO_cg_curve3,'Color','m')%wd的这个就是瞎写的 很差很差 实验结果一下子就看出来了 他里面的论文全是0 根本不可能
hold on
semilogy(GWO_cg_curve4,'Color','c')
hold on
semilogy(GWO_cg_curve5,'Color','m')%wd的这个就是瞎写的 很差很差 实验结果一下子就看出来了 他里面的论文全是0 根本不可能

hold on
semilogy(GWO_cg_curve6,'Color','k')%wd的这个就是瞎写的 很差很差 实验结果一下子就看出来了 他里面的论文全是0 根本不可能

title('Objective space')
xlabel('Iteration');
ylabel('Best score obtained so far');

axis tight
grid on
box on
legend('GWO','PSO','GWO1','wdGWO','cGWO','mGWO','epdGWO')

display(['The best solution obtained by GWO is : ', num2str(Best_pos)]);
display(['The best optimal value of the objective funciton found by GWO is : ', num2str(Best_score)]);

        



