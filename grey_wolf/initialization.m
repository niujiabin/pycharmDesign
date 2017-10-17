function [Positions] = initialization( SearchAgents_no,dim,ub,lb )

Boundary_no = size(ub,2);%��ȡ�������� ���ؾ��������

if Boundary_no == 1
    Positions = rand(SearchAgents_no,dim).*(ub-lb)+lb;
end

%һ�д���һ����
if Boundary_no > 1
    for i=1:dim
        ub_i = ub(i);
        lb_i = lb(i);
        Positions(:,i) = rand(SearchAgents_no,1).*(ub_i-lb_i)+lb_i;
    end
end

end
