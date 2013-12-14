clear all;
close all;
cd outputData

cpp_output
figure(1); plot(t_cpp,a_cpp,'linewidth',2);xlabel('Time (s)','fontsize',16); ylabel('Ground acceleration (a_g)','fontsize',16); grid minor;
set(gca,'fontsize',16);
for j=1:4
    figure(2);subplot(4,1,4-j+1);plot(t_cpp,u_cpp(j,:),'linewidth',2);%legend([num2str(j),'-th floor'],'location','EastOutside');
    if j==1
    title([num2str(j),'st floor'],'fontsize',16);
    else if j==2
            title([num2str(j),'nd floor'],'fontsize',16);
        else if j==3
            title([num2str(j),'rd floor'],'fontsize',16);
        else
            title([num2str(j),'th floor'],'fontsize',16);
            end
        end
    end
    set(gca,'fontsize',16);
    if j==1
        xlabel('Time (s)','fontsize',16);
    end
    if j==2
        ylabel('Displacement (m)','fontsize',16);
    end
    figure(3);subplot(4,1,4-j+1);plot(t_cpp,ud_cpp(j,:),'linewidth',2);legend(strcat('dof ',num2str(j)),'location','southeast');
    set(gca,'fontsize',16);
    if j==1
        xlabel('Time (s)','fontsize',16);
    end
    if j==2
        ylabel('Velocity (m/s)','fontsize',16);
    end
    figure(4);subplot(4,1,4-j+1);plot(t_cpp,udd_cpp(j,:),'linewidth',2);
    %legend([num2str(j),'-th floor'],'location','EastOutside');
    if j==1
    title([num2str(j),'st floor'],'fontsize',16);
    else if j==2
            title([num2str(j),'nd floor'],'fontsize',16);
        else if j==3
            title([num2str(j),'rd floor'],'fontsize',16);
        else
            title([num2str(j),'th floor'],'fontsize',16);
            end
        end
    end
    grid minor;
    set(gca,'fontsize',16);
    if j==1
        xlabel('Time (s)','fontsize',16);
    end
    if j==2
        ylabel('Acceleration (m/s^2)','fontsize',16);
    end
    figure(5);subplot(4,1,4-j+1);plot(t_cpp,ru_cpp(j,:),'linewidth',2);legend(strcat('dof ',num2str(j-1),'-',num2str(j)));
    set(gca,'fontsize',16);
    if j==1
        xlabel('Time (s)','fontsize',16);
    end
    if j==2
        ylabel('Relative displacement(m)','fontsize',16);
    end
    figure(6+j-1);plot(ru_cpp(j,:),resfor_cpp(j,:),'linewidth',2); xlabel('Relative displacement (m)','fontsize',16); ylabel('Restoring force (N)','fontsize',16); grid on;
     set(gca,'fontsize',16);
     legend(strcat('floor ',num2str(j)),'location','southeast');
end
for j=1:9
figure(j);print(gcf,'-dpng',strcat('hysteretic_eg_cpp_',num2str(j)))
end
cd ..
