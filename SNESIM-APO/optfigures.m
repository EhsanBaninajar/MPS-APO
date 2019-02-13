%this part show all the results as a plot
%% figure 1 Theta
figure(1);
subplot(3,1,1)
xlabel('Iterations','FontSize',10)
ylabel('Max Cond')
hold on
    title(['Maximum number of conditioning data=',num2str(round(theta(1)*100,0))])
    plot(k,round(theta(1)*100,0),'or'),ylim([0 100]),xlim([0 inf])

subplot(3,1,2)
xlabel('Iterations','FontSize',10)
ylabel('Cmin')
hold on
    title(['minimum number of replicates=',num2str(round(theta(2)*100,0))])
    plot(k,round(theta(2)*100,0),'or'),ylim([0 100]),xlim([0 inf])

subplot(3,1,3)
xlabel('Iterations','FontSize',10)
ylabel('Nb Multigrids')
hold on
    title(['number of multiple grids=',num2str(round(theta(3)*10,0))])
    plot(k,round(theta(3)*10,0),'or'),ylim([0 10]),xlim([0 inf])
%% figure 2 MSE function
Current=[k,MSE_globalmain];
figure(2);
% title('MSE Function value')
xlabel('Iterations')
ylabel('MSE')
hold on
xlim([0 inf]),ylim([0 inf])
if stage==0
    scatter(k,MSE_globalmain,50,'MarkerEdgeColor','b','MarkerFaceColor',[0 0.5 0.5])
else
    if MSE_globalmain<=(Ttradeoff+100)/100*funcbeststage0
        scatter(k,MSE_globalmain,50,'MarkerEdgeColor','b','MarkerFaceColor',[0 0.5 0.5])
    else
        scatter(k,MSE_globalmain,50,'MarkerEdgeColor','b','MarkerFaceColor','r')
    end
end
% if k==0
%     Best=Current;
% end
% if Current(1,2)<Best(1,2)
%     line([Current(1),Best(1)],[Current(2),Best(2)])
%     Best=Current;
% end
if counter==5 && stage==0 || stage==1
    line([k,k+1],[(Ttradeoff+100)/100*funcbeststage0,(Ttradeoff+100)/100*funcbeststage0],'color','red')
end
%% figure 3 Time
figure(3);
% title('Time')
xlabel('Iterations')
ylabel('Computation Time')
hold on
xlim([0 inf]),ylim([0 inf])
scatter(k,Timemain,50,'MarkerEdgeColor','b','MarkerFaceColor',[0 0.5 0.5]);