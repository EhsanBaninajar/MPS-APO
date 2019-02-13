%this part show all the results as a plot
%% figure 1 Theta
figure(1);
subplot(2,2,1)
xlabel('Iterations','FontSize',10)
ylabel('tilesize')
hold on
    title(['Optimal size for the patch size=',num2str(round(theta(1)*100,0))])
    plot(k,round(theta(1)*100,0),'or'),ylim([0 100]),xlim([0 inf])

subplot(2,2,2)
xlabel('Iterations','FontSize',10)
ylabel('overlap')
hold on
    title(['Optimal size for the overlap region=',num2str(round(theta(2)*100*theta(1),0))])
    plot(k,(round(theta(2)*100*theta(1),0)),'or'),ylim([0 50]),xlim([0 inf])

subplot(2,2,3)
xlabel('Iterations','FontSize',10)
ylabel('Number of replicas')
hold on
    title(['Standard Number of replicas=',num2str(round(theta(3)*10,0))])
    plot(k,round(theta(3)*10,0),'or'),ylim([0 20]),xlim([0 inf])
    
% subplot(2,2,4)
% xlabel('Iterations','FontSize',10)
% ylabel('Conditioning weight')
% hold on
%     title(['Conditioning weight=',num2str(theta(4))])
%     plot(k,theta(4),'or'),ylim([0 1]),xlim([0 inf])
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