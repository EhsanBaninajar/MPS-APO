% In case of random gapsize
% gapsize_min=10;
% gapsize_max=20;
% gapsize=randi([gapsize_min,gapsize_max]);

TiVar=1; % Which Ti variable will be used to calculate te contrast?(TiVar)
gapid=0;
Ti_sort=Ti;
gapdistance=10;
gaptoppercent=20; % the top percentage of gap variability
%% preallocating
std_gap=zeros((Xsize_Ti-gapsize-2*gapdistance+1)*(Ysize_Ti-gapsize-2*gapdistance+1),1);
%% gap std database
if Zsize_Ti==1
    for i=1+gapdistance:(Xsize_Ti-gapsize-gapdistance+1)
        for j=1+gapdistance:(Ysize_Ti-gapsize-gapdistance+1)
            gapid=gapid+1;
            std_gap(gapid,:)=std2(Ti_sort(i:(i+gapsize-1),j:(j+gapsize-1),:,TiVar));
        end
    end
end
[sorted_std,sorted_std_index]=sortrows(-1.*std_gap(:,1));
top_gaps=fix(((Xsize_Ti-gapsize-gapdistance+1)*(Ysize_Ti-gapsize-gapdistance+1)*gaptoppercent)/100);
sorted_std_index=sorted_std_index(1:top_gaps);