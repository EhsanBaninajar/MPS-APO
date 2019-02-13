% In case of random gapsize
% gapsize_min=10;
% gapsize_max=20;
% gapsize=randi([gapsize_min,gapsize_max]);

%preallocating
std_gap=zeros((Xsize_Ti-gapsize+1)*(Ysize_Ti-gapsize+1),1);
%% gap std database
% Which Ti variable will be used to calculate te contrast?(TiVar)
TiVar=1;
gapid=0;
Ti_sort=Ti;

if Zsize_Ti==1
    for i=1:(Xsize_Ti-gapsize+1)
        for j=1:(Ysize_Ti-gapsize+1)
            gapid=gapid+1;
            std_gap(gapid,:)=std2(Ti_sort(i:(i+gapsize-1),j:(j+gapsize-1),:,TiVar));
        end
    end
else
    for i=1:(Xsize_Ti-gapsize+1)
        for j=1:(Ysize_Ti-gapsize+1)
            for k=1:(Zsize_Ti-gapsize+1)
                gapid=gapid+1;
                Ti_sort3d=Ti_sort(i:(i+gapsize-1),j:(j+gapsize-1),k:(k+gapsize-1),TiVar);
                std_gap(gapid,:)=std(Ti_sort3d(:));
            end
        end
    end
end
[sorted_std,sorted_std_index]=sortrows(-1.*std_gap(:,1));
top_gaps=fix(((Xsize_Ti-gapsize+1)*(Ysize_Ti-gapsize+1)*20)/100);
sorted_std_index=sorted_std_index(1:top_gaps);