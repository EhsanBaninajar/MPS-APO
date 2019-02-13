function WrGrid(grid,filename,namevar)

% See DS user guide for functions documentation

% Written by Ehsanollah Baninajar, 2015

fid=fopen(filename,'wt');

if fid==-1
    disp('cannot open file')
    fclose(fid);
    return
end

x=size(grid,1);
y=size(grid,2);
z=size(grid,3);
nvar=size(grid,4);
grid=reshape(grid,x*y*z,nvar);

fprintf(fid,'%i %i %i\n',x,y,z);
fprintf(fid,'%i\n',nvar);
for i=1:nvar
    if iscell(namevar)==1
        str=namevar{i};
    else
        str=namevar;
    end
    fprintf(fid,'%s\n',str);
end
fclose(fid);

save('Grid_body.txt','grid','-ascii', '-double', '-tabs');
str=['copy ',filename,'+Grid_body.txt ',filename];
system(str);
delete('Grid_body.txt');