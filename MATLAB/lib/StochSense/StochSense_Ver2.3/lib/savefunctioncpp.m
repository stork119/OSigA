function savefunctioncpp(rhs,name)


sizef = size(rhs);

dim1=sizef(1);

dim2=sizef(2);

R=rhs


file = fopen(name,'w');

% writing the jacobian
fprintf(file,'%s\n\n','name');
%add the force options
%fprintf(file, '%%Model jacobian dy/dy\n\n');
fprintf(file,'%s\n');
for i=1:dim1  
    for j=1:dim2
        fprintf(file,['m(',num2str(i-1),',', num2str(j-1),') = ']);
        fprintf(file,' %s',char(R(i,j)));
        fprintf(file,';\n');
    end
end   
fclose(file);
