% code copied from this project:
% https://github.com/Bounciness/Volume-and-Sampling
% https://link.springer.com/content/pdf/10.1007/s12532-015-0097-z.pdf

% we convert it to the .ine format, so that we can read it.

for dim=3:20
    fileID = fopen(sprintf('birk%d.ine',dim),'w');
    A=-eye(dim^2);
    C=zeros(2*dim,dim^2);
    for i=1:dim
        %row constraints
        for j=(i-1)*dim+1:i*dim
          C(i,j)=1; 
        end
        
        %col constraints
        for j=i:dim:dim^2
           C(dim+i,j)=1; 
        end
    end
    
    nullC=null(C);
    %let v be any solution to Cx=d
    v=ones(dim^2,1)/dim;
    
    new_A=A*nullC;
    new_b=zeros(dim^2,1)-A*v;
    
    fprintf(fileID,'birk%d.ine\n',dim);
    fprintf(fileID,'H-representation\n');
    fprintf(fileID,'begin\n');
    fprintf(fileID,' %d %d real\n',size(new_A,1),size(new_A,2)+1);
    for i=1:size(new_A,1)
        fprintf(fileID,'%.20e ', new_b(i));
        for j=1:size(new_A,2)
            fprintf(fileID,'%.20e ', new_A(i,j));
        end
        fprintf(fileID,'\n');
    end
    fprintf(fileID,'end\n');
    fprintf(fileID,'input_incidence\n');
    fclose(fileID);
end

fprintf('done.\n');