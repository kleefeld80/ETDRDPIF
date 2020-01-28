% coarse to fine grid extractor for two dimensional domain with neuman
% boudary conditions. The fine grid is the reference solution for computing
% errors.
% Emmanuel Asante-Asamani
%11/05/14

function sol = finetocoarse(n,nf,g,solref)
% INPUTS n: the number of nodes along x direction of coarse grid
%        nf: the number of nodes along the x direction of fine grid
%        solref: the reference solution from which grid is to be extracted
%        g: the refinement level of the fine grid

% OUTPUTS sol: the extracted solution from the fine grid

Stemp = zeros(n);
SOL = reshape(solref,nf,nf);
SOL = SOL'; 
for i = 1:n
    for j = 1:n
        switch g
            case 1
               Stemp(i,j)=SOL(i,j);
            case 2
                Stemp(i,j)=SOL(2*i-1,2*j-1);
            case 3
                Stemp(i,j) = SOL(4*i-3,4*j-3);
            case 4
                Stemp(i,j) = SOL(8*i-7,8*j-7);
            case 5
                Stemp(i,j) = SOL(16*i-15,16*j-15);
            case 6
                Stemp(i,j) = SOL(32*i-31,32*j-31);
            case 7
                Stemp(i,j) = SOL(64*i-63,64*j-63);
            case 8
                Stemp(i,j) = SOL(128*i-127,128*j-127);
        end                        
    end
end
Stemp = Stemp';
sol = reshape(Stemp,n^2,1);