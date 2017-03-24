function [La Ls] = resGen(FDO,a,e,N,dx,per)

idx = 1/dx;

La = zeros(N);
Ls = La;

CA = zeros(7,3);
CS = CA;

CA(1:7,1) = [ 0 0  -1 0  1  0 0]/2;
CA(1:7,2) = [ 0 1  -8 0  8 -1 0]/12;
CA(1:7,3) = [-1 9 -45 0 45 -9 1]/60;

CS(1:7,1) = [ 0  0  -1   2  -1  0  0];
CS(1:7,2) = [ 0  1  -4   6  -4  1  0];
CS(1:7,3) = [-1  6 -15  20 -15  6 -1];

CA =        a*idx*CA;
CS = e*abs(a)*idx*CS;

FDO2 = FDO/2;

if per == 1
    for j = 0:N-1
        
            La(j+1, mod(j-3,N)+1 ) = CA(1,FDO2);
            La(j+1, mod(j-2,N)+1 ) = CA(2,FDO2);
            La(j+1, mod(j-1,N)+1 ) = CA(3,FDO2);
            La(j+1, mod(j-0,N)+1 ) = CA(4,FDO2);
            La(j+1, mod(j+1,N)+1 ) = CA(5,FDO2);
            La(j+1, mod(j+2,N)+1 ) = CA(6,FDO2);
            La(j+1, mod(j+3,N)+1 ) = CA(7,FDO2);

            Ls(j+1, mod(j-3,N)+1 ) = CS(1,FDO2);
            Ls(j+1, mod(j-2,N)+1 ) = CS(2,FDO2);
            Ls(j+1, mod(j-1,N)+1 ) = CS(3,FDO2);
            Ls(j+1, mod(j-0,N)+1 ) = CS(4,FDO2);
            Ls(j+1, mod(j+1,N)+1 ) = CS(5,FDO2);
            Ls(j+1, mod(j+2,N)+1 ) = CS(6,FDO2);
            Ls(j+1, mod(j+3,N)+1 ) = CS(7,FDO2);

    end
else
    
    % Interior
    for j = FDO2:N-FDO2-1
      
            La(j+1, mod(j-3,N)+1 ) = CA(1,FDO2);
            La(j+1, mod(j-2,N)+1 ) = CA(2,FDO2);
            La(j+1, mod(j-1,N)+1 ) = CA(3,FDO2);
            La(j+1, mod(j-0,N)+1 ) = CA(4,FDO2);
            La(j+1, mod(j+1,N)+1 ) = CA(5,FDO2);
            La(j+1, mod(j+2,N)+1 ) = CA(6,FDO2);
            La(j+1, mod(j+3,N)+1 ) = CA(7,FDO2);

            Ls(j+1, mod(j-3,N)+1 ) = CS(1,FDO2);
            Ls(j+1, mod(j-2,N)+1 ) = CS(2,FDO2);
            Ls(j+1, mod(j-1,N)+1 ) = CS(3,FDO2);
            Ls(j+1, mod(j-0,N)+1 ) = CS(4,FDO2);
            Ls(j+1, mod(j+1,N)+1 ) = CS(5,FDO2);
            Ls(j+1, mod(j+2,N)+1 ) = CS(6,FDO2);
            Ls(j+1, mod(j+3,N)+1 ) = CS(7,FDO2);
           
        
    end
    
    % Boundaries
    if FDO > 2
        for j = [2 N-1]
           
                
                La(j, j-1) = CA(3,1);
                La(j, j  ) = CA(4,1);
                La(j, j+1) = CA(5,1);
                
                Ls(j, j-1) = CS(3,1);
                Ls(j, j  ) = CS(4,1);
                Ls(j, j+1) = CS(5,1);
                
           
        end
   
        if FDO > 4
            for j = [3 N-2]
                    
                    La(j, j-2) = CA(2,2);
                    La(j, j-1) = CA(3,2);
                    La(j, j)   = CA(4,2);
                    La(j, j+1) = CA(5,2);
                    La(j, j+2) = CA(6,2);
                    
        
                    Ls(j, j-2) = CS(2,2);
                    Ls(j, j-1) = CS(3,2);
                    Ls(j, j)   = CS(4,2);
                    Ls(j, j+1) = CS(5,2);
                    Ls(j, j+2) = CS(6,2);
                    
            end
        end
    end
    
    
    for j = N
         La(j,j)   =  a*idx;
         La(j,j-1) = -a*idx;
    end
   
    
end

%Lspace = (La + Ls);