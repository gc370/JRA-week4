function [matrixuimplicit,errormatriximplicit] = week4parabolicimplicitq2withMLfortimederivative(xbeginning,xend,tbeginning,tend,g1,g2,eta,numberofpointsinx,numberofpointsint,actualsol)

    % ---- Inputs ---- %

    %   xbeginning        -> beginning point of the domain of x
    %   xend              -> end point of the domain of x
    %   tbeginning        -> beginning point of the domain of t (normally 0)
    %   tend              -> end point of the domain of x
    %   g1                -> @(t)   Neumann function at 0
    %   g2                -> Boundary value at the end point of the domain of x
    %   eta               -> function @(x) for intial conditions at t = tbeggining
    %   numberofpointsinx -> chosen number of nodal points such that numberofpointsinx + 1 is the total nodes including boundary points in space
    %   numberofpointsint -> chosen number of nodal p   oints such that numberofpointsint + 1 is the total nodes including boundary points in time
    %   f                 -> @(x,t) such that is the right hand side of the equation
    %   actualsol         -> @(x,t) the correct solution


    % ---- Mesh ---- %

    dx = (xend-xbeginning) / (numberofpointsinx);
    dt = (tend-tbeginning) / (numberofpointsint);

    pointx = zeros(1,numberofpointsinx+1);
    pointt = zeros(1,numberofpointsint+1);
    
            for i = 1:numberofpointsinx+1
                   pointx(i) = xbeginning + (i-1)*dx;
            end

            for j = 1:numberofpointsint+1
                   pointt(j) = tbeginning + (j-1)*dt;
            end


    % --- Construction of matrices of linear coefficients as in week 4 report (A and M) --- %

    matrixt = zeros(numberofpointsinx,numberofpointsinx);
    matrixm = dx*eye(numberofpointsinx,numberofpointsinx);

    % --- Here mass lumping has been used for matrix m --- %
    
    matrixt(1,1) = 1/dx;
    matrixt(1,2) = -1/dx;
    
    matrixt(numberofpointsinx,numberofpointsinx-1) = -1/dx;
    matrixt(numberofpointsinx,numberofpointsinx) = 2/dx;
    
    % ---- build the loop for the tri diagonal elements ---- %
    
            for i = 2: numberofpointsinx-1
    
                    matrixt(i,i-1) =-1/dx;
                    matrixt(i,i) =2/dx;
                    matrixt(i,i+1) =-1/dx;
        
            end
            
            
   % ---- build the matrix of results matrixu, and insert the boundary values and intial conditions. ---- %         
            
   matrixuimplicit = zeros(numberofpointsint+1,numberofpointsinx+1);
   
   % ---- Insert first row (t=0) using eta (intial condition) ---- %
   
           for i = 1:numberofpointsinx+1
                   matrixuimplicit(1,i) = eta(pointx(i)); 
           end        
           
           
  

   
   % ---- Start building the vector f for each iteration ---- %
   fvec = zeros(numberofpointsinx,1);
   matrixtobeinverted = ((matrixm) + (dt*matrixt));
   invertedmatrix = inv(matrixtobeinverted);
   
          for j = 2:numberofpointsint+1

                fvec(1)= -dt*g1(pointt(j)) ;
           
                    newresultsfornextiteration = invertedmatrix*((matrixm *(transpose(matrixuimplicit(j-1,1:numberofpointsinx)) + fvec)));        

                    for i=1:numberofpointsinx
            
                        matrixuimplicit(j,i) = newresultsfornextiteration(i);  
       
                    end
            
                matrixuimplicit(j,numberofpointsinx+1) = g2;
          end 
          
  % ---- Build matrix of correct solutions ---- %        
          
  correctsols = zeros(numberofpointsint+1,numberofpointsinx+1);
  
        for j=1:numberofpointsint+1
            for i = 1:numberofpointsinx+1
                correctsols(j,i) = actualsol(pointx(i),pointt(j));         
            end
        end
          
  % ---- Error matrix ---- %
  
  errormatriximplicit = abs(matrixuimplicit - correctsols);
          
end