function [matrixuimplicitwithoutml,errormatriximplicitwithoutml] = week4parabolicimplicitdirichlet(xbeginning,xend,tbeginning,tend,g1,g2,eta,numberofpointsinx,numberofpointsint,f,actualsol)

    % ---- Inputs ---- %

    %   xbeginning        -> beginning point of the domain of x
    %   xend              -> end point of the domain of x
    %   tbeginning        -> beginning point of the domain of t (normally 0)
    %   tend              -> end point of the domain of x
    %   g1                -> Boundary value at the beggining point of the domain of x
    %   g2                -> Boundary value at the end point of the domain of x
    %   eta               -> function @(x) for intial conditions at t = tbeggining
    %   numberofpointsinx -> chosen number of nodal points such that numberofpointsinx + 1 is the total nodes including boundary points in space
    %   numberofpointsint -> chosen number of nodal points such that numberofpointsint + 1 is the total nodes including boundary points in time
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

    matrixa = zeros(numberofpointsinx-1,numberofpointsinx-1);
    matrixm = zeros(numberofpointsinx-1,numberofpointsinx-1);
    
    matrixa(1,1) = 2/dx;
    matrixa(1,2) = -1/dx;
    matrixa(numberofpointsinx-1,numberofpointsinx-2) = -1/dx;
    matrixa(numberofpointsinx-1,numberofpointsinx-1) = 2/dx;
    
    matrixm(1,1) = (2*dx)/3;
    matrixm(1,2) = dx/6;
    matrixm(numberofpointsinx-1,numberofpointsinx-2) = dx/6;
    matrixm(numberofpointsinx-1,numberofpointsinx-1) = (2*dx)/3;
    
    % ---- build the loop for the tri diagonal elements ---- %
    
            for i = 2: numberofpointsinx-2
    
                    matrixa(i,i-1) =-1/dx;
                    matrixa(i,i) =2/dx;
                    matrixa(i,i+1) =-1/dx;
        
                    matrixm(i,i-1) =dx/6;
                    matrixm(i,i) =(2/3)*dx;
                    matrixm(i,i+1) =dx/6;
        
            end
            
            
   % ---- build the matrix of results matrixu, and insert the boundary values and intial conditions. ---- %         
            
   matrixuimplicitwithoutml = zeros(numberofpointsint+1,numberofpointsinx+1);
   
   % ---- Insert first row (t=0) using eta (intial condition) ---- %
   
           for i = 1:numberofpointsinx+1
                   matrixuimplicitwithoutml(1,i) = eta(pointx(i)); 
           end        
           
           
   % ---- Build Euler's explicit scheme FEM ---- %        

   
   % ---- Start building the vector f for each iteration ---- %
   fvec = zeros(numberofpointsinx-1,1);
   matrixtobeinverted = (matrixm + (dt*matrixa));
   invertedmatrix = inv(matrixtobeinverted);
   
          for j = 2:numberofpointsint+1

                fvec(1)= dt*dx*f(pointx(2),pointt(j)) + (g1*(dt/dx)) ;
                fvec(numberofpointsinx-1)= dt*dx*f(pointx(numberofpointsinx),pointt(j)) + (g2*(dt/dx));

                    for i = 2: numberofpointsinx-2 
                        fvec(i) = dt*dx*f(pointx(i+1),pointt(j));
                    end
           
                    newresultsfornextiteration = invertedmatrix*((matrixm *(transpose(matrixuimplicitwithoutml(j-1,2:numberofpointsinx)) + fvec)));        

                    for i=1:numberofpointsinx-1
            
                        matrixuimplicitwithoutml(j,i+1) = newresultsfornextiteration(i);  
       
                    end
            
                matrixuimplicitwithoutml(j,1) = g1;
                matrixuimplicitwithoutml(j,numberofpointsinx+1) = g2;
          end 
          
  % ---- Build matrix of correct solutions ---- %        
          
  correctsols = zeros(numberofpointsint+1,numberofpointsinx+1);
  
        for j=1:numberofpointsint+1
            for i = 1:numberofpointsinx+1
                correctsols(j,i) = actualsol(pointx(i),pointt(j));         
            end
        end
          
  % ---- Error matrix ---- %
  
  errormatriximplicitwithoutml = abs(matrixuimplicitwithoutml - correctsols);
          
end
