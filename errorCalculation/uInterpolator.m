%For a given solution of the wave problem: u, dudt where u and dudt are vectors
%The class uInterpolator
%Creates an interpolation object that makes it possible to evaluate the u
%dudt and gradu everywhere over the background mesh.
classdef uInterpolator < handle
    properties(Access=private)
        %dealuneytriangulation object.
        dt;
        %a cell where each entry corresponds to a triangle in the mesh.
        %ABC{j} contains the vector which should be used to evaluate u on
        %triangle j, that is if x,y is located on triangle j we have:
        %abc=ABC{j}
        %u(x,y)=abc(1)+abc(2)*x+abc(3)*y
        ABC;
        %Works in the same way but for dudt.
        ABCdudt;
    end
    
    methods(Access=public)
        %Constructor, input parameters are delauney triangulation which can
        %be obtainted from cutMesh.dt, vector uSol and dudtSol, of size
        %sum(relevant) where relevant is the boolean vector specifying
        %which nodes in the background mesh that were used in the
        %calculation. uConst is what value the other values in the mesh 
        %should be set to.
        function[this]=uInterpolator(dTri,uSol,dudtSol,relevant,uConst)
            %unWrap u and dudt to the size of the number of nodes in the
            %background mesh.
            u=uConst*ones(size(dTri.Points,1),1);
            u(relevant)=uSol;
            dudt=zeros(size(dTri.Points,1),1);
            dudt(relevant)=dudtSol;
            %Save delauneytriangulation.
            this.dt=dTri;
            %Preallocate cells.
            numTri=size(dTri,1);
            this.ABC=zeros(3,numTri);
            this.ABCdudt=zeros(3,numTri);
            %Loop over all cells and calculate ABC and ABCdudt for each
            %triangle.
            for j=1:numTri
                tri=dTri.ConnectivityList(j,:);
                Xtri=dTri.Points(tri,:);
                uEncl=u(tri);
                dudtEncl=dudt(tri);
                %Solve a system LHS*ABC=uEnc to get the polynom on each
                %triangle.
                LHS=[ones(3,1) Xtri(:,1) Xtri(:,2)];
                this.ABC(:,j)=LHS\uEncl;
                this.ABCdudt(:,j)=LHS\dudtEncl;
            end
        end
        
        %Function to evaluate u at points x y.
        function[u]=evaluate(this,x,y)
            u=zeros(size(x));
            for i=1:size(x,1)
                for j=1:size(x,2)
                    %find which triangle is enclosing this triangle.
                    indexEncTri=this.dt.pointLocation(x(i,j),y(i,j));
                    abc=this.ABC(:,indexEncTri);
                    %interpolated u.
                    u(i,j)=abc(1)+abc(2)*x(i,j)+abc(3)*y(i,j);
                end
            end
        end
        
        %Function to evaluate dudt at points x y.
         function[dudt]=evaluatedudt(this,x,y)
            dudt=zeros(size(x));
            for i=1:size(x,1)
                for j=1:size(x,2)
                    %find which triangle is enclosing this triangle.
                    indexEncTri=this.dt.pointLocation(x(i,j),y(i,j));
                    abc=this.ABCdudt(:,indexEncTri);
                    %interpolated u.
                    dudt(i,j)=abc(1)+abc(2)*x(i,j)+abc(3)*y(i,j);
                end
            end
        end
        
        
        %Function to evaluate grad_u at points x y.
        %x and y must be vectors.
        function[gradu]=evaluategrad(this,x,y)
            %Force column vectors
            x=x(:);y=y(:);
            gradu=zeros(size(x,1),2);
            for i=1:size(x,1)
                %find which triangle is enclosing this triangle.
                indexEncTri=this.dt.pointLocation(x(i),y(i));
                abc=this.ABC(:,indexEncTri);
                %interpolated u.
                gradu(i,:)=[abc(2), abc(3)];
            end
        end
        
    end
end