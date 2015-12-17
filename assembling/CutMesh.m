%Returns and object representing a mesh cut by the incomining boundary
%polygon XB.
classdef CutMesh<handle
    properties
        %delauney triangulation object.
        dt;
        %Mesh-size.
        h;
        %matlab-cell of polygons representig the elements cut by the boundary.
        %Xcut{j} contains a polygon describing the cut element.
        %Xcut{j} is a m-by-2 matrix, m is either 4 or 3
        %since the cut with a triangle in 2d has either 4 or 3 corners.
        %Each column of Xcut defines the x and y coordinates of the
        %polygon. The first 2 points are located on the boundary.
        Xcut;
        %Faces of the elements that are intersected or cut by the boundary.
        assocFaces;
        %logical array defining which nodes are relevant for the
        %calculation
        relevant;
        %Nodes that are on the outer boundary of the structured mesh.
        outerNodes;
        %Boolean describing if we are trying to solve an inner or outer
        %problem.
        haveInnerProb;
    end
    
    properties(Access=private)
        %Logical array describing which triangles are inside the boundary
        %XB.
        inside;
        %Logical array describing which triangles are intersected by the 
        %boundary XB.
        intersected;
    end
    
    methods(Access=public)
        %Constructor, xLim,yLim are bounds for the structured background
        %mesh. n is a mesh-parameter. XB is a m-by-2 matrix describing the boundary
        %as a polygon. Each column of XB describes contains the x and y
        %coordinates of this polygon.
        %haveInnerProbelm is a boolean describing if 
        function[this]=CutMesh(xLim,yLim,n,XB,haveInnerProblem)
            %force approximately the same grid spacing.
            fracnxny=(xLim(2)-xLim(1))/(yLim(2)-yLim(1));
            ny=n;
            nx=floor(fracnxny*ny);
            createBackgroundMesh(this,xLim,yLim,nx,ny);
            checkIntersection(this,XB,haveInnerProblem);
            createBoundaryAssociatedFaces(this,haveInnerProblem);
            calculateRelevant(this,haveInnerProblem);
            calculateOuterNodes(this)
            this.haveInnerProb=haveInnerProblem;            
        end
                
        
        %Returns the elements that are intersected by the boundary.
        function[tri]=getTriIntersected(this)
            tri=this.dt.ConnectivityList(this.intersected,:);
        end
        
        %Returns the triangles relevant for the calculation which are not
        %intersected by the boundary.
        function[tri]=getTriInside(this)
            if(this.haveInnerProb)
                tri=this.dt.ConnectivityList(this.inside,:);
            else
                outside=(~this.inside)&(~this.intersected);
                tri=this.dt.ConnectivityList(outside,:);
            end
        end
        
        function[]=plot(this)
            hold on;
           for j=1:length(this.Xcut)
               xCut=this.Xcut{j};
               plot(xCut(:,1),xCut(:,2),'b-');
           end
           hold off;
        end
    end
    
    methods(Access=private)
        %Returns the nodes which are relevant for the calculation.
        function[]=calculateRelevant(this,haveInnerProblem)
            numNodes=size(this.dt.Points,1);
            this.relevant=false(numNodes,1);
            if(haveInnerProblem)
                %This corresponds to solving for the nodes of the triangles that are
                %inside or intersected by the boundary.
                relevantTriangles=this.inside|this.intersected;
            else%outerProblem
                %If we have an outer problem all nodes that are not inside
                %XB are relevant.
                relevantTriangles=~this.inside;
            end
            nodesInsideOrIntersected=unique(this.dt.ConnectivityList(relevantTriangles,:));
            this.relevant(nodesInsideOrIntersected)=true;
        end
        
        %Identifies the nodes that are on the outer boundary of the
        %background mesh.
        function[]=calculateOuterNodes(this)
            outerEdges=this.dt.freeBoundary;
            this.outerNodes=false(size(this.dt.Points,1),1);
            this.outerNodes(unique(outerEdges(:)))=true;
            %Cut to relevant
            this.outerNodes=this.outerNodes(this.relevant);
        end
        
        %Creates the delauneyTriangulation that is the background mesh.
        function[]=createBackgroundMesh(this,xlim,ylim,nx,ny)
            x=linspace(xlim(1),xlim(2),nx);
            y=linspace(ylim(1),ylim(2),ny);
            [Xm,Ym]=meshgrid(x,y);
            this.dt=delaunayTriangulation(Xm(:), Ym(:));
            %Mesh parameter.
            [~,r]=this.dt.incenter();
            this.h=2*max(r);
        end
        
        %Calculates the cell-array Xcut defining the cut-cells.
        function[]=checkIntersection(this,XB,haveInnerProblem)
            Xcenter=mean(XB);
            %Preallocate
            numElements=size(this.dt,1);
            this.inside=false(numElements,1);
            this.intersected=false(numElements,1);
            this.Xcut=cell(numElements,1);
            %Loop over all elements and check intersection.
            for j=1:numElements;
                Xelem=this.dt.Points(this.dt.ConnectivityList(j,:),:);
                %Check how many nodes are inside the boundary.
                nodesInside=inpolygon(Xelem(:,1),Xelem(:,2),XB(:,1),XB(:,2));
                numNodesInside=sum(nodesInside);
                %If all nodes are inside the boundary the element is inside
                %the boundary.
                this.inside(j)=numNodesInside==3;
                %An element is intersected if some but not all nodes are
                %inside the boundary.
                this.intersected(j)=numNodesInside>0 && ~(this.inside(j));
                if(this.intersected(j))
                    %Returns the cut element.
                    [this.Xcut{j}]=this.getSlicedElement(Xelem,XB,...
                        nodesInside,Xcenter,haveInnerProblem);
                end
            end
            %Make Xcut the size of the number of intersected elements.
            this.Xcut=this.Xcut(this.intersected);
        end
        
        
        
        %This function returns the faces known as FG in article by burman hansbo
        %inparameters are delauney triangulation, logical array with triangles that
        %are intersected and logical array specifying triangles that are inside the
        %fictious boundary.
        function[]=createBoundaryAssociatedFaces(this,haveInnerProblem)
            %Faces are only indentified by which points they connect to so I can create
            %to triangulations which have the same Points as dt but only consist
            %elements cut by the boundary. This should give a warning: disable this.
            warning off;
            %Make triangulation from triangles that are intersected by the boundary.
            trIsec=triangulation(this.dt.ConnectivityList(this.intersected,:),this.dt.Points);
            %Create a triangulation which has the outerboundary that I want
            %to remove from trIsec. This is different for an inner and outer
            %problem.
            if(haveInnerProblem)
                trFreeBoundary=triangulation(this.dt.ConnectivityList(this.intersected|this.inside,:),this.dt.Points);
            else
                trFreeBoundary=triangulation(this.dt.ConnectivityList(this.inside,:),this.dt.Points);
            end
            warning on;
            %The boundary faces of FG are the boundary faces of the
            %trFreeBoundary
            outerFaces=trFreeBoundary.freeBoundary;
            %take into account that one the boundary faces could be ordered
            %[1 2] and not [2 1].
            outerFaces=[outerFaces;fliplr(outerFaces)];
            %FG are the faces of trIsec setminus the outerfaces.
            this.assocFaces=setdiff(trIsec.edges,outerFaces,'rows');
            % % plotting for debugging.
%             close all;
%             figure;
%             triplot(this.dt,'b');
%             hold on;
%             circleHandle=ezplot('x.^2 +y.^2-1');
%             set(circleHandle,'linewidth',2);
%             for j=1:size(outerFaces,1)
%                 plot(this.dt.Points(outerFaces(j,:),1),this.dt.Points(outerFaces(j,:),2),'m',...
%                     'linewidth',2)
%             end
%             figure;
%             triplot(this.dt,'b');
%             hold on;
%             circleHandle=ezplot('x.^2 +y.^2-1');
%             set(circleHandle,'linewidth',2);
%             for j=1:size(this.assocFaces,1)
%                 plot(this.dt.Points(this.assocFaces(j,:),1),this.dt.Points(this.assocFaces(j,:),2),'r','linewidth',2)
%             end
        end
    end
    
    methods(Static,Access=private)
        %Calculates the Xcut for the incoming Xelem by taking the
        %intersection with the boundary XB. nodesInside is an array with
        %size of size(Xelem,2) definig which nodes are inside the boundary.
        %Xcenter is the center of XB.
        function[Xcut]=getSlicedElement(Xelem,XB,nodesInside,Xcenter,haveInnerProblem)
            XClosed=[Xelem;Xelem(1,:)];
            %Intersection of 2 polygons.
            [xIsec,yIsec] = polyxpoly(XB([1:end 1],1),XB([1:end 1],2),...
                XClosed(:,1), XClosed(:,2),'unique');
            Xisec=[xIsec,yIsec];
            XisecClosed=[Xisec;Xcenter];
            %Orientation should be counterclockwise for an inner problem
            %and clockwise for an outer problem. should concatenate inside
            %nodes for the outer problem and inside nodes for an inner
            %problem.
            if(haveInnerProblem)
                %toCat is which nodes should be concatenated with Xisec.
                %This is the inner nodes for an inner problem and the outer
                %nodes for an outer problem.
                toCat=nodesInside;
                flipOrientation=ispolycw(XisecClosed(:,1),XisecClosed(:,2));
            else
                toCat=~nodesInside;
                flipOrientation=~ispolycw(XisecClosed(:,1),XisecClosed(:,2));
            end
            if(flipOrientation)
                Xisec=flipud(Xisec);
            end
            %concatenate nodes to Xcut.
            Xcut=[Xisec;Xelem(toCat,:)];
            %Xelem(toCat) might not be ordered correctly, take convex cull
            %in order to make it so. This assumes that convhull will not
            %change the order of Xisec, check this condition at the end.
            order=convhull(Xcut(:,1),Xcut(:,2));
            order=order(1:end-1);
            Xcut=Xcut(order,:);
            if(order(1)~=1 || order(2)~=2)
                error('First points are not intersection points');
            end
            %Plotting for debugging:
            %             clf;
            %             hold on;
            %             plot(Xisec(:,1),Xisec(:,2),'ro');
            %             plot(Xcut(:,1),Xcut(:,2),'k');
            %             yLim=ylim();
            %             xLim=xlim();
            %             plot(XB(:,1),XB(:,2));
            %             ylim(yLim);
            %             xlim(xLim);
            %             plot(Xelem(:,1),Xelem(:,2),'b');
            %
        end
        
    end
end