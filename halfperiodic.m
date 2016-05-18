classdef halfperiodic < handle
    methods(Access=public)
        function[this]=halfperiodic(n)
            if(nargin<1)
                n=22;
            end
            addpath analyticFunctions/plane;
            addpath assembling/;
            addpath errorCalculation;
            this.xLim=[-1.5,1.5];
            this.createMesh(this.xLim,n);
            this.nodes=this.create_sets(this.cutMesh);
            this.createSystem();            

        end
        
        function[t,u,dudt]=solve(this,waveNumber,endTime,nMax,saveSolution)
            [u0,du0dt]=this.setupInitialConditions(waveNumber);
            domainWidth=diff(this.xLim);
            [dt,nSteps]=getnSteps(this.cfl,this.waveSpeed,endTime,nMax,domainWidth);
            [u,dudt,t]=timeStepRK(u0,du0dt,this.A,this.M,this.L,...
                nSteps,dt,this.waveSpeed);
            u=halfperiodic.convertToFull(u,this.nodes);
            dudt=halfperiodic.convertToFull(dudt,this.nodes);
            this.saveIfNecessary(saveSolution,t,u,dudt,dt)
        end
    end
    
    properties(Access=public)
        xLim;
        yLim;
        cutMesh;
        nodes;
        waveSpeed=1;
        cfl=0.4;
        M;
        A;
        L;
    end
    
    methods(Access=private)
        
        function[]=saveIfNecessary(this,saveSolution,t,u,dudt,dt)
            if(saveSolution)
                uConst=0;
                saveName=['periodic_' '.mat'];
                c=this.waveSpeed;
                cutMesh=this.cutMesh;
                xLim=this.xLim;
                yLim=this.yLim;
                save(saveName,'t','u','dudt','cutMesh','xLim','yLim','c',...
                    'dt','uConst');
            end
            
        end
        function[]=createMesh(this,xLim,n)
            yMeshLim=xLim;
            H=diff(xLim)/(n-1);
            epsilon=.5;
            xPolygon=[-2;2];
            this.yLim=[yMeshLim(1)+epsilon*H,yMeshLim(2)-epsilon*H];
            [xB,yB]=meshgrid(xPolygon,this.yLim);
            order=convhull(xB,yB);
            XB=[xB(:),yB(:)];
            XB=XB(order,:);
            haveInnerProblem=true;
            this.cutMesh=CutMesh(xLim,yMeshLim,n,XB,haveInnerProblem);
        end
        
        function[]=createSystem(this)
            gD=@(x,y) zeros(size(x));
            f=@(x,y) zeros(size(x));
            dirichletInner=false;
            dirichletOuter=false;%Not used only for an outer problem.
            [MGlobal,AGlobal]=assemble(this.cutMesh,f,gD,dirichletInner,gD,dirichletOuter);
            this.M=halfperiodic.makeMatrixPeriodic(MGlobal,this.nodes);
            this.A=halfperiodic.makeMatrixPeriodic(AGlobal,this.nodes);
            this.L=@(t) zeros(size(this.A,1),1);
        end
        
        function[u0,du0dt]=setupInitialConditions(this,waveNumber)
            t=0;
            uAnaly=@(x,y) planeWave(x,y,t,waveNumber,this.waveSpeed);
            dudtAnaly=@(x,y) dplaneWavedt(x,y,t,waveNumber,this.waveSpeed);
            u0=uAnaly(this.cutMesh.dt.Points(:,1),this.cutMesh.dt.Points(:,2));
            du0dt=dudtAnaly(this.cutMesh.dt.Points(:,1),this.cutMesh.dt.Points(:,2));
            u0=halfperiodic.removeRightBoundaryDofs(u0,this.nodes);
            du0dt=halfperiodic.removeRightBoundaryDofs(du0dt,this.nodes);
        end
    end
    
    methods(Access=private,Static)    
        function[uFull]=convertToFull(u,nodes)
            nInternalNodes=length(nodes.internal);
            nBoundaryNodes=length(nodes.leftBoundary);
            nNodesInPeriodic=nInternalNodes+nBoundaryNodes;
            nNodesInFull=nInternalNodes+2*nBoundaryNodes;
            %Make sure u is the right size.
            assert(size(u,1)==nInternalNodes+nBoundaryNodes);
            internalIndices=1:nInternalNodes;
            boundaryIndices=(nInternalNodes+1):nNodesInPeriodic;
            uFull=zeros(nNodesInFull,size(u,2));
            uFull(nodes.internal,:)=u(internalIndices,:);
            uFull(nodes.leftBoundary,:)=u(boundaryIndices,:);
            uFull(nodes.rightBoundary,:)=u(boundaryIndices,:);
        end
        
        function[u]=removeRightBoundaryDofs(u,nodes)
            nInternalNodes=length(nodes.internal);
            nBoundaryNodes=length(nodes.leftBoundary);
            assert(length(u)==nInternalNodes+2*nBoundaryNodes);
            u=[u(nodes.internal);u(nodes.leftBoundary)];
        end
        
        function[matrix]=makeMatrixPeriodic(globalMatrix,nodes)
            nInternal=length(nodes.internal);
            nBoundary=length(nodes.leftBoundary);
            nDofs=nInternal+nBoundary;
            matrix=sparse(nDofs,nDofs);
            upperIndices=1:nInternal;
            lowerIndices=(nInternal+1):nDofs;
            matrix(upperIndices,upperIndices)=globalMatrix(nodes.internal,nodes.internal);
            matrix(upperIndices,lowerIndices)=globalMatrix(nodes.internal,nodes.leftBoundary)+...
                globalMatrix(nodes.internal,nodes.rightBoundary);
            matrix(lowerIndices,upperIndices)=globalMatrix(nodes.leftBoundary,nodes.internal)+...
                globalMatrix(nodes.rightBoundary,nodes.internal);
            matrix(lowerIndices,lowerIndices)=...
                globalMatrix(nodes.leftBoundary,nodes.leftBoundary)+...
                globalMatrix(nodes.leftBoundary,nodes.rightBoundary)+...
                globalMatrix(nodes.leftBoundary,nodes.rightBoundary)+...
                globalMatrix(nodes.rightBoundary,nodes.rightBoundary);
        end
        
        function[nodes]=create_sets(cutMesh)
            h=cutMesh.h;
            triangulation=cutMesh.dt;
            eps=.1*h;
            freeBoundaryEdges=triangulation.freeBoundary;
            boundaryNodes=unique(freeBoundaryEdges(:));
            boundaryPoints=triangulation.Points(boundaryNodes,:);
            %Find nodes on the right boundary
            isOnRightBoundary=(max(boundaryPoints(:,1))-eps) < boundaryPoints(:,1);
            rightBoundary=boundaryNodes(isOnRightBoundary);
            %Find nodes on left boundary
            isOnLeftBoundary=boundaryPoints(:,1)< (min(boundaryPoints(:,1))+eps);
            leftBoundary=boundaryNodes(isOnLeftBoundary);
            %Find all nodes in between left and right boundary.
            internal=halfperiodic.getInternalNodes(triangulation,leftBoundary,rightBoundary);
            nodes=struct('leftBoundary',leftBoundary,'rightBoundary',rightBoundary,...
                'internal',internal);
        end
        
        function[internalNodes]=getInternalNodes(triangulation,leftBoundary,rightBoundary)
            allNodes=1:size(triangulation.Points,1);
            toKeep=true(size(triangulation.Points,1),1);
            toKeep(rightBoundary)=false;
            toKeep(leftBoundary)=false;
            internalNodes=allNodes(toKeep);
        end
    end
end
