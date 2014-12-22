%Wrapper for simplyfying plotting over the cutMesh. uConst is the level
%which the nodes that are not part of the calculation should take.
function[]=surfWrap(cutMesh,u,uConst)
uFull=uConst*ones(size(cutMesh.dt.Points,1),1);
uFull(cutMesh.relevant)=u;
trisurf(cutMesh.dt.ConnectivityList,cutMesh.dt.Points(:,1),cutMesh.dt.Points(:,2),uFull)%;,'edgecolor','none');
end