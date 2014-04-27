clear all

Shape = 1; %[1 0045];
DomainSize = 4;
ref = 2;
powerRef = 2;
[TRItemp, Nodestemp, Toptemp, Bottomtemp, Lefttemp, Righttemp, InnerBound] = getDiscreteGeometry(Shape, DomainSize, ref, powerRef);

%Allow ranomization of node order to investigate sparse pattern
Nodes = Nodestemp;
TRI = TRItemp;
Top = Toptemp;
Bottom =Bottomtemp;
Left = Lefttemp;
Right = Righttemp;

for i=1:length(InnerBound)
    InnerBoundarytemp(i) = find(and(abs(Nodes(:,1)-InnerBound(i,1))<1e-10,abs(Nodes(:,2)-InnerBound(i,2))<1e-10));
end
InnerBoundary = InnerBoundarytemp;

if false %True for randomization, false to leave mesh alone
    Choices = 1:length(Nodes);
    for i=1:length(Nodes)
        Chooser = max(1,round(rand(1)*length(Choices)));
        Chosen = Choices(Chooser);
        Nodes(i,:) = Nodestemp(Chosen,:);
        TRI(TRItemp==Chosen) = i;
        Top(Toptemp==Chosen) = i;
        Bottom(Bottomtemp==Chosen) = i;
        Left(Lefttemp==Chosen) = i;
        Right(Righttemp==Chosen) = i;
        InnerBoundary(InnerBoundarytemp==Chosen) = i;
        Choices = Choices([1:Chooser-1 , Chooser+1:end]);
    end
end

Nodes = reshape(Nodes',1,numel(Nodes))';

dlmwrite('JJBElems',reshape(TRI',1,numel(TRI)), 'delimiter', '\n');
dlmwrite('JJBCold',unique([Left' Top' Right'])', 'delimiter', '\n');
dlmwrite('JJBHot',InnerBoundary', 'delimiter', '\n');
save('JJBNodes','Nodes','-ascii','-double');