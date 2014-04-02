Shape = 1; %[1 0045];
DomainSize = 4;
ref = 2;
powerRef = 2;
[TRItemp, Nodestemp, Toptemp, Bottomtemp, Lefttemp, Righttemp, InnerBoundarytemp] = getDiscreteGeometry(Shape, DomainSize, ref, powerRef);

%Allow ranomization of node order to investigate sparse pattern
Nodes = Nodestemp;
TRI = TRItemp;
Top = Toptemp;
Bottom =Bottomtemp;
Left = Lefttemp;
Right = Righttemp;
InnerBoundary = InnerBoundarytemp;

if false %True for ranomization, false to leave mesh alone
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
dlmwrite('JJBLeft',Left, 'delimiter', '\n');
dlmwrite('JJBRight',Right, 'delimiter', '\n');
save('JJBNodes','Nodes','-ascii','-double');