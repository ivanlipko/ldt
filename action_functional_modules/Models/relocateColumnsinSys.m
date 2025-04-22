%% 
% sys = relocateColumnsinSys(sys, relocColumns);
% 
function sys = relocateColumnsinSys(sys, cols)
    sys.A = sys.A(cols,cols);
%     sys.B = sys.B(cols,cols);
	if(size(sys.B,2) > 1)
        sys.B = sys.B(cols,:);
    else
        sys.B = sys.B(cols);
    end
    
    if size(sys.C,1) == size(sys.C,2)
        sys.C = sys.C(cols,cols);
    else
        sys.C = sys.C(:,cols);
    end
%     sys.D = sys.D(cols,cols);
	if(size(sys.D,2) > 1)
        sys.D = sys.D(cols,:);
    else
        if size(sys.D,1) == size(sys.A,1)
            sys.D = sys.D(cols);
        else
        end
    end
    sys.StateName = sys.StateName(cols);
    sys.StateUnit = sys.StateUnit(cols);
end
