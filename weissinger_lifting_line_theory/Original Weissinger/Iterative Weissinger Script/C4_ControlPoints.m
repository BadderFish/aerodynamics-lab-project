function CPoint = C4_ControlPoints(config,internalMesh)
% This function finds the control points, located at 1/4 of the chord of
% the wing, at every section 

CPoint = cell(config.NBodies, 1);

for iBody = 1:config.NBodies
    CPoint{iBody, 1} = zeros(config.SemiSpanwiseDiscr(iBody), 3);
    for i = 1:config.SemiSpanwiseDiscr(iBody)
        CPoint_LE = (internalMesh{iBody,1}{1,i}.LERoot + internalMesh{iBody,1}{1,i}.LEtip)/2;
        CPoint_TE = (internalMesh{iBody,1}{end,i}.TERoot + internalMesh{iBody,1}{end,i}.TEtip)/2;
    
        versor = (CPoint_LE - CPoint_TE) / norm(CPoint_LE-CPoint_TE);
    
        CPoint{iBody, 1}(i,:) = CPoint_TE + (3/4 * norm(CPoint_LE - CPoint_TE) ) * versor;
        
    end 
    CPoint{iBody, 1} = [CPoint{iBody, 1}; flip([CPoint{iBody, 1}(:,1), -CPoint{iBody, 1}(:,2), CPoint{iBody,1}(:,3)])];
end

end

