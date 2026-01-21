%% WEISSINGER METHOD --> Iterative process
close all
clear
clc

%% Angle of attack vector
alpha_vect = 0:1:5;
config.CL = zeros(length(alpha_vect),1);
config.CD_ind = zeros(length(alpha_vect),1);
M = 20;
N = 10;
config.alpha_ind = zeros(length(alpha_vect),2*M);


for i = 1:length(alpha_vect)

    % Printf Alpha value
    fprintf('--- For alpha = %.3f deg ---\n',alpha_vect(i));

    %% Test Case: Cessna 172
    U_Inf_Mag = 62;          % [m/s]
    alpha     = 3;           % [deg]
    beta      = 0;           % [deg]

    U_Inf = [cosd(beta) sind(beta) 0] .* U_Inf_Mag;
    rho   = 1.225;

    config.NBodies = 2;

    config.RootChord     = [1.487,1.3];   % [m]
    config.DihedralAngle = [1.75,0];      % [deg]
    config.SweepAngle    = [0,0];         % [deg]
    config.TaperRatio    = [1,0.615];
    config.Span          = [10.89,3.5];   % [m]

    config.LEPosition_X = [0,4.5*cosd(alpha_vect(i))];        % [m]
    config.LEPosition_Y = [0,0];
    config.LEPosition_Z = [0,-0.5*(1 + sind(alpha_vect(i)))];

    config.RotationAngle_X = [0,0];
    config.RotationAngle_Y = [alpha_vect(i),alpha_vect(i)-2];
    config.RotationAngle_Z = [0,0];

    % Discretization
    config.SemiSpanwiseDiscr = [M,M];
    config.ChordwiseDiscr    = [N,N];

    %% Preliminary computations
    config.SemiSpan = config.Span ./ 2;

    config.Surface = 2 * ( ...
        config.SemiSpan .* ...
        config.RootChord .* ...
        (1 + config.TaperRatio) ./ 2 );

    config.SurfaceProjected = config.Surface .* cosd(config.DihedralAngle);

    config.TipChord = config.RootChord .* config.TaperRatio;

    config.MAC = (2/3) .* config.RootChord .* ...
        ( (1 + config.TaperRatio + config.TaperRatio.^2) ./ ...
          (1 + config.TaperRatio) );

    config.AspectRatio = config.Span.^2 ./ config.Surface;

    %% Geometry structures
    ControlPoints     = cell(config.NBodies,1);
    InducedPoints     = cell(config.NBodies,1);
    Normals           = cell(config.NBodies,1);
    InfiniteVortices  = cell(config.NBodies,1);
    Vortices          = cell(config.NBodies,1);
    internalMesh      = cell(config.NBodies,1);
    WingExtremes      = cell(config.NBodies,1);

    for iBody = 1:config.NBodies
        [ ControlPoints{iBody}, ...
          InducedPoints{iBody}, ...
          Normals{iBody}, ...
          InfiniteVortices{iBody}, ...
          Vortices{iBody}, ...
          internalMesh{iBody}, ...
          WingExtremes{iBody} ] = createStructure(config,iBody);
    end

    %% Matrix initialization
    NPanelsTot = 2 * config.SemiSpanwiseDiscr * config.ChordwiseDiscr';
    matrixA   = zeros(NPanelsTot);
    knownTerm = zeros(NPanelsTot,1);

    %% Influence matrix construction
    rowIndex = 0;

    for iBody = 1:config.NBodies
        for iChord = 1:config.ChordwiseDiscr(iBody)
            for iSpan = 1:2*config.SemiSpanwiseDiscr(iBody)

                rowIndex = rowIndex + 1;
                columnIndex = 0;

                CP = ControlPoints{iBody}{iChord,iSpan}.Coords;
                n  = Normals{iBody}{iChord,iSpan}.Coords;

                for jBody = 1:config.NBodies
                    for jChord = 1:config.ChordwiseDiscr(jBody)
                        for jSpan = 1:2*config.SemiSpanwiseDiscr(jBody)

                            columnIndex = columnIndex + 1;

                            % Root semi-infinite vortex
                            E1 = InfiniteVortices{jBody}{jChord,jSpan}.Root.toInfty;
                            E2 = InfiniteVortices{jBody}{jChord,jSpan}.Root.onWing;
                            U  = vortexInfluence(CP,E1,E2);

                            % Bound vortex
                            E1 = Vortices{jBody}{jChord,jSpan}.Root;
                            E2 = Vortices{jBody}{jChord,jSpan}.Tip;
                            U  = U + vortexInfluence(CP,E1,E2);

                            % Tip semi-infinite vortex
                            E1 = InfiniteVortices{jBody}{jChord,jSpan}.Tip.onWing;
                            E2 = InfiniteVortices{jBody}{jChord,jSpan}.Tip.toInfty;
                            U  = U + vortexInfluence(CP,E1,E2);

                            matrixA(rowIndex,columnIndex) = dot(U,n);
                        end
                    end
                end
            end
        end
    end

    %% Known term
    rowIndex = 0;
    for iBody = 1:config.NBodies
        for iChord = 1:config.ChordwiseDiscr(iBody)
            for iSpan = 1:2*config.SemiSpanwiseDiscr(iBody)

                rowIndex = rowIndex + 1;
                n = Normals{iBody}{iChord,iSpan}.Coords;
                knownTerm(rowIndex) = -dot(U_Inf,n);

            end
        end
    end

    %% Solve system
    Solution = linsolve(matrixA,knownTerm);

    Gamma = cell(config.NBodies,1);
    rowIndex = 0;

    for iBody = 1:config.NBodies
        Gamma{iBody} = zeros( ...
            config.ChordwiseDiscr(iBody), ...
            2*config.SemiSpanwiseDiscr(iBody) );

        for iChord = 1:config.ChordwiseDiscr(iBody)
            for iSpan = 1:2*config.SemiSpanwiseDiscr(iBody)

                rowIndex = rowIndex + 1;
                Gamma{iBody}(iChord,iSpan) = Solution(rowIndex);

            end
        end
    end

    %% Lift
    [L_2D,L,delta_b] = Lift_fun(config,rho,U_Inf_Mag,Gamma,0);
    config.CL(i) = L / (0.5*rho*U_Inf_Mag^2*config.Surface(1));
    fprintf('CL = %.3f\n',config.CL(i));

    %% Induced drag
    CPoint = C4_ControlPoints(config,internalMesh);

    [alpha_ind,Di2D,Di] = Induced_Drag_fun( ...
        config,Normals,L_2D,CPoint,InfiniteVortices,U_Inf,Gamma,0);
    % config.alpha_ind(i,:) = alpha_ind;

    config.CD_ind(i) = Di / (0.5*rho*U_Inf_Mag^2*config.Surface(1));
    fprintf('CDi = %.3f\n',config.CD_ind(i));
end

save Wessinger_complete.mat config