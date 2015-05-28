function [g1L, g2L, g3L, g4L] = sdlfmXsdlfmvKernGradientICBlock(lfmKern1, ...
    lfmKern2, tInit1, tInit2, kyy, kyv, kvy, kvv, i, j, generalConst, ...
    generalConstGrad, gkyy1, gkyy2, gkyy3, gkyy4, gkyv1, gkyv2, ...
    gkyv3, gkyv4, gkvy1, gkvy2, gkvy3, gkvy4, gkvv1, gkvv2, gkvv3, gkvv4)

% SDLFMXSDLFMVKERNGRADIENTICBLOCK Partial derivatives initial conditions
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN


g1L = cell(2);
g2L = cell(2);
g3L = cell(2);
g4L = cell(2);

typeParam = {'sdlfm','sdlfmv'};

% Compute for tInit1(1) and tInit2(2)

[g1L{1,2}, g2L{1,2}, g3L{1,2}] = sdlfmXsdlfmvKernGradientBlock(lfmKern1, lfmKern2, tInit1(1), ...
    tInit2(2), kyy, kyv, kvy, kvv, i, j, generalConst, generalConstGrad, 1);

g1Plus = sdlfmvXsdlfmvKernComputeBlock(lfmKern1, lfmKern2, tInit1(1), ...
    tInit2(2), kyy, kyv, kvy, kvv, i, j, generalConst);

g2Plus = sdlfmaXsdlfmKernComputeBlock(lfmKern2, lfmKern1, tInit2(2), ...
    tInit1(1), kyy, kvy, kyv, kvv, j, i, generalConst');

if (i==1 && j==1) || (i==1) && (j~=1)

    if i~=1 || j~=1
        [gkyy1IC, gkyy2IC, gkyy3IC, gkyy4IC] = sdlfmXsdlfmKernGradientIC(lfmKern1, lfmKern2, ...
            tInit1(1), tInit2(2), gkyy1{i,j}, gkyy2{i,j}, gkyy3{i,j}, gkyy4{i,j}, gkyv1{i,j}, gkyv2{i,j}, ...
            gkyv3{i,j}, gkyv4{i,j}, gkvy1{i,j}, gkvy2{i,j}, gkvy3{i,j}, gkvy4{i,j}, gkvv1{i,j}, gkvv2{i,j}, ...
            gkvv3{i,j}, gkvv4{i,j}, 1, typeParam);
        for k=1:3
            g1L{1,2}{k} = g1L{1,2}{k} + gkyy1IC{k};
            g2L{1,2}{k} = g2L{1,2}{k} + gkyy2IC{k};
        end
        cLength = length(gkyy3IC);
        g3L{1,2}(1:cLength) = g3L{1,2}(1:cLength) + gkyy3IC;
        g3L{1,2}(i)   = g3L{1,2}(i) + g1Plus;
        g3L{1,2}(j+1) = g2Plus;
        g4L{1,2} = gkyy4IC;
    else
        g3L{1,2}(i)   = g3L{1,2}(i) + g1Plus;
        g3L{1,2}(j+1) = g2Plus;
        g4L{1,2} = sdlfmKernMeanCovPartial(lfmKern1(1), lfmKern2(1), tInit1(1), ...
            tInit2(2), 1, typeParam);
    end
end

% Compute for tInit1(2) and tInit2(1)

[g1L{2,1}, g2L{2,1}, g3L{2,1}] = sdlfmXsdlfmvKernGradientBlock(lfmKern1, lfmKern2, tInit1(2), ...
    tInit2(1), kyy, kyv, kvy, kvv, i, j, generalConst, generalConstGrad, 1);


g1Plus = sdlfmvXsdlfmvKernComputeBlock(lfmKern1, lfmKern2, tInit1(2), ...
    tInit2(1), kyy, kyv, kvy, kvv, i, j, generalConst);

g2Plus = sdlfmaXsdlfmKernComputeBlock(lfmKern2, lfmKern1, tInit2(1), ...
    tInit1(2), kyy, kvy, kyv, kvv, j, i, generalConst');

if (i==1 && j==1) || (i~=1) && (j==1)
    if i~=1 || j~=1
        [gkyy1IC, gkyy2IC, gkyy3IC, gkyy4IC] = sdlfmXsdlfmKernGradientIC(lfmKern1, lfmKern2, ...
            tInit1(2), tInit2(1), gkyy1{i,j}, gkyy2{i,j}, gkyy3{i,j}, gkyy4{i,j}, gkyv1{i,j}, gkyv2{i,j}, ...
            gkyv3{i,j}, gkyv4{i,j}, gkvy1{i,j}, gkvy2{i,j}, gkvy3{i,j}, gkvy4{i,j}, gkvv1{i,j}, gkvv2{i,j}, ...
            gkvv3{i,j}, gkvv4{i,j}, 1, typeParam);
        for k=1:3
            g1L{2,1}{k} = g1L{2,1}{k} + gkyy1IC{k};
            g2L{2,1}{k} = g2L{2,1}{k} + gkyy2IC{k};
        end
        cLength = length(gkyy3IC);
        g3L{2,1}(1:cLength) = g3L{2,1}(1:cLength) + gkyy3IC;
        g3L{2,1}(j)   =  g3L{2,1}(j) + g2Plus;
        g3L{2,1}(i+1) =  g1Plus;
        g4L{2,1} = gkyy4IC;
    else
        g3L{2,1}(j)   =  g3L{2,1}(j) + g2Plus;
        g3L{2,1}(i+1) =  g1Plus;
        g4L{2,1} = sdlfmKernMeanCovPartial(lfmKern1(1), lfmKern2(1), tInit1(2), ...
            tInit2(1), 1, typeParam);
    end
end

% Compute for tInit1(2) and tInit2(2)

[g1L{2,2}, g2L{2,2}, g3Local] = sdlfmXsdlfmvKernGradientBlock(lfmKern1, lfmKern2, tInit1(2), ...
    tInit2(2), kyy, kyv, kvy, kvv, i, j, generalConst, generalConstGrad, 1);

g1Plus = sdlfmvXsdlfmvKernComputeBlock(lfmKern1, lfmKern2, tInit1(2), ...
    tInit2(2), kyy, kyv, kvy, kvv, i, j, generalConst);

g2Plus = sdlfmaXsdlfmKernComputeBlock(lfmKern2, lfmKern1, tInit2(2), ...
    tInit1(2), kyy, kvy, kyv, kvv, j, i, generalConst');

if i~=1 || j~=1
    [gkyy1IC, gkyy2IC, gkyy3IC, gkyy4IC] = sdlfmXsdlfmKernGradientIC(lfmKern1, lfmKern2, ...
        tInit1(2), tInit2(2), gkyy1{i,j}, gkyy2{i,j}, gkyy3{i,j}, gkyy4{i,j}, gkyv1{i,j}, gkyv2{i,j}, ...
        gkyv3{i,j}, gkyv4{i,j}, gkvy1{i,j}, gkvy2{i,j}, gkvy3{i,j}, gkvy4{i,j}, gkvv1{i,j}, gkvv2{i,j}, ...
        gkvv3{i,j}, gkvv4{i,j}, 1, typeParam);
    if i>=2 && j>=2
        g1L{2,2}{1} = g1L{2,2}{1} + gkyy1IC{1};
        g2L{2,2}{1} = g2L{2,2}{1} + gkyy2IC{1};
        g1L{2,2}{2} = [gkyy1IC{2} g1L{2,2}{2}];
        g2L{2,2}{2} = [gkyy2IC{2} g2L{2,2}{2}];
        g1L{2,2}{3} = [gkyy1IC{3} g1L{2,2}{3}];
        g2L{2,2}{3} = [gkyy2IC{3} g2L{2,2}{3}];
    else
        for k=1:3
            g1L{2,2}{k} = g1L{2,2}{k} + gkyy1IC{k};
            g2L{2,2}{k} = g2L{2,2}{k} + gkyy2IC{k};
        end
    end
    cLength = length(g3Local);
    g3L{2,2} = zeros(1, cLength+1);
    g3L{2,2}(1:cLength) = g3L{2,2}(1:cLength) + g3Local;
    cLength = length(gkyy3IC);
    g3L{2,2}(1:cLength) = g3L{2,2}(1:cLength) + gkyy3IC;
    g3L{2,2}(i+1) = g3L{2,2}(i+1) + g1Plus;
    g3L{2,2}(j+1) = g3L{2,2}(j+1) + g2Plus;
    g4L{2,2} = gkyy4IC;
else
    g3L{2,2} = g3Local;
    g3L{2,2}(i+1) = g1Plus;
    g3L{2,2}(j+1) = g3L{2,2}(j+1) + g2Plus;
    g4L{2,2} = sdlfmKernMeanCovPartial(lfmKern1(1), lfmKern2(1), tInit1(2), ...
        tInit2(2), 1, typeParam);
end

