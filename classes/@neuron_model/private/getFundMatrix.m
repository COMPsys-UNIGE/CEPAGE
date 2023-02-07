function fundMatrix = getFundMatrix(object,nPoints,limitCycle)

    fundMatrix = cell(nPoints,1);
    auxMatrices = cell(nPoints,1);
    
    T = limitCycle.T;
    X = limitCycle.X;
    
    nx = object.getnx();
    
    period = (T(end)-T(1));
    
    deltaPhi = period/nPoints;
    
    CI = eye(nx);
    xPoint = interp1(T,X,deltaPhi*(0:nPoints-1));
    J = object.getJacobian(zeros(size(xPoint,1),1),xPoint');
    
    for i=1:nPoints
      [~,tmp] = ode45(@funForODE,[0 deltaPhi],CI(:),odeset,J{i});
      auxMatrices{i} = reshape(tmp(end,:),nx,nx)';
    end

    
    parfor i=1:nPoints
        tmp2 = eye(nx);
        for j=[i-1:-1:1 , nPoints:-1:i]
            tmp2 = tmp2*auxMatrices{j};
        end
        fundMatrix{i} = tmp2;
        
    end

    function dx = funForODE(t,x,J)
        nx_ = sqrt(numel(x));
        x = reshape(x,nx_,nx_)';
        dx = (J*x)';
        dx = dx(:);
    end
    
end