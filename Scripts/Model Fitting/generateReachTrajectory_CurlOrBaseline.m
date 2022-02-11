function SOL = generateReachTrajectory_CurlOrBaseline(params)

x0 = params.x0;
trialLength = params.trialLength;

A = eye(size(params.Asys)) + params.Asys*params.dt + params.curlGain*params.Fsys*params.dt;
Aest = eye(size(params.Asys)) + params.Asys*params.dt + params.estGain*params.Fsys*params.dt;
B = params.B;
Q = params.Q;
R = params.R;
Phi = params.Phi;
 
%Terminal Conditions
S = zeros(size(params.Asys,1),size(params.Asys,2),trialLength);
S(:,:,trialLength) = Phi;
L = zeros(size(params.B,2),size(params.Asys,1),trialLength-1);

% Iterate Backwards
for k = trialLength-1:-1:1
    
    % Ricatti Equation
    S(:,:,k) = Q + Aest'*S(:,:,k+1)*Aest - ...
        (Aest'*S(:,:,k+1)*B)*inv(R + B'*S(:,:,k+1)*B)*(B'*S(:,:,k+1)*Aest);
    % Feedback Gain Matrix
    L(:,:,k) = inv(B'*S(:,:,k+1)*B + R)*B'*S(:,:,k+1)*Aest;

end

% Iterate Forward, Find Trajectory
x = zeros(trialLength,size(params.Asys,1));
x(1,:) = x0';
u = zeros(trialLength-1,2);
P = zeros(size(params.Asys,1),size(params.Asys,2),trialLength);

for k = 2:trialLength
    
    u(k-1,:) = -L(:,:,k-1)*x(k-1,:)';

    x(k,:) = ( A*x(k-1,:)' + B*u(k-1,:)' )';
    
    P(:,:,k) = A*P(:,:,k-1)*A';
    
end


SOL.x = x;
SOL.u = u;
SOL.time = [0:params.dt:(trialLength-1)*params.dt]';


