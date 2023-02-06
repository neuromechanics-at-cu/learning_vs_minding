%% Figure 1: Simple Toy Model
% Demonstrate that equivalent trajectories can be produced by varying terms
% q,r, and alpha


% Manipulate cost ratios by changing value of q
q1 = [0.5,1,2];

% Define constants for Dynamical Model
a1 = 1;
r1 = 1;
r2 = 0.01:0.01:5; % Sanity check 
r2 = r1*r2;
b = 1;
t = 0:.2:4;
a2 = 0:0.01:1.2; % Levels of learning

% Intialize matrices
q2 = zeros(length(r2),length(a2),length(q1));
q2prop = zeros(length(r2),length(a2),length(q1));
QRratio = zeros(length(r2),length(a2),length(q1));
Jq = zeros(length(r2),length(a2),length(q1));
Jr = zeros(length(r2),length(a2),length(q1));
Jtot = zeros(length(r2),length(a2),length(q1));

leglabels = {num2str(q1(1)/r1),num2str(q1(2)/r1),num2str(q1(3)/r1)};

for n = 1:length(q1)
    % Describe/Solve Dynamical Model
    C = sqrt(a1^2 + b^2*q1(n)/r1);
    p1 = a1*r1/b^2 + r1/b^2*sqrt(a1^2 + b^2*q1(n)/r1);
    k1 = b*p1/r1;
    x1 = 10*exp((a1-b*k1)*t);

    % Plot Cost Ratios
    for ii = 1:length(a2)
        for jj = 1:length(r2)
            q2(jj,ii,n) = ((a1 + C -a2(ii))^2 - a2(ii)^2)*r2(jj)/b^2;
        end
        QRratio(:,ii,n) = q2(:,ii,n)'./r2;
    end
    
    %Calculate total costs/ratio of costs
    for ii = 1:length(a2)
        for jj = 1:length(r2)
            Jq(jj,ii,n) = sum(q2(jj,ii,n)*x1.^2);
            Jr(jj,ii,n) = sum(r2(jj)*t.^2);
            Jtot(jj,ii,n) = Jq(jj,ii,n) + Jr(jj,ii,n);
        end
    end
end

%Plot Results
figure(101)
plot(a2,squeeze(QRratio(1,:,:)))
xlabel('Proportion learned')
ylabel('Ratio of tracking cost weight to effort cost weight (Q/R)')
xlim([0 1]) % For paper, only plot to 100% learning
legend(leglabels)

figure(102)
plot(a2,squeeze(Jq(1,:,:)./Jr(1,:,:)))
xlabel('Proportion learned')
ylabel('Ratio of tracking cost to effort cost (J_Q/J_R)')
xlim([0 1]) % For paper, only plot to 100% learning
legend(leglabels)
