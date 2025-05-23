function dydt = t3bpdydt(t,y,sys)
% implementation of ode
% sys holds system parameters
% get coordinates
for body = 1:3
    u(body,1:3) = transpose(y(3*(body-1)+1:3*body,1));
end

r = zeros(3,3,3);
e = zeros(3,3,3);

% compose vectors, unit vectors and distances
r(1,2,:) = u(2,:)-u(1,:); % this from m[1] to m[2]
R(1,2)   = sqrt(sum(r(1,2,:).^2)); % magnitude of vector
e(1,2,:) = r(1,2,:)/R(1,2);        % unit vector length
r(2,1,:) = -r(1,2,:);     % .. and back
e(2,1,:) = -e(1,2,:);
R(2,1)   = R(1,2);

r(1,3,:) = u(3,:)-u(1,:); % this from m[1] to m[2]
R(1,3)   = sqrt(sum(r(1,3,:).^2));
e(1,2,:) = r(1,2,:)/R(1,2);        % unit vector length
r(3,1,:) = -r(1,3,:);     % .. and back
e(3,1,:) = -e(1,3,:);
R(3,1)   = R(1,3);

r(2,3,:) = u(3,:)-u(2,:); % this from m[1] to m[2]
R(2,3)   = sqrt(sum(r(2,3,:).^2));
e(2,3,:) = r(2,3,:)/R(2,3);        % unit vector length
r(3,2,:) = -r(2,3,:);     % .. and back
e(3,2,:) = -e(2,3,:);
R(3,2)   = R(2,3);

%% accelerations
% mass 1
dydt(9+1:9+3,1)=sys.theta(2)/R(1,2)*e(1,2,:)+sys.theta(3)/R(1,3)*e(1,3,:);
% mass 2
dydt(9+4:9+6,1)=sys.theta(1)/R(2,1)*e(2,1,:)+sys.theta(3)/R(2,3)*e(2,3,:);
% mass 3
dydt(9+7:9+9,1)=sys.theta(1)/R(3,1)*e(3,1,:)+sys.theta(2)/R(3,2)*e(3,2,:);

%% velocities
dydt(1:9,1) = y(10:18,1);

%%
waitbar(t / sys.tEnd);

end

