%% This program is concerned with the problem of planning a 2D motion for a robot from point (0, 0) to point (1.2, 1.6).
%% The robot might have to avoid obstacles, assume circular, with radius = 0.28
%% The controls for the robot are its speed v and its angle (direction) theta
%% The optimal control problem is of the form:
%% J = integral( 2nd derivative of x + 2nd derivative of y)
%% subject to:  derivative of x = v * cos(theta)
%%              derivative of y = v * sin(theta)
%% as well as obstacle constraints 
%% The software Casadi has been used to solve this
clear all
close all
clc


tiledlayout(2,2); %% Showing the approach for four different cases

%% No obstacles
nexttile;
xx = Get_path(0);
plot(xx(1,:),xx(2,:))
xlim([0 1.8])
ylim([0 1.8])

ang=0:0.01:2*pi; 
xp = 0.28 * cos(ang);
yp = 0.28 * sin(ang);

%% One obstacle
nexttile;
xx = Get_path(1);
plot(xx(1,:),xx(2,:))
hold on
xlim([0 1.8])
ylim([0 1.8])
plot(0.4+xp,0.5+yp)
plot(xx(1,:),xx(2,:))
hold off

%% Two obstacles
nexttile;
xx = Get_path(2);
plot(xx(1,:),xx(2,:))
hold on
xlim([0 1.8])
ylim([0 1.8])
plot(0.4+xp,0.5+yp)
plot(0.8+xp,1.5+yp)
plot(xx(1,:),xx(2,:))
hold off

%% Three obstacles
nexttile;
xx = Get_path(3);
plot(xx(1,:),xx(2,:))
hold on
grid on
xlim([0 1.8])
ylim([0 1.8])
plot(0.4+xp,0.5+yp)
plot(0.8+xp,1.5+yp)
plot(1+xp,0.8+yp)
plot(xx(1,:),xx(2,:))
hold off

function xx = Get_path(M)
    % CasADi v3.4.5
    addpath('C:\Users\Hritika\Desktop\Study\SOP\NEW\Casadi\casadi-windows-matlabR2016a-v3.5.5')
    import casadi.*

    T = 0.2; %[s]
    N = 10; % prediction horizon
    
    v_max = 0.3; v_min = -v_max; % range of control variable v
    theta_max = pi/2; theta_min = -theta_max; % range of control variable theta 

    x = SX.sym('x'); y = SX.sym('y'); % define states
    states = [x;y]; 
    n_states = length(states);

    v = SX.sym('v'); theta = SX.sym('theta'); % define control variables
    controls = [v;theta]; n_controls = length(controls);
    rhs = [v*cos(theta);v*sin(theta)]; % system r.h.s

    f = Function('f',{states,controls},{rhs}); % nonlinear mapping function f(x,u)
    U = SX.sym('U',n_controls,N); % Decision variables (controls)
    P = SX.sym('P',n_states + n_states);
    % parameters (which include at the initial state of the robot and the reference state)

    X = SX.sym('X',n_states,(N+1));
    % A vector that represents the states over the optimization problem.

    obj = 0; % Objective function
    g = [];  % constraints vector

    Q = zeros(2,2); Q(1,1) = 2;Q(2,2) = 2.1;% % weighing matrices (states)

    vel_ = zeros(1,2); % variable storing previous velocity for each given state
    st  = X(:,1); % initial state
    g = [g;st-P(1:2)]; % initial condition constraints

    for k = 1:N
        st = X(:,k);  con = U(:,k);
        vel = [con(1)*cos(con(2)) con(1)*sin(con(2))]; % to get [v*cos(theta),v*sin(theta)]
        acc = (vel - vel_); vel_=vel(1); %To calculate the derivative of velocity for cost
        obj = obj+(st-P(3:4))'*Q*(st-P(3:4))+ (acc(1)^2 + acc(2)^2)/0.5;%+ acc'*R*acc %+ con'*R*con; % calculate obj
        st_next = X(:,k+1); %Objective function has additional term to ensure it goes to goal
        f_value = f(st,con);
        st_next_euler = st+ (T*f_value); % compute for next step
        g = [g;st_next-st_next_euler]; % compute constraints
    end
    % Add constraints for collision avoidance
    % M is the number of obstacles
    if M == 1 
        for k = 1:N+1
            g = [g; -sqrt((X(1,k)-0.4)^2+(X(2,k)-0.5)^2) + 0.3];
        end
    elseif M == 2 
          for k = 1:N+1
            g = [g; -sqrt((X(1,k)-0.4)^2+(X(2,k)-0.5)^2) + 0.3];
            g = [g; -sqrt((X(1,k)-0.8)^2+(X(2,k)-1.5)^2) + 0.3];
          end  
    elseif M == 3
        for k = 1:N+1
            g = [g; -sqrt((X(1,k)-0.4)^2+(X(2,k)-0.5)^2) + 0.3];
            g = [g; -sqrt((X(1,k)-0.8)^2+(X(2,k)-1.5)^2) + 0.3];
            g = [g; -sqrt((X(1,k)-1)^2+(X(2,k)-0.8)^2) + 0.3];
          end
    end
    
    xs = [1.2 ; 1.6]; % Reference posture.

    % make the decision variable one column  vector
    OPT_variables = [reshape(X,2*(N+1),1);reshape(U,2*N,1)];

    nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P);

    opts = struct;
    opts.ipopt.max_iter = 100; % maximum number of iterations
    opts.ipopt.print_level =0; % so that output is not too detailed
    opts.print_time = 0;
    opts.ipopt.acceptable_tol =1e-8; % acceptable relative convergence tolerance
    opts.ipopt.acceptable_obj_change_tol = 1e-6; % acceptance stopping criterion based on objective function change.

    solver = nlpsol('solver', 'ipopt', nlp_prob,opts); % NLP solver defined

    args = struct;
    args.lbg(1:2*(N+1)) = 0; % equality constraints
    args.ubg(1:2*(N+1)) = 0; % equality constraints

    args.lbg(2*(N+1)+1 : (N+1)+ (M+1)*(N+1)) = -inf; % inequality constraints
    args.ubg(2*(N+1)+1 : (N+1)+ (M+1)*(N+1)) = 0; % inequality constraints

    args.lbx(1:2:2*(N+1),1) = -2; %state x lower bound
    args.ubx(1:2:2*(N+1),1) = 2; %state x upper bound
    args.lbx(2:2:2*(N+1),1) = -2; %state y lower bound
    args.ubx(2:2:2*(N+1),1) = 2; %state y upper bound

    args.lbx(2*(N+1)+1:2:2*(N+1)+2*N,1) = v_min; %v lower bound
    args.ubx(2*(N+1)+1:2:2*(N+1)+2*N,1) = v_max; %v upper bound
    args.lbx(2*(N+1)+2:2:2*(N+1)+2*N,1) = theta_min; %theta lower bound
    args.ubx(2*(N+1)+2:2:2*(N+1)+2*N,1) = theta_max; %theta upper bound
    %----------------------------------------------
    % ALL OF THE ABOVE IS JUST A PROBLEM SET UP


    % THE SIMULATION LOOP SHOULD START FROM HERE
    %-------------------------------------------
    t0 = 0;
    x0 = [0 ; 0 ];    % initial condition.

    xx(:,1) = x0; % xx contains the history of states
    t(1) = t0; 

    u0 = zeros(N,2);        % two control inputs for each robot
    X0 = repmat(x0,1,N+1)'; % initialization of the states decision variables

    sim_tim = 30; % Maximum simulation time

    % Start MPC
    mpciter = 0;
    xx1 = [];
    u_cl=[];

    % the main simulaton loop... it works as long as the error is greater
    % than 10^-6 and the number of mpc steps is less than its maximum
    % value.
    x_prev = [0 0];
    main_loop = tic;
    while(abs(norm((x0-xs),2)) > 1e-3 && mpciter < sim_tim / T)
        args.p   = [x0;xs]; % set the values of the parameters vector
        % initial value of the optimization variables
        args.x0  = [reshape(X0',[],1);reshape(u0',2*N,1)];
        sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
            'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);
        u = reshape(full(sol.x(2*(N+1)+1:end))',2,N)'; % get controls only from the solution

        xx1(:,1:3,mpciter+1)= reshape(full(sol.x(1:3*(N+1)))',3,N+1)'; % get solution TRAJECTORY
        u_cl= [u_cl ; u(1,:)];
        t(mpciter+1) = t0;

        x_prev = [x0(1) x0(2)];
        % Apply the control and shift the solution
        [t0, x0, u0] = shift(T, t0, x0, u,f);
        xx(:,mpciter+2) = x0;
        X0 = reshape(full(sol.x(1:2*(N+1)))',2,[])'; % get solution TRAJECTORY
        % Shift trajectory to initialize the next step
        X0 = [X0(2:end,:);X0(end,:)];
        mpciter;
        mpciter = mpciter + 1;
    end;

    main_loop_time = toc(main_loop);
    ss_error = norm((x0-xs),2)
    average_mpc_time = main_loop_time/(mpciter+1)
end


function [t0, x0, u0] = shift(T, t0, x0, u,f)
    st = x0;
    con = u(1,:)';
    f_value = f(st,con);
    st = st+ (T*f_value);
    x0 = full(st);

    t0 = t0 + T;
    u0 = [u(2:size(u,1),:); u(size(u,1),:)];
end



