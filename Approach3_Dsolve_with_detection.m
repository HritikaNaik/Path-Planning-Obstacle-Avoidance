%% This program is concerned with the problem of planning a 2D motion for a robot from point (0, 0) to point (5, 2).
% The robot might have to avoid obstacles, 
% The controls for the robot is its accelaration u
% The optimal control problem is of the form:
% $$J = \int_{0,0}^{5,2} \frac{1}{2}*$$
% subject to:  $$ \frac{d x}{dt} = v * cos(theta)$$
%              $$ \frac{d y}{dt} = v * sin(theta)$$
% Each time it comes accros an obstacle, it recomputes its trajectory, with repect to the neaarest point of the obstacle.
clc
close all;
%% Create environment
load exampleMaps.mat
map = binaryOccupancyMap(simpleMap,1);

%% Create Obstacle or goal or both as [x y label]
%% Uncomment based on case
% No obstacles
% objects = [5, 2, 1];
% One obstacle
% objects = [0.4,2,2;  5, 2, 1];
% Three obstacles
 objects = [0.4,2,2; 1.4,2.8,2; 4.1,3.7,2; 5, 2, 1];
% compound obstacle
% objects = [0.4,1.8,2; 0.4,1.9,2; 0.3,1.9,2; 0.3,1.8,2; 5, 2, 1];
       
%% Create object Detector sensor
detector = ObjectDetector;
detector.fieldOfView = 3*pi/4;
detector.maxRange = 0.3;

%% Create visualizer
viz = Visualizer2D;
attachObjectDetector(viz,detector);
viz.objectColors = [1 0 0;0 1 0;0 0 1];
viz.objectMarkers = 'so^';

%% Simulation parameters
sampleTime = 0.05;             % Sample time [s]
initPose = [0; 0; pi/4];        % Initial pose (x y theta)
tVec = 0:sampleTime:2;         % Time array

pose = zeros(3,numel(tVec)); % Pose matrix
pose(:,1) = initPose; 

i = 0;j = 0;
syms x1 x2 X t u sys


%% Simulation loop
r = rateControl(1/sampleTime);
for idx = 1:numel(tVec)   
    time = idx * sampleTime; % Time since start of simulation
    
    % Update object detector and visualization
    detections = detector(pose(:,idx),objects);
    
    viz(pose(:,idx),objects) % Display current position
    ylim([0 4.2]);
    xlim([0 5.1]);
    pause(0.25)
    
    pos = pose(:,idx); %  Current position of robot
     
    if i == 0 %Beginning of program
        sys = Find_path(pos(1),pos(2),time,[0,0],0)
        i = 1;
    end
    
    if (~isempty(detections) && j == 0) % Detecting this object for the first time
        B = size(detections);  % How many objects at once
        obj = [0 0];
        minim = 2; % Assume least distance to obstacle is the range itself
        for id = 1:B
            C = detections(id,:);
            if (C(3) == 2 && detections(1)<minim) % Which is the closest of all detected objects
                obj(1) = pos(1) + C(1)*cos(pos(3)); % Position of obstacle
                obj(2) = pos(2) + C(1)*sin(pos(3));
                a_p = pos % Save current value
                a_o = obj
                minim = detections(1);
                pos(1) = pos(1) + C(1)*cos(C(2)+(pi/6))/2;
                if C(2) > 0 % Obstacle is at a positive angle
                    pos(2) = pos(2) + C(1)*sin(C(2)+(pi/6))/2;
                else
                    pos(2) = pos(2) + C(1)*sin(C(2)-(pi/6))/2;
                end
            end
        end
        if (minim < 2) % Closest obstacle
            sys = Find_path(pos(1),pos(2),time,obj,1) % Generate new trajectory
        end
        j = 1;
    end
    
    if (isempty(detections) && j == 1) % Once away from obstacle, original values returned
        j = 0;
        a = pos;
    end
    
    %% Computing New Positions 
    x = subs(sys(1),time);
    y = subs(sys(2),time);
    
    
    if (x-pos(1)) == 0    
        ang  = pi/2;
    else
        ang = atan((y-pos(2))/(x-pos(1)));
    end
    pose(:,idx+1) = [x,y,ang].'; % Change in position
    %end
    
    if (abs(5 - x) < 0.05 && abs(2 - y) < 0.05)  % When its close enough to the goal
        break
    end
    
    waitfor(r);
end

viz(pose(:,idx+1),objects)
ylim([0 4.2]);
xlim([0 5.1]);

function equations = Find_path(a,b,ti,obj,k)
    % State equations
    syms x y p1 p2 u;
    Dx = y;
    Dy = -x + u;
    
    % Cost function inside the integral
    syms g;
    g = 0.5*u^2 - k*(1 - ((x-obj(1))/0.05 + (y-obj(2))/0.05));
    
    % Hamiltonian
    syms p1 p2 H;
    H = g + p1*Dx + p2*Dy;
    
    % Costate equations
    Dp1 = -diff(H,x);
    Dp2 = -diff(H,y);
    
    % Solve for control u
    du = diff(H,u);
    sol_u = solve(du, u);
    
    % Substitute u to state equations
    Dy = subs(Dy, u, sol_u);
    
    % Convert symbolic objects to strings for using 'dsolve'
    eq1 = strcat('Dx=',char(Dx));
    eq2 = strcat('Dy=',char(Dy));
    eq3 = strcat('Dp1=',char(Dp1));
    eq4 = strcat('Dp2=',char(Dp2));
    
    % Use boundary conditions to determine the coefficients
    %    Original Case: (a) x(0)=y(0)=0; x(2) = 5; y(2) = 2;
    con1 = 'x('+string(ti)+') = '+string(a);  % 'x(0) = 0';
    con2 = 'y('+string(ti)+') = '+string(b);  % 'y(0) = 0';
    con3 = 'x(2) = 5';
    con4 = 'y(2) = 2';
    
    sol = dsolve(eq1,eq2,eq3,eq4,con1,con2,con3,con4);% Solution

    syms X Y
    X = sol.x;
    Y = sol.y;
    equations = [X;Y];
end
