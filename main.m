clear;
clc;

%%% Here we use the general algorithm ()
%% INPUT DATA
N_cases = 3;

n_iter = 15;

gamma = 0.5;

phi = [0.01, 0.1, 1];
far_point = [230, 180, 180];
END = [250, 200, 200];
Np = [1000, 1200, 2000];
bb = [30, 20, 30];
eps = 0.00001./bb;
params.phi = phi;
params.END = END;
params.far_point = far_point;
params.Np = Np;
params.bb = bb;
params.eps = eps;
params.gamma=gamma;

% phi = 0.01;       % phi = 0.1;        % phi = 1;          % phi = 10;
% far_point = 230;	% far_point = 180;  % far_point = 180;  % far_point = 130;
% END = 250;        % END = 200;        % END = 200;        % END = 150;
% Np = 1000;        % Np = 1200;        % Np = 2000;        % Np = 8000;
% bb = 30;          % bb = 20;          % bb = 30;          % bb = 10;
% eps = 0.00001/bb;	% eps = 0.00001/bb; % eps = 0.00001/bb; % eps = 0.00001/bb;

%in this loop we define the integration interval X
for j=1:N_cases

    A = END(j);
    N1 = floor(1-(Np(j)-1)*bb(j)*log(eps(j))/(A-bb(j)-bb(j)*log(eps(j))));
    N2 = Np(j)-N1;
    k = (A-bb(j))/N2;
    m = bb(j) - k*N1;
    xp_0 = zeros(1,length(N1));
    for q=1:N1
        xp_0(q) = eps(j)*exp(-log(eps(j))*(q-1)/(N1-1));
    end
    xp_0 = xp_0 - eps(j);
    xp = bb(j).*xp_0;
    log_step = xp(end)-xp(end-1);
    Xp_0 = linspace(xp(end)+log_step,A,N2);
    
    Xp = [xp Xp_0];
    Xm = -fliplr(Xp);
    struct(j).X = [Xm(1:end-1), Xp];

    struct(j).DISCR = zeros(2, n_iter+1);
    struct(j).DIFF = zeros(4, n_iter+1);
    struct(j).SOL = zeros(4, n_iter+1,length(struct(j).X));
end


TIME = [0 0 0];
%this loop runs the algorithm
for j=1:N_cases
    time = cputime;
    [struct(j).DISCR, struct(j).DIFF, struct(j).SOL] = ...
        benchmark(phi(j), gamma, n_iter, far_point(j), END(j), struct(j).X);
    disp(['***************** Case ',num2str(j),' done ***************'])
    TIME(j) = time - cputime;
    
end

params.time=TIME;

% Create folder name string
folderName = sprintf('results/gamma_%g', gamma);

% Check if the folder exists, create if not
if ~exist(folderName, 'dir')
    mkdir(folderName);
end

% Save the file in the specified folder
save(fullfile(folderName, 'results.mat'), 'params', 'struct');



