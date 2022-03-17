%*************************************************************************%
%                                                                         %
%                   (Grid) Convergence CFD Verification                   %
%                                                                         %
%*************************************************************************%

%*************************************************************************%
%                                                                         %
% Author:                                                                 %
% Llorente-Lázaro, Víctor Javier  ::  ETSIAE - UPM                        %
%                                                                         %
% Date:                                                                   %
% February 18, 2022                                                       %
%                                                                         %
% Version 1.0:                                                            %
% - Constant refinement ratio                                             %
%                                                                         %
% References:                                                             %
% [1] Roache, P.J., Verification and Validation in Computational Science  %
%     and Engineering, Hermosa Publishers, Albuquerque, New Mexico (1998) %
% [2] Cadafalch, J., Pérez-Segarra, C.D., Cònsul, R., and Oliva, A.       %
%     (2002) Verification of Finite Volume Computation on Steady State    %
%     Fluid Flow and Heat Transfer. Journal of Fluid Engineering          %
%     124(1):11-21                                                        %
% [3] https://www.grc.nasa.gov/www/wind/valid/tutorial/spatconv.html      %
%                                                                         %
%*************************************************************************%

%*************************************************************************%
%                                                                         %
% n: data set (>=3)                                                       %
% h = [h(1) ... h(n)]: geometric discretization parameter representative  % 
%                      of the grid spacing.                               %
%                      h(1) is the finer grid                             %
%                      h(n) is the coarser grid                           %
% phi = [phi(1) ... phi(n)]: numerical solution.                          %
% p = [p(1) ... p(n-2)]: order of accuracy of the numerical scheme.       %
% r = [r(1) ... r(n-1)]: refinement ratio                                 %
% GCI = [GCI(1) ... GCI(n-1)]: Grid Convergence Index (error estimation)  % 
%                                                                         %
%*************************************************************************%

clear all
close all
clc

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             Post-processing                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Case: examples or problem
problem = 'example 2';

% Define grid spacings and some observed quatity corresponging 
% to each grid spacing
switch problem
    case 'example 1'
%       Solutions 
        phi(1) = 0.971023;
        phi(2) = 0.970500;
        phi(3) = 0.968540;        
        phi(4) = 0.961780;
%       Maximum solutions
        phi_max(1) = 1;
        phi_max(2) = 1;
        phi_max(3) = 1;
        phi_max(4) = 1;
%       Grids    
        h(1) = 0.5;
        h(2) = 1;
        h(3) = 2;
        h(4) = 4;   
    case 'example 2'
%       Solutions 
        phi(1) = 0.701487;
        phi(2) = 0.701509;
        phi(3) = 0.701598;        
        phi(4) = 0.701933;
        phi(5) = 0.702980;
        phi(6) = 0.703646;
%       Maximum solutions
        phi_max(1) = 1;
        phi_max(2) = 1;
        phi_max(3) = 1;
        phi_max(4) = 1;
        phi_max(5) = 1;
        phi_max(6) = 1;
%       Grids    
        h(1) = 3.125e-3;
        h(2) = 6.25e-3;
        h(3) = 0.0125;
        h(4) = 0.025;
        h(5) = 0.05;
        h(6) = 0.1;
    case 'example 3'
%       Solutions 
        phi(1) = 0.715007185936;
        phi(2) = 0.715007007122;
        phi(3) = 0.715006709099;        
        phi(4) = 0.714961528778;
        phi(5) = 0.714575171471;
%       Maximum solutions
        phi_max(1) = 1;
        phi_max(2) = 1;
        phi_max(3) = 1;
        phi_max(4) = 1;
        phi_max(5) = 1;
%       Grids    
        h(1) = 6.25e-3;
        h(2) = 0.0125;
        h(3) = 0.025;
        h(4) = 0.05;
        h(5) = 0.1;  
    case 'problem'
    otherwise
        warning('No problem introduce')
end

% Dimensionality of the problem: one, two, three
dimension = 'two';
switch dimension
    case 'one'
        d = 1;
    case 'two'
        d = 2;
    case 'three'
        d = 3;
    otherwise
        warning('No dimension introduce')
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Collection of data                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of data sets
n = size(phi,2);
% Check data sets
if (n < 3)
    warning('The number of data must be three or more')
    return
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               h-refinement                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Refinement ratio 
for i = 1:n-1
    r(i) = h(i+1)/h(i);
end
% Small parameter
epsilon = 1.0e-5;
% Average refinement ratio
r_avg = mean(r);
% Is refinement ratio constant?
if (abs(r(1:n-1) - r_avg) > epsilon)
    warning('Grids with a non-constant refinement ratio')
    fprintf(r)
    return
end
% h should be ordering from the finer grid to the coarser grid (r > 1)
if (r <= 1)
    warning('Refinement ratio should be great than 1')
    fprintf(r)
    return
end
% Control volume of the grid
ConVol = h.^(1/d);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Nodal Classification                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Normalization of the solution
for i = 1:n    
    norm_phi(i) = phi(i)/abs(phi_max(i));   
end
% Difference in the normalized solution
for i = 1:n-1
    D_norm_phi(i) = norm_phi(i) - norm_phi(i+1);
end
% Small parameter
epsilon = 1.0e-30;
% Class of node
% Richardson node: nodes where e = C0*h^p can be calculated. These
%                  Richardson nodes do not necessarily fulfill all the 
%                  requirements for the generalized Richardson 
%                  extrapolation. In fact, the solution may be outside 
%                  the asymptotic range (h not small enough).
% Converged node:  the condition of converged nodes is ill-defined because 
%                  it can also be accomplished by inflection points in the 
%                  solution where all three solutions cross through the 
%                  same point, which are obviously not converged nodes.
% Oscillatory node: zig-zag solution
for i = 1:n-2
    aux = D_norm_phi(i)*D_norm_phi(i+1);
    if (aux > epsilon)
        node{i} = 'Richardson';
    elseif (abs(aux) < epsilon) && (aux > -epsilon)
        node{i} = 'Converged';
    else
        node{i} = 'Oscillatory';
    end
end
% Location of the Richardson node in the arrays
i_R = find(strcmp(node,'Richardson'));
% Location of the Converged node in the arrays
i_C = find(strcmp(node,'Converged'));
% Location of the Oscillatory node in the arrays
i_O = find(strcmp(node,'Oscillatory'));
% Location of the Richardson nodes and Converged nodes
i_RC = [i_R i_C];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Calculation of the observed p                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Diference of the solution
for i = 1:n-1
    D_phi(i) = phi(i+1) - phi(i);
end
% Local order of convergence
for i = 1:n-2
    p(i) = log(D_phi(i+1)/D_phi(i))/log(r_avg);
end
% Global order of convergence at the Richardson nodes
p_avg = mean(p(i_R));
% Standard deviation of p from the global p: measure of how close the 
% solutions are to the asymptotic range
p_std = std(p(i_R));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Calculation of the observed GCI                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Factor of Safety 
if (n >= 3)
    Fs = 1.25;
else
    Fs = 3;
end
% GCI computation
for i = 1:n-2
%   Theoretical absolute error    
    e_abs = abs(D_phi(i));
%   Theoretical relative error    
    e_rel = e_abs/abs(phi(i));
%   Richardson correction
    factor = 1/abs(r_avg^p(i) - 1);
%   Grid Convergence Index    
    GCI(i) = Fs*e_rel*factor;
end
% At the converged nodes the GCI is assumed to be 0
GCI(i_C) = 0;
% Overall volume occupied by all the Richardson nodes and converged nodes
ConVol_total = sum(ConVol(i_RC));
% Global GCI
GCI_avg = dot(ConVol(i_RC),GCI(i_RC))/ConVol_total;
% Standard deviation of GCI
GCI_std = std(GCI(i_RC));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            Asymptotic range                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Checking for asymptotic range
for i = 1:n-3
    ratio(i) = r_avg^p(i)*(GCI(i)/GCI(i+1));
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                Plotting                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Observation vs grid size
subplot(2,2,1)
semilogx(h,phi,'-o')
xlabel('Grid size')
ylabel('Observation')
xlim([h(1) h(n)])
grid on
% Order of convergence vs grid size
subplot(2,2,2)
semilogx(h(1:n-2),p,'-o')
hold on
x = linspace(h(1),h(n),21);
y = p_avg*ones(size(x));
y_1 = (p_avg + p_std)*ones(size(x));
y_2 = (p_avg - p_std)*ones(size(x));
semilogx(x,y,'r-.',x,y_1,'k-.',x,y_2,'k-.')
xlabel('Grid size')
ylabel('Order of convergence')
xlim([h(1) h(n)])
grid on
% GCI vs grid size
subplot(2,2,3)
semilogx(h(1:n-2),GCI*100,'-o')
hold on
y = GCI_avg*100*ones(size(x));
y_1 = (GCI_avg + GCI_std)*100*ones(size(x));
y_2 = (GCI_avg - GCI_std)*100*ones(size(x));
semilogx(x,y,'r-.',x,y_1,'k-.',x,y_2,'k-.')
xlabel('Grid size')
ylabel('GCI(%)')
xlim([h(1) h(n)])
grid on
% Asymptotic range vs grid size
subplot(2,2,4)
semilogx(h(1:n-3),ratio,'-o')
hold on
y = 1*ones(size(x));
semilogx(x,y,'r-.');
xlabel('Grid size')
ylabel('Asymptotic ratio')
xlim([h(1) h(n)])
grid on
%
set(gcf,'color','w')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Array dimension fitting                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% refinement ratio
r(n) = 0;
% node classification 
node{n-1} = '          ';
node{n  } = '          ';
% order of convergence
p(n-1:n) = 0;
% GCI index
GCI(n-1:n) = 0;
% asymptotic range
ratio(n-2:n) = 0;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              Print results                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('--- Performs CFD verification --- \n')
fprintf('\n')
fprintf('Number of data sets = %d\n',n)
fprintf('\n')
fprintf('Grid        Refinement    Solution    Order of       Class of       GCI(%%)      Asymptotic \n')
fprintf('size        ratio                     convergence    node                       ratio (~1)  \n')
fprintf('--------    ----------    --------    -----------    -----------    --------    ----------  \n')
for i = 1:n
    fprintf('%f    %f      %f    %f       %s     %f    %f\n',h(i),r(i),phi(i),p(i),node{i},GCI(i)*100,ratio(i))
end
fprintf('\n')
fprintf('Main results:\n')
fprintf('+ The global order of convergence is %f (std = %f)\n',p_avg,p_std)
fprintf('+ The global GCI is %f%% (std = %f%%)\n',GCI_avg*100,GCI_std*100)
fprintf('\n')
fprintf('--- End --- \n')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END PROGRAM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%