%% M A S T E R  D O C
% Courtesy of: SEBASTIAN ABISLEIMAN

%% Integration in 1D and 2D, Definite and Indefinite
%% Derivations
%% Gradients
%% Logical Indexing, Sort
%% Linear System of Equations
%% Plotting 1D and 2D Functions
%% Plotting surf, contour, mesh, quiver (vector field), histogram, colorbar
%% Random Numbers rand (interval between 0 and 1), randn (normally distributed numbers), randi (integers)
%% Mean, Max, Min, Std, Mod
%% Vectorization (way to express a formula to express not just scalar, but also array) !!!
%% Center of Mass 
%% Vector Manipulation
%% Monte Carlo Estimation
%% Simple Interpolation & Extrapolation
%% Linear Fit
%% Displaying Results
%% For & While Loops
%% If & Switch Bifurcations
%% Fibonacci Sequence
%% Variable Types
%% Changing Characters
%% Cell Structure
%%
%%
%%
%%
%%

%% FORMATING
%  format long; pi 		gives it with lots of decimals
%  format short; pi		gives it with few decimals

%% E TO THE POWER
% exp(1) = e^1

%% STEPS
% x = [0:1:20]
% x = [start:size of step:end]
% x = linspace(0,20,20)
% x = linspace(start, end, how many steps)
% use ?:? to claim all of that row/column

%% VECTORIZATION
% just put a . before the command
% when multiplying a vector times a scalar 

%% TENSOR PRODUCT
%%y*x instead of x*y

%% PLOTTING 1D
%%figure
%%plot(x,y,'o-')
%%xlabel(?label?)
%%legend(?x name?, ?y name?)
%%ylim([0, 10]) sets a range for the y axis
%%grid on

%% SUBPLOT
% subplot(1,4,1)
% plot(etc.)
% allows for multiple plots in same figure
% 1st input, number of rows
% 2nd input, number of columns
% 3rd input, number of graph going to plot

%% CONVERTING TO COLUMN (AKA TRANSPOSE)
% x = [1 2 3]? = [1;2;3]
% apostraphe transposes

%% EXPONENT
% 10e3 = 10*10^3

%% E TO THE POWER
% exp(whatever function)

%% SUM
% sum(x)

%% CUMPRODUCT
% cumprod(x,dimension)
% takes the product by term in the dimension specified
% dimension 1 is down across rows, dimension 2 is across or across columns

%% CUMSUM
% cumsum(x)
% takes the sum by term in the dimension specified
% dimension 1 is down across rows, dimension 2 is across or across columns

%% ONES AND ZEROS
% ones(4,1)	
% zeros(4,1)
% creates matrices of ones and zeros respectively
% first number rows, second number columns

%% TRAPEZOIDAL RULE
% trapz(x,y)
% takes trapezoidal integral of y function with respect to x 
% y and x must be same length
% x can also be a column vector and y an array with same length as x

%% CALCULATING TIME
% tic
% function
% t1 = toc		
% calculates time of function calculation; time btw tic and toc
%% FIRST DERIVATIVE
% dydx = diff(y)./diff(x)		*PLOTS ONE LESS VALUE*
% dydx = gradient(y,x)
%
% *NOTE: Y IS FIRST, X IS SECOND AS INPUTS*

%% GRADIENT IN 3D
% [dzdx,dzdy] = gradient(Z,dx,dy);
% [Gx,Gy] = gradient(Z,dx,dy);
% this creates partial derivatives of Z with respect to dx and dy
% USED TO MAKE QUIVER WORK

%% MESHGRID
% [X,Y]=meshgrid(x,y);
% creates matrices X and Y from vectors x and y
% use for 2D plotting

%% SURF PLOT
% [A,B] = meshgrid(a,b);
% Z = A.^2 + B.^2;
% surf(A,B,Z)
% A, B, and Z are matrices
% makes 2D plot

%% CONTOUR
% contour(Z)
% makes 1D plot with contour lines of a matrix Z
% 
% contour(X,Y,Z)
% must use meshgrid prior (see above)
% makes 2D plot with contour lines 
% plot of matrix Z using X and Y matrices from meshgrid
% contour(X,Y,Z,N)
% N determines number of lines wanted

%% MESH
% mesh(X,Y,Z)
% must use meshgrid prior (see above)
% makes 2D plot like surf but is like a mesh or a fishing net
% plot of matrix Z using X and Y matrices from meshgrid

%% QUIVER/GRADIENT PLOT
% [dzdx,dzdy] = gradient(Z,dx,dy);  need gradient for quiver
% [X,Y] = meshgrid(x,y);            need meshgrid to create matrices X & Y
% quiver(X,Y,dzdx,dzdy)             actual plot
% plots a gradient or a field of vectors
% much like a plot of electric field

%% COLOR BAR
% colorbar
% displays a colorbar

%% BAR 
% bar(y) or bar(x,y) 
% creates bars for each value
%

%% HISTOGRAM
% histogram(x,m)
% plots a histogram of vector or matrix x
% m represents number of bins wanted
% histogram(x,edges)
% edges is a vector descring the edges of the bins

%% INTEGRAL
% int = cumsum(y*dx)		Indefinite integral
% int = sum(y*dx)			Riemann sum
% int = trapz(x,y)          Trapezoidal method integral

%% LOGICAL INDEXING
% Use to find values that meet conditions as well as location of those values
% x>3       gives 0 and 1 as a logical statement
%           1 is true and 0 is false
% x(x>3)    gives values of x that satisfy condition other possible conditions:
%
%  <= (less than or equal to)
%  >= (greater than or equal to)
%  == (equal to)
%  | (or)
%  & (and)

%% SORT
% B = sort(A)           sorts the elements of vector A in ascending order.
% B = sort(A,dim)       returns the sorted elements of matrix A along dimension dim
%                       dim = 1 across; dim = 2 down
% B = sort(___,direction)   sorts elements by direction
%                           direction  = 'ascending' or 'descending'

%% RESIZING PLOT
% Resize = plot (Size*x, Size*y)

%% ROTATING PLOT
% alpha = deg2rad(Angle);
% Rot = [cos(alpha) -sin(alpha); sin(alpha) cos(alpha)];
% Rotate = plot(Rot*x, Rot*y) 
% input Angle want to rotate by

%% TRANSLATING PLOT
% Translate = plot(Distance+x, Distance+y)
% translates the function plotted by Distance added or subtracted

%% LINEAR SYSTEMS OF EQUATIONS
% A= matrix of coefficients from system of equations
% should be organized as each equation coefficients is one row in matrix
% b = right side of equations; what the equations equal
% trying to solve for vertical vector of variable values
%
% Example:
% A = [3 2 -1; -1 3 2; 3 -1 -1];
% b = [10; 5; -1];
%
% Inverse method: 
% inv(A)*b
%
% Gaussian method: 
% A\b 
%
% Cramer's Method (determinants): 
% A1 = [b A(:,2:end)];
% A2 = [A(:,1) b A(:,3)];
% A3 = [A(:,1:2) b];
% X1 = det(A1)/det(A);
% X2 = det(A2)/det(A);
% X3 = det(A3)/det(A);
% x4 = [X1; X2; X3]
%
% Reduced Row Echelean Method: 
% X=rref([A b]); 
% X=X(:,end)

%% RANDOM NUMBERS
% rand(5,1)
% Generates random numbers btw 0 and 1 as a matrix
% 1st input is number of rows, 2nd input is number of columns
%
% randn(5,1)
% Generates random numbers as a matrix
% Numbers are in a NORMAL DISRIBUTION like a gaussian curve
% 1st input is number of rows, 2nd input is number of columns
%
% randi(imax,n)
% Creates matrix of size n with integer numbers from 1 to imax

%% MAX VALUE
% S = max(x)
% x being a vector, will return single largest value of x
% S = max(x,2)
% if x is a matrix use, will return column or row vector depending on 2nd input
% for 2nd input 1 takes max down, 2 takes max across
% [Y,I] = max(X) will return indices for largest values of vector I

%% MIN VALUE
% S = min(x)
% x being a vector, will return single smallest value of x
% S = min(x,2)
% if x is a matrix use, will return column or row vector depending on 2nd input
% for 2nd input 1 takes min down, 2 takes min across
% [Y,I] = max(X) will return indices for min values of vector I

%% MEAN VALUE
% S = std(x)
% x being a vector, will return single value
% S = std(x,2)
% if x is a matrix use, will return column or row vector depending on 2nd input
% for 2nd input 1 takes mean down, 2 takes mean across

%% STANDARD DEVIATION
% S = std(x)
% x being a vector, will return single value
% S = mean(x,0,2)
% if x is a matrix use, will return column or row vector depending on 3rd input
% 2nd input determines weight vector(usually mark as 0 if no weight vector)
% 3rd input 1 takes std down, 2 takes std across

%% MODULUS
% b = mod(a,m)
% returns the remainder after division of a by m
% where a is the dividend and m is the divisor

%% SQAURE ROOT
% sqrt(X) is the square root of the elements of X

%% VECTORIZATION
% just put a . before the command

%% CENTER OF MASS
% x = [b/2 0; -b/2 0; 0 h]';
% m = [1 1 16];
% 
% CenterMass = sum(x.*m,2)/sum(m)
% x being vector of coordinates of points
% m being the weight of the points (if no weight use 1)
% 
% Centering Molecule at the center of mass
% y = x-X;
% SEE WATER CENTER OF MASS PROBLEM FROM NOTES IF NEED MORE INFO

%% VECTOR MANIPULATION
% Defning and plotting vector
% x = [2;5];
% figure
% quiver(0,0,x(1),x(2),0)
% or can plot like this
% plot([X1 X2],[Y1 Y2],[Z1 Z2])
% 
% STRETCHING
% I = eye(2);       Identitty Matrix
% y = 2*I*x;        Stretches by factor of 2
% also
% L = [2 0; 0 3];   
% y = L*x;          Stretches by factor of 2 along x and 3 along y
% 
% REFLECTION
% L = [-1 0; 0 1];  Reflects across y axis by making x values negative
% y = L*x;
% 
% ROTATION
% alpha = deg2rad(90)       Input angle in degrees
% R = [cos(alpha) -sin(alpha); sin(alpha) cos(alpha)]
% y = R*x
% 
% TRANSLATION
% D = [Distance;Distance];
% y = D+x;

%%  MONTE CARLO ESTIMATION 1D
% Defining function to estimate
% R = pi;              Represents square area around portion graph analysing
% x = linspace(0,pi,100);
% dx = pi/100;
% y = sin(x);
% 
% Generating Random Points
% Npoints = 10000;          Represents # of rando #s or aka accuracy  
% X = pi*rand(1,Npoints);   making rando data along the x-axis that is from
%                           0 to pi
% Y = rand(1,Npoints);      making rando data along the y-axis
% 
% Checking if points are in or out of function
% check = Y <= sin(X);      given an x is the y value less than or equal to
%                           the function
% Nwithin = sum(check);     number of points inside function
% 
% Plotting
% plot(X(check),Y(check),'r.')      plots points inside func
% plot(X(~check),Y(~check),'b.')    plots point outside func with ~
% tilde ~ is the *not* operator
% 
% ACTUAL ESTIMATE
% AreaMC = R*Nwithin/Npoints;       
% Remember R is square area around portion graph analysing

%% MONTE CARLO ESTIMATION 2D
% Defining function to estimate
% 
% V = 100;                  total volume in cube encompassing function
% dx = 0.05;
% dy = 0.05;
% x = -5:dx:5;
% y = -5:dy:5;
% [X,Y] = meshgrid(x,y);
% Z = exp(-(X.^2+Y.^2));
%
% Generating Random Points
% Npoints = 100000;         Represents # of rando #s or aka accuracy
% X = 10*rand(1,Npoints)-5; mult and subt to get domain to be from -5 to 5
% Y = 10*rand(1,Npoints)-5; 
% Z = rand(1,Npoints);
%
% Checking if points are in or out of function
% check = Z < exp(-(X.^2+Y.^2));
% Nwithin = sum(check);
% 
% Plotting
% plot3(X(check),Y(check),Z(check),'r.')
% plot3(X(~check),Y(~check),Z(~check),'b.') % tilde ~ is not operator
% 
% ACTUAL ESTIMATE
% VolumeMC = V*Nwithin/Npoints;
% Remember V is square area around portion graph analysing

%% SIMPLE INTERPOLATION
% Defining function
% x = 0:6;
% y = sin(x);
%
% xq = 0:0.1:2*pi;                  represents query points
% yq1 = interp1(x,y,xq,'linear');   more accurate; linear approximation
% yq2 = interp1(x,y,xq,'spline');   less accurate; spline approximation

%% SIMPLE EXTRAPOLATION
% Estimating R by minimizing the quadratic error between I^2 and p
% can extrapolate further points and interpolate missing internal points
% given vectors I and P that only have 10 values
% R3 = (I.^2)'\P';          we model a system of equations
% X = 0:13;                 we go further than given points
% plot(X,X.^2*R3,'r-');     we evaluate for up to 13 points
%                           Go to JoulesLaw Notes for more info

%% LINEAR FIT
% 
% Plotting Original Data
% figure(1)
% scatter(x,y)              x and y vectors
% Data fitting
% p = polyfit(x,y,1);       1 determines order function want to fit
%                           creates polynomial
% Overlaying fitted line
% hold on
% X = 0:5;
% Y = polyval(p,X);         evaluates polynomial at X value
% plot(X,Y,'r-');
% Mean Quadratic Error
% a = p(1);
% b = p(2);
% mqe = mean((y - (a*x+b)).^2);

%% DISPLAYING RESULTS
% fprintf('My result is %1.4f\n', AreaMC)  
% 
% prints text in command line
% \n means new line
% 1.4 means 1 decimal place before decimal point and 4 decimal places after the decimal point
% fprintf('Relative Error  : %f%%\n', AreaMC) 
% two %% mean a percentage symbol when printing
% %f is expecting the result of a function
% %s is expecting a string
% Can also use disp to display text
% disp('Done!')

%% FOR LOOP
% must declare variables working with before loop (can set empty vector or variable to zero)
%
% s = 0;
% for n = 1:1:10         executes 10 steps between 1 and 10 in steps of 1
%    s = s + n;          declares a new s value for each iteration
%                        that is called on again in the next iteration
% end

%% WHILE LOOP
% executes the enclosing command 'while' a specified condition is true
% 
% count = 1;
% check = 1;
% while check == 1                  while this condition is true                  
%     x = randi(10);
%     if mod(x,2) == 1              checking for odd numbers bc number /2
%         count = count + 1;        if results in remainder of 1 is odd
%     else 
%         check = 0;                if not an odd number the while conditio
%     end                           gets set to false and the loop ends
% end
% 
% fprintf('First Even Number : %i\n', x)
% 
% fprintf('Number of iterations until even number is drawn : %i\n', count)

%% IF BIFURCATION
% 'if' a condition is true then execute the enclosing code
% 'elseif' if the previous condition was false and this is true execute the
%          enclosing code
% 'else' if all previous conditions are not met execute code
% 
% #Example#
% n = input('Input a number : ');
%
% if mod(n,2) == 0                      %is the number inputted even?
%     
%     fprintf('Number %i is even\n',n)  %then
%     
% else fprintf('Number %i is odd\n',n)  %otherwise
%     
% end

%% SWITCH BIFURCATION
% 
% switch b              the variable b is considered
%     case 1            in the case that b is 1 execute code
%         tot = 0.75;
%     case 2            in the case that b is 2 execute code
%         tot = 1.25;
%     case 3            in the case that b is 3 execute code
%         tot = 1.65;
%     otherwise         in any other case execute code
%         tot = 1.65 + 0.3*(b-3);
% end

%% SOUND
% soundsc(x,fs)     x is vector of audio signal/frequencies
%                   fs is sample rate
% soundsc(x)        gives a default sample rate of 8192
% scales the values of audio signal y to fit in the range from ?1.0 to 1.0
% 
% #Example#
% fs = 8192;   % samples per second
% d = 0.5;     % seconds
% 
% x = randn(1,d*fs);            Make white noise and play it
% t = linspace(0,d,d*fs);
% 
% figure
% plot(t,x)
% xlabel('Time')
% ylabel('Amplitude')
% 
% soundsc(x,fs)                 
%
% Creating a Percussion Rhythm
% y = exp(-t./0.05);            Decaying Envelope
% p = y.*x;
% 
% P = [p 0.5*p];
% P = [P P P P];
% 
% figure
% plot(P)
% 
% soundsc(P,fs);

%% FIBONACCI SEQUENCE
% 
% N = input('What number in the Fibonacci sequence ? ');
% 
% Pre-allocating Vector
% 
% X = zeros(1,N);
% ang = zeros(1,N);
% 
% Initital Values
% 
% X(1) = 1;
% X(2) = 2;
% 
% Computing Sequence
% 
% for n = 1:N-2
%     X(n+2) = X(n) + X(n+1);
% end
% 
% theta = (1:N)*pi/2;
% 
% figure
% polar(theta,X)
% 
% hold on
% 
% itheta = linspace (0,max(theta),10*N);
% Y = interp1(theta,X,itheta,'pchip','extrap');
% polar(itheta,Y)
% 

%% VARIABLE TYPES
% dec2bin           use to get binary numbers for letters that are in decimal system
% 
% a = 'hello'       apostraphes indicate that hello is of character variable type
%                    
% double(a)         goes from letters to decimal
% dec2bin(a)        goes from decimal to binary
% 
% Sin * Mantissa * 2^(E - a)   need S, M, and a
% 
% Precision Types
% single : 32bit
% double : 64bit
% 
% 1 byte = 8bits
% 
% string(a)         makes a string of characters to become 1 unit
%
% struc('word')     string vertical catonation 
%
% b = logical(a)    can force matlab to only use 1 bit precision
% 
% double, single, logical, char     different variable types

%% CHANGING CHARACTERS
% s = input('Enter text in small caps : ','s');
%
% S = char(s - 32);
% 
% fprintf('Here is the transformed text : %s\n',S)    % the %s says the function should expect a string

%% CELLS/STRUCTURE
%
% for making a matrix with words
% cell(1,2)         makes cell 1 row and 2 column long
% cellplot(c)       plots a visual of the cell
% 
% Year = [book.year]        returns the values as a vector
% [Year,ind] = sort(Year,'descending')
% Author = [char(book.author)]
%
% Initializing Structure
% book  = struct('author', '','title','','year',zeros(1),'publisher','');     %creates structure
% 
% Entering Book Information
% answer = 'y';
% count = 0;      %book count
% 
% while answer == 'y'
%     count = count + 1;
%     fprintf('Book number: %i\n',count)
%    
%     book(count).author = input('Enter the author : ','s');
%     book(count).title = input('Enter the title : ','s');
%     book(count).year = input('Enter the year : ');
%     book(count).publisher = input('Enter the publisher : ','s');
%     
%     answer = input('Do you want to enter another book (y/n)?','s');
%     clc
%     
% end
% 
% fprintf('You have entered %i books\n',count)

%% IMAGE COUNTING
%
% I = imread('rice.png');
% 
% Plotting Image
% figure
% imshow(I)
% 
% Intensity Histogram
% figure
% histogram(I(:))
% xlabel('Intensity')
% ylabel('Count')
% 
% Determine Threshold of Intesity
% threshold = 130;
% 
% Binary Image
% U = I>threshold;
% 
% Plotting Binary Image
% figure 
% imshow(U);
% axis on
% 
% Size of representative grain
% S = 10*20; % in pixels
% 
% Grain Count
% N = sum(U(:))/S;

%% SYMBOLIC CALCULUS
% Defining symbolic function
% syms x
% y = sin(x);
% 
% Plotting symbolic function 
% fplot(x,y,[0 2*pi])
% 
% Indefinite Integral
% int(y,x)      of y by x
% 
% Definite Integral between 0 and pi
% int(y,x,0,pi)
% 
% First Order Derivative 
% diff(y,1,x)
% 
% Second Order Derivative
% diff(y,2,x)
% 
% Taylor expansion around x = 0
% T3 = taylor(y,x,0,'order',3)
% T5 = taylor(y,x,0,'order',5)
% 
% Plotting function and its Taylor Approximation
% figure
% fplot(x,y,[-pi pi]);
% hold on
% fplot(x,T3,[-pi pi]);
% fplot(x,T5,[-pi pi]);
% legend('Original function','3rd Order Approximation','5th Order Approximation',...
%     'location','northwest')
%
%% SYMBOLIC LINEAR SYSTEM OF EQUATIONS
% %% Defining Variables

A = [5 6 -3; 3 -3 -2; 2 -4 -12];
b = [10;14;24];


x = A\b


% Solving the linear system symbolically

As = sym([5 6 -3; 3 -3 -2; 2 -4 -12]);
bs = sym([10;14;24]);
xs = sym('xs',[3 1])


% First solution to symbolic system

xs1 = As\bs

% Second solution to symbolic system

S = solve(As*xs==bs)

xs2 = [S.xs1; S.xs2; S.xs3]

% Conversion of symbolic solution to double precision

double(xs2)

% Conversion of double precision solution to symbols

sym(x)

% Longer way of solving the system symbolically

syms x y z   % defining symbolic variables

eq1 = 5*x + 6*y - 3*z == 10;
eq2 = 3*x - 3*y - 2*z == 14;
eq3 = 2*x - 4*y - 12*z == 24;

% S = solve(eq1,eq2,eq3,x,y,z)
%
%xs3 = [S.x; S.y; S.z]
%
%
%
%%

%% NEED MORE HELP?
% help functionname
% doc functionname

