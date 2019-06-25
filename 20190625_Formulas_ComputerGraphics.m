%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Formulas - Computer Graphics	%
% Gal Argov Sofer		%
% 25/06/2019			%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% World Coordinates
0<=V<=1

% Graphic Window
setviewport(x1, y1, x2, y2);

% Raster
putpixel(x1, y1, color);

% Draw Line
line(x1, y1, x2, y2, color);

% Point
(x1, y1) (x2, y2);

% DDA
range = max(|x2-x1|, |y2-y1|);
x=x1;
dx=|x2-x1|/range;
dy=|y2-y1|/range;

for(i=range;i>0;i--) {
	putpixel(Round(x), Round(y));
	x=x+dx;
	y=y+dy;
}

% Bresenheim
if (|x2-x1| > |y2-y1|) {
	err = |y2-y1| / |x2-x1| -0.5;
	xp=Round(x1);
	yp=Round(y1);
}

putpixel(xp,yp);
if (err > 0)
{
	yp = yp + 1;
	err = err -1;
}
xp = xp + 1;
err = err + |y2-y1| / |x2-x1|;

if (|x2-x1| > |y2-y1|) {
	err = |y2-y1| / |x2-x1| -0.5;
	xp=x1;
	yp=y1;
}

putpixel(xp, yp);
if (err > 0) {
	yp = yp + 1;
	err = err -1;
}
xp = xp + 1;
err = err + |y2-y1| / |x2-x1|;

err = |y2-y1| / |x2-x1| -0.5;
err = err + |y2-y1| / |x2-x1|;
errp = 2 * |x2-x1| * err;

errp = 2 * |y2-y1| - |x2-x1|;
if (errp > 0) {
	yp = yp + 1;
	errp = err -2 * |x2-x1|;
}

% Linear Function
y = mx + b;
y = ((y2-y1)/(x2-x1))*(x-x1)+y1;
Ax+Bx+C=0;
b=-C/B;
m=-A/B;
p={x,y,1}; 																	% vector row
C=[A B C]^T; 																% vector col
{x,y,1} * [A B C]^T = 0;

% P2, P1, L are vectors
P2=P1_tL;
(0<=t<=1);

x = x(P1) + t * x(L);
y = y(P1) + t * y(L);

% P and C is vector
C = [A B C]^T; 																	% Parameters
P = {x, y, 1}; 																	% Point
P * C < 0;

y=f(x);
x=f(t); y=f(t);

% Circle
x^2 + y^2 = r^2;
x = r * cos(t);												% trigonometric function parameter 
y = r * sin(t);												% trigonometric function parameter 
(x-x_c)^2 + (y-y_c)^2 = r^2;
y = y_c + sqrt(r^2 - (x-x_c)^2 );
y = y_c - sqrt(r^2 - (x-x_c)^2 );

% Parabola
y = x^2 / 4 * p;
x = 2 * p * t;
y = p * t^2;

% Polynom
y = a * x^3 + b * x^2 + c * x + d;

% Recursive calculation (Circles) - Calculate sin(dTheta) cos(dTheta) just ones.
x = x_c + r * cos(theta);
y = y_c + r * sin(theta);
x_(n+1)= x_c + r * cos(theta + dTheta);										% dTheta is the step
x_(n+1)= x_c + r * cos(theta) * cos(dTheta) - r * sin(theta) * sin(dTheta);	% dTheta is the step
x_(n+1)= x_c + x_n * cos(dTheta) - y_n * sin(dTheta);						% dTheta is the step
y_(n+1)= y_c + y_n * cos(dTheta) + x_n * sin(dTheta); 						% dTheta is the step

% Bresenhaim
y ^ 2 = r ^ 2 - (x_i + 1) ^ 2;
d_1 = (y_i) ^ 2 - y ^ 2 = (y_i) ^ 2 - r ^ 2 + (x_i + 1) ^ 2;
d_2 = (y) ^ 2 - (y_i + 1) ^ 2 = r ^ 2 - (x_i + 1) ^ 2 - (y_i + 1) ^ 2;
p_i = d_1 - d_2 = 2 * (x_i + 1) ^ 2 + (y_i) ^ 2 + (y_i - 1) ^ 2 - 2 * r ^ 2;
d_1 = (y_i) ^ 2 - r ^ 2 + (x_i + 1) ^ 2;
d_2 = r ^ 2 - (x_i + 1) ^ 2 - (y_i + 1) ^ 2;
p_i = d_1 - d_2 = 2 * (x_i + 1) ^ 2 + (y_i) ^ 2 + (y_i - 1) ^ 2 - 2 * r ^ 2;

if (d_1 < d_2) 																% means p_i < 0
	y = y_i;
else
	y = y_i - 1;

% calculate recursive  parameter p with p_(i + 1) with p_i
p_(i+1) = 2(x_(i+1) + 1) ^ 2 + (y_(i+1) - 1) ^ 2 -2 * r ^ 2;
p_(i+1) = p_i + 4 * x_i + 6 + 2 * (y_(i+1)^2 - (y_i)^2 -2*(y_(i+1) - y_i);
p_i = d_1-d_2 = 2 * (x_i + 1) ^ 2 + y_i ^ 2 + (y_i - 1) ^ 2 -2 * r * 2;
p_1 = 3 - 2 * r;

% Algo
(x_1,y_1)=(0,r);
p_1 = 3 - 2 * r;
for (i = 1; x_i < y_i ; i = i + 1)
{
	x_(i+1)= x_c + x_n * cos(dTheta) - y_n * sin(dTheta);
	y_(i+1)= y_c + y_n * cos(dTheta) + x_n * sin(dTheta);
	x_(i+1) = x_1 +1;
	if (p_i < 0)
		y_i = y_(i+1);
	else if (p_i >= 0)
		y_i = y_(i+1) - 1;

	if (p < 0)
			p_(i+1) = p_i + 4 * x_i + 6;
	if (p >= 0)
		p_(i+1) = p_i + 4 * (x_i -y_i) + 10;
}

%Examples from presentation
Procedure bres_circle (x_center, y_center, radius : integer);
Var p, x, y : integer; Procedure plot_circle_points;
begin
set_pixel (x_center + x , y_center + y); 
set_pixel (x_center – x , y_center + y); 
set_pixel (x_center + x , y_center - y); 
set_pixel (x_center – x , y_center - y); 
set_pixel (x_center + y , y_center + x); 
set_pixel (x_center – y , y_center + x); 
set_pixel (x_center + y , y_center - x); 
set_pixel (x_center – y , y_center - x);
end; //{plot_circle_points}

begin //{bres_circle} 
	x :=0;
	y := radius;
	p := 3 - 2 * radius; 
	while x < y do begin
		plot_circle_points;
		if p < 0 then p := p + 4 * x + 6 
		else begin
			p := p + 4 * (x - y) + 10;
			y := y - 1 
		end;//{if p not < 0}
		x := x + 1 
	end; //{while x < y}
	if x = y then plot_circle_points 
end; //{bres_circle}

procedure boundary_fill ( x, y, fill_color, boundary, :integer); 
var present_color : integer;
begin
	present_color := inquire_color (x, y);
	if (present_color<>boundary and present_color<>fill_color)then begin
		set_pixel (x, y, fill_color);
		boundary_fill (x+1, y, fill_color, boundary); 
		boundary_fill (x-1, y, fill_color, boundary); 
		boundary_fill (x, y+1, fill_color, boundary); 
		boundary_fill (x, y-1, fill_color, boundary);
	end; // if present color 
end; //bondary fill

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Bezier curves
x(t) = a_x * t ^ 3 + b_x * t ^ 2 + c_x * t + d;
y(t) = a_y * t ^ 3 + b_y * t ^ 2 + c_y * t + d;
(0 <= t <= 1)

x(t) = T * C * x; 					% T and C are vectors
T = [t3 t2 t1 1];
Cx = [a b c d];

x'(0) = 3 * (x_2 - x_1);
x'(1) = 3 * (x_4 - x_3);

% Parameter
x(0) = x_1;
x(1) = x_4;

$ Polynom
x'(t) = 3 * a_x * t ^ 2 + 2 * b_x * t + c_x;
x(t) = a_x * t ^ 3 + b_x * t ^ 2 + c_x * t + d_x;					% (a)
x'(0) = 3 * (x_2 - x_1);											% (b)
x'(1) = 3 * (x_4 - x_3);											% (b)
x(0) = x_1;															% (c)
x(1) = x_4;															% (c)
x'(t) = 3 * a_x * t ^ 2 + 2 * b_x * t + c_x;						% (d)
x(0) = d_x;															% (e)
x(1) = A_x + b_x + c_x + d_x;										% (e)
x'(0) = c_x;														% (e)
x'(1) = 3 * a_x + 2 * b_x + c_x;									% (e)
c_x = 3 * x_2  - 3 * x_1;											% (f)
3 * a_x + 2 * b_x + c_x = 3 * x_4 - 3 * x_3;						% (f)
d_x = x_1;															% (f)
a_1 + b_x + c_x + d_x = x_4;										% (f)

a_x = -x_1 + 3 * x_2 - 3* x_3 + x_4;
b_x =  3 * x_1 - 6 * x_2 + 3 * x_3;
c_x = -3 * x_1 + 3 * x_2;
d_x = x_1;

% Matrix Representation: C = Mb * P_x
Mb = [-1 3 -3 1; 3 -6 3 0; -3 3 0 0; 1 0 0 0];
P_x = [x_1; x_2; x_3; x_4];
x(t) = T * Mb * P_x;
x(t) = [t3 t2 t1 1] * [-1 3 -3 1; 3 -6 3 0; -3 3 0 0; 1 0 0 0] * [x_1; x_2; x_3; x_4];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Transformation 2D

% Scaling
X' = Sx * X;
Y' = Sy * Y;
[X' Y'] = [X Y] *[Sx 0; 0 Sy];											% Matrix Representation

% Rotation
X' = r * cos(phi + theta);
Y' = r * sin(phi + theta);
X = r * cos(phi);
Y = r * sin(phi);
X' = r * cos(phi) * cos (theta) - r * sin(phi) * sin(theta);
Y' = r * sin(phi) * cos (theta) + r * cos(phi) * sin(theta);
X' = X * cos(theta) - Y * sin(theta);
Y' = X * sin(theta) + Y * cos(theta);
[X' Y'] = [X Y] * [cos(theta) sin(theta); -sin(theta) cos(theta)];		% Matrix Representation


% Translation (NOT Linear)
X' = X + T_x;
Y' = Y + T_y;

%Linear if:
((L(x+y) = L(x) + L(y)) && (L(cX) = cL(X))) for each scalar c

% Affine transformation
A(x) = L(x) + t;
P = [x y 1];
T_1 = [a(1,1) a(1,2); a(2,1) a(2,2); a(3,1) a(3,2)];
T_1 = [a(1,1) a(1,2) 0; a(2,1) a(2,2) 0; a(3,1) a(3,2) w];

[x' y' 1] = [x y 1 ] * [1 0 0; 0 1 0; T_x T_y 1];							% Translation
[x y 1] = [u v 1] * [cos(phi) sin(phi) 0; -sin(phi) cos(phi) 0; 0 0 1];		% Rotation
[x y 1] = [u v 1] * [S_u 0 0; 0 S_v 0; 0 0 1];								% Scaling

% Matrix
T = [1 0 0; 0 1 0; T_x T_y 1];												% Translation
S = [S_x 0 0; 0 S_y 0; 0 0 1];												% Scaling
R = [cos(phi) sin(phi) 0; -sin(phi) cos(phi) 0; 0 0 1];						% Rotation

% Mirror y
M_y = [-1 0 0; 0 1 0; 0 0 1];

% Mirror x
M_x = [1 0 0; 0 -1 0; 0 0 1];

% Reflection xy
M_xy = [-1 0 0; 0 -1 0; 0 0 1];

%Shearing
x' = x;
y' = y + b * x;

x' = x + a * y;
y' = y;

Sh_x = [1 0 0; a 1 0; 0 0 1];
Sh_y = [1 b 0; 0 1 0; 0 0 1];

% Invertible
T^-1 = [1 0 0; 0 1 0; -T_x -T_y 1];
S^-1 = [1/S_x 0 0; 0 1/S_y 0; 0 0 1; 0 0 1];
R^-1 = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];

% Scaling from point
M = T(-X_c, -Y_c) * S(S_x, S_y) * T(X_c, Y_c);
M = [1 0 0; 0 1 0; -T_x -T_y 1] * [S_x 0 0; 0 S_y 0; 0 0 1; 0 0 1] * [1 0 0; 0 1 0; T_x T_y 1];
M = [1 0 0; 0 1 0; -x_c -y_c 1] * [S_x 0 0; 0 S_y 0; 0 0 1; 0 0 1] * [1 0 0; 0 1 0; x_c y_c 1];
M = [S_x 0 0; 0 S_y 0; x_c(1-S_x) y_c(1-S_y) 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Window
[Wxl, Wxr, Wyt, Wyb];

% View Port
[Vxl, Vxr, Vyt, Vyb];

% Scissoring - Image Oriented
% Clipping - Object Oriented

Mv = Tw * Swv * Tv;
Tw = [1 0 0; 0 1 0; -Wxl -Wyb 1];
Tv = [1 0 0; 0 1 0; Vxl Vyb 1];
Swv = [S_x 0 0; 0 S_y 0; 0 0 1];
S_x = (Vxr - Vxl) / (Wxr - Wxl);
S_y = (Vyt - Vyb) / (Wyt - Wyb);

% Equation Mapping
Mwv = [S_x 0 0; 0 S_y 0; Vxl-S_x*Wxl Vyb-S_y*Wyb 1];
X'=S_x*(X-Wxl)+Vxl;
Y'=S_y*(Y-Wyb)+Vyb;

% Cohen-Sutherland Algorithm
% Cyrus-Beck Clipping
