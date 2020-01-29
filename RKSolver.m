function h_numerical = RKSolver(x_lower, x_upper, x_steps, dograph)
%   Rung Kutta Solver for Lane Emden Equations


x0 = x_lower;                   %Lower Bound of interval
xn = x_upper;                   %End of Interval on value line x
n = x_steps;                    %Number of steps
dx = (xn-x0)/(n-1);                 %Increments dividing 0 to pi equally with n steps steps

h(1)    = 1;                    %Initial Condition
f(1)    = 0;                    
x(1)    = x0;

h_dot   = @(f) f;
f_dot   = @(x,h,f) -(2/x)*f - h;
f_dot_singularity = @(h) -(1/3); 

x=0

%If statement to mitigate singularity problem
i = 1;
if x == 0 
        
    h_1 = h_dot(f(i));
    f_1 = f_dot_singularity(h(i));
        
    h_2 = h_dot(f(i)+dx*f_1/2);
    f_2 = f_dot_singularity(h(i)+h_1*dx/2);
        
    h_3 = h_dot(f(i)+dx*f_2/2);
    f_3 = f_dot_singularity(h(i)+h_2*dx/2);
        
    h_4 = h_dot(f(i)+f_3*dx);
    f_4 = f_dot_singularity(h(i)+h_3*dx);
        
    h(i+1) = h(i) + (h_1 + 2*h_2 + 2*h_3 + h_4)*dx/6;     
    f(i+1) = f(i) + (f_1 + 2*f_2 + 2*f_3 + f_4)*dx/6;
        
    i = i+1;
    x= x+dx;
end
                        
for i = 2:n-1
     
    h_1 = h_dot(f(i));
    f_1 = f_dot(x,h(i), f(i));
    
    h_2 = h_dot(f(i)+dx*f_1/2);
    f_2 = f_dot(x+dx/2, h(i)+h_1*dx/2, f(i)+f_1*dx/2);
    
    h_3 = h_dot(f(i)+dx*f_2/2);
    f_3 = f_dot(x+dx/2, h(i)+h_2*dx/2, f(i)+f_2*dx/2);
    
    h_4 = h_dot(f(i)+f_3*dx);
    f_4 = f_dot(x+dx, h(i)+h_3*dx, f(i)+f_3*dx);
     
    h(i+1) = h(i) + (h_1 + 2*h_2 + 2*h_3 + h_4)*dx/6;     
    f(i+1) = f(i) + (f_1 + 2*f_2 + 2*f_3 + f_4)*dx/6;    
    x = x+dx;

end
%store the results in h_numerical (for consistency to the question sheet)
h_numerical = h;
format long


%if the function parameters require a graph, we produce a graph. 

if dograph == true
    xt = linspace(x0, xn,length(h));
    
    plot(xt,h, ':b')
    grid on
    
    xlabel('x')
    xticks([0, pi/4, pi/2, 3*pi/4 pi])
    xticklabels({'0','\pi/4','\pi/2','3\pi/4','pi'})
    title('Numerical Computation of Lane Emden Equation')

    %we label the y axis with our equation name
    ylabel('Lane-Emden, Equation 1')
end


end

