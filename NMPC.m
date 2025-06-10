function u=NMPC(CSTR,x0,N,xr,w1,w2,w3,flag)
    %the idea is to model the problem as a nonlinear program which has dimensions and constraints depending on N: 
    %our variables will be the N+1 x, and the N u (2N+2+N)

    %parameters
    Ts=0.2;
    teta=20;
    K=300;
    M=5;
    xf=0.3947;
    xc=0.3816;
    alpha=0.117;

    %linear equality constraints for initial conditions (Aeq*x=x0)
    Aeq=eye(2,3*N+2); %extracts (x(1:2))
    beq=x0;

    %upper bounds and lower bounds
    ub=zeros(3*N+2,1);
    lb=zeros(3*N+2,1);
    state_u=[1;0.66];
    input_u=2;
    state_l=[0;0];
    input_l=0.1;
    %the state and input constraints are true at every time step so we need
    %to stack them in the right way to have lb<=x<=ub
    for i=0:N
        ub(2*i+1:2*i+2)=state_u; %state constraints
        lb(2*i+1:2*i+2)=state_l;  
    end
    for i=1:N
        ub(2*N+2+i)=input_u; %input constraints
        lb(2*N+2+i)=input_l;
    end

    %nonlinear equalities for system evolution x(i+1)=cstr(x(i),u(i): R^(3N+2)--> R^(2N)
    function f = make_fun(A,Ts, teta, K, xf, M, alpha, xc)
        c = @(x) x(3:4) - CSTR(x(1:2), Ts, teta, K, xf, M, alpha, xc, x(5));
        f = @(x) c(A * x);
    end
    f_vec = cell(N,1); %container for our N constraints for the evolution of the states

    for i = 1:N
      A_i = zeros(5, 3*N + 2);
      A_i(1:4, 2*i-1:2*i+2) = eye(4); %selecting the xi and x(i+1) state
      A_i(5, 2*N + 2 + i) = 1; %selecting the u(i) input

     % now we ensure that each f_vec has its own A_i and don't point to the same matrix
     %since it changes during the loop matlab uses the current value of the
     %parameter, which we don't want
      f_vec{i} = make_fun(A_i,Ts, teta, K, xf, M, alpha, xc); 
    end                                                       
                                                              
    %vectorizing the function
    Feq = @(x) expand_feq(x, f_vec);
    %we need to call another non-anonymous function because matlab behaved
    %weirdly: it copied two times the output so I just cut it in half
    function y = expand_feq(x, f_vec)
        y_cells = cellfun(@(f) f(x), f_vec, 'UniformOutput', false); %now each cell is a function handle acting on x
        y = vertcat(y_cells{:});
        y=y(1:2*N)'; 
    end

    %objective
    tracking=@(x) w1*sum((x-ones(N,1)*xr).^2); %it will receive the vectors of x_2s
    terminal=@(x) w3*(x-xr)^2; %it will receive only the last x_2
    if flag==1
        input=@(x) w2*sum((x(2:end)-x(1:end-1)).^2); %it receives the entire input sequence
    else  %giving the possibility to change objective 
        input=@(x) w2*sum(x.^2);
    end
    %selecting the second variable of each x until the N'th
    B=zeros(N,3*N+2);
    for i=1:N
        B(i,2*i)=1;
    end

    objective= @(x) tracking(B*x)+terminal(x(2*N+2))+input(x(2*N+3:end));

    deq=zeros(2*N,1);
    %size(deq) == size(Feq(zeros(3*N+2)))

    %now we can define the optimizer
    opts = optiset('solver','ipopt','display','iter');
    Opt=opti('ndec',3*N+2,'fun',objective,'bounds',lb,ub, 'Aeq',Aeq,'beq',beq,'nlmix',Feq,deq,zeros(2*N,1),'options',opts);
    in=0.5*ones(3*N+2,1);  %choosing as start a random point (not feasible)
    [x,~,~,~] = solve(Opt,in);
    u=x(2*N+3); %extract the input u0 as the 2N+2+1 variable
end