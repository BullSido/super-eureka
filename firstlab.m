clear;
clc;
%7th variant

%ready to hand in
%
%
%config_values
%given values
r = 0.66;
K=202;
delta = r/K;
obj_name = "Paramecium";
start_point = 0;
end_point = 10;
number_of_steps = 31;
starting_population = 8;

%first
fneeds_plot_approximate = false;
fneeds_plot_precise = false;
fneeds_xl_table = false;

%second
second_r = delta;
second_starting_population = 8;

sneeds_plot_approximate = false;
sneeds_plot_precise = false;
sneeds_errors_calculated = false;
sneeds_plot = false;

%third
third_starting_population = 8;
talpha = 1;
tbeta = 5;
ttau = 20;

tneeds_plot_approximate = false;
tneeds_plot_precise = false;
tneeds_errors_calculated = false;
tneeds_plot = false;

%fourth
fth_starting_population = 8;
fthalpha = 1;
fthbeta = 15;
fthtau = 20;
fthd = 0.1;

fth_needs_plot = false;
fth_needs_errors_calculated = false;

%fivth
fv_starting_population = 8;
fvtau = 20;
fvalpha = 2;
fvbeta = 5;
fvdelta = 1 / K ;
fvd = 1/K;

fv_needs_plot = true;
fv_needs_errors_calculated = true;
%end_of_config_values
%
%
%

step_length=(end_point-start_point)/number_of_steps;
Time = start_point :step_length : end_point; 


%first task
%dy/dt = r * y - delta * y^2
%approximation
if fneeds_plot_approximate || fneeds_xl_table
    population_vector=zeros(1,number_of_steps+1);
    population_vector(1) = starting_population;
    for i = 1:number_of_steps
        population_vector(i+1) = population_vector(i) + step_length * population_vector(i) * (r - delta * population_vector(i) );
    end
end
%precise
if fneeds_plot_precise || fneeds_xl_table
    precise_population = zeros(1,number_of_steps+1);
    Coef = (r-delta*starting_population)*exp(r*start_point)/(starting_population*r);
    for i=1:number_of_steps+1
    precise_population(i)=r/(delta+Coef * r * exp(-r * Time(i)));
    end
end

if fneeds_plot_precise || fneeds_plot_approximate
    figure();
    hold on
        title("first");
        xlabel("time");
        ylabel(obj_name);
        if fneeds_plot_approximate
            plot(Time, population_vector, "m*-");
        end
        if fneeds_plot_precise
            plot(Time,precise_population,"bo-");
        end
    hold off
end

if fneeds_xl_table
    absolute_error = abs(precise_population - population_vector);
    relative_error = absolute_error./precise_population;
    points = 1:(number_of_steps+1);

    Names={'Time t','Exact solution','Numerical solution','Abs','Otn'}; 
    xlswrite('LR02.xls',points,'a2:a32');
    xlswrite('LR02.xls',Names,'b1:f1');
    xlswrite('LR02',Time','b2:b32');
    xlswrite('LR02.xls',population_vector','c2:c32');
    xlswrite('LR02.xls',precise_population','d2:d32');
    xlswrite('LR02.xls',absolute_error','e2:e32'); 
    xlswrite('LR02.xls',relative_error','f2:f32');
end


%second task
%dy/dt = r * y^2
starting_population = second_starting_population;

%approximate
if sneeds_plot_approximate || sneeds_errors_calculated || sneeds_plot
    sec_approx_population = zeros(1,number_of_steps + 1);
    sec_approx_population(1) = starting_population;
    for i = 1:number_of_steps
        sec_approx_population(i+1) = sec_approx_population(i) + step_length * sec_approx_population(i)  * sec_approx_population(i) * second_r ;
    end
end

%precise
if sneeds_plot_precise || sneeds_errors_calculated || sneeds_plot
    sec_prec_population = zeros(1,number_of_steps+1);
    Coef = 1/starting_population + second_r * start_point;
    for i=1:number_of_steps+1
    sec_prec_population(i)= - 1 / (second_r * Time(i) - Coef);
    end
end

if sneeds_errors_calculated
    sec_absolute_error = abs(sec_prec_population - sec_approx_population);
    sec_relative_error = sec_absolute_error./sec_prec_population;
    figure();
    hold on
        title("second absolute error");
        xlabel("time");
        ylabel(obj_name);
        plot(Time, sec_absolute_error, "m*-");
    hold off
    figure();
    hold on
        title("second relative error");
        xlabel("time");
        ylabel(obj_name);
        plot(Time, sec_relative_error, "m*-");
    hold off
end

if sneeds_plot
    figure();
    hold on
        title("second");
        xlabel("time");
        ylabel(obj_name);
        plot(Time, sec_approx_population, "m*-");
        plot(Time,sec_prec_population,"bo-");
    hold off
end

%third task
%dy/dt = alpha * beta * y^2 / (beta + tau * y)
starting_population = third_starting_population;

if tneeds_plot_precise || tneeds_errors_calculated || tneeds_plot
    thr_prec_population = zeros(1,number_of_steps+1);
    syms y(t);
    condition = y(start_point) == starting_population;
    ode = diff(y,t) == talpha * tbeta * y^2 / (tbeta + ttau * y);
    thr_prec_sym(t) = dsolve(ode, condition);
    for i = 1:number_of_steps+1
        thr_prec_population(i) = thr_prec_sym(Time(i));
    end
end

if tneeds_plot_approximate || tneeds_errors_calculated || tneeds_plot
    thr_approx_population = zeros(1,number_of_steps + 1);
    thr_approx_population(1) = starting_population;
    for i = 1:number_of_steps
        thr_approx_population(i+1) = thr_approx_population(i) + step_length * (talpha * tbeta * thr_approx_population(i) * thr_approx_population(i) / (tbeta + ttau * thr_approx_population(i))) ;
    end
end

if tneeds_plot
    figure();
    hold on
        title("third");
        xlabel("time");
        ylabel(obj_name);
        plot(Time, thr_approx_population, "m*-");
        plot(Time, thr_prec_population, "bo-");
   hold off
end

if tneeds_errors_calculated
    thr_absolute_error = abs(thr_prec_population - thr_approx_population);
    thr_relative_error = thr_absolute_error./thr_prec_population;
    figure();
    hold on
        title("third absolute error");
        xlabel("time");
        ylabel(obj_name);
        plot(Time, thr_absolute_error, "m*-");
    hold off
    figure();
    hold on
        title("third relative error");
        xlabel("time");
        ylabel(obj_name);
        plot(Time, thr_relative_error, "m*-");
    hold off
end


%fourth task
%dy/dt = alpha * beta * y^2 / (beta + tau * y) - dy
starting_population = fth_starting_population;

%approx
if fth_needs_errors_calculated || fth_needs_plot
    fth_approx_population = zeros(1,number_of_steps + 1);
    fth_approx_population(1) = starting_population;
    for i = 1:number_of_steps
        fth_approx_population(i+1) =  fth_approx_population(i) + step_length * (fthalpha * fthbeta * fth_approx_population(i) * fth_approx_population(i) / (fthbeta + fthtau * fth_approx_population(i)) - fthd * fth_approx_population(i)) ;
    end
end

%precise
if fth_needs_errors_calculated || fth_needs_plot
    fth_prec_population = zeros(1,number_of_steps+1);
    syms y(t);
    condition = y(start_point) == starting_population;
    ode = diff(y,t)  == (((fthalpha * fthbeta - fthd * (fthbeta + fthtau)) * y ^ 2)/ (fthbeta + fthtau * y) );
    fth_prec_sym(t) = dsolve(ode, condition);
    for i = 1:number_of_steps+1
        fth_prec_population(i) = fth_prec_sym(Time(i));
    end
end

if fth_needs_plot
    figure();
    hold on
        title("fourth");
        xlabel("time");
        ylabel(obj_name);
        plot(Time, fth_approx_population, "m*-");
        plot(Time, fth_prec_population, "bo-");
   hold off
end

if fth_needs_errors_calculated
    fth_absolute_error = abs(fth_prec_population - fth_approx_population);
    fth_relative_error = fth_absolute_error./fth_prec_population;
    figure();
    hold on
        title("fourth absolute error");
        xlabel("time");
        ylabel(obj_name);
        plot(Time, fth_absolute_error, "m*-");
    hold off
    figure();
    hold on
        title("fourth relative error");
        xlabel("time");
        ylabel(obj_name);
        plot(Time, fth_relative_error, "m*-");
    hold off
end


%fivth
starting_population = fv_starting_population;

%approx
if fv_needs_errors_calculated || fv_needs_plot
    fv_approx_population = zeros(1, number_of_steps+1);
    fv_approx_population(1) = starting_population;
    for i=1:number_of_steps
        fv_approx_population(i+1) = fv_approx_population(i) + step_length * fvalpha * fvbeta * fv_approx_population(i) * fv_approx_population(i) / (fvbeta + fvtau * fv_approx_population(i)) - fvd * fv_approx_population(i) - fvdelta * fv_approx_population(i) * fv_approx_population(i);
    end
end

%precise
if fv_needs_errors_calculated || fv_needs_plot
    [Time, fv_prec_population] = ode45(@(t,x)fvalpha*(fvbeta * x^2)/(fvbeta+fvtau*x)-fvd-fvdelta*x^2,start_point:step_length:end_point, starting_population);
    fv_prec_population = transpose(fv_prec_population);
end

if fv_needs_plot
    figure();
    hold on
        title("fivth");
        xlabel("time");
        ylabel(obj_name);
        plot(Time, fv_approx_population, "m*-");
        plot(Time, fv_prec_population, "bo-");
   hold off
end

if fv_needs_errors_calculated
    fv_absolute_error = abs(fv_prec_population - fv_approx_population);
    fv_relative_error = fv_absolute_error./fv_prec_population;
    figure();
    hold on
        title("fivth absolute error");
        xlabel("time");
        ylabel(obj_name);
        plot(Time, fv_absolute_error, "m*-");
    hold off
    figure();
    hold on
        title("fivth relative error");
        xlabel("time");
        ylabel(obj_name);
        plot(Time, fv_relative_error, "m*-");
    hold off
end