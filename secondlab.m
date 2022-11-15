clc;
clear;

%starting conditions
starting_phi = pi/2;
friction_coeff = 1;
mass = 10;

g = 9.81;
L = 6; %length of rope


%decide which
first = true;

%solved equation
if first
    w = sqrt(g/L);
else
    w = sqrt(g/L - friction_coeff^2 / (4 * mass^2));
end

T = 2 * pi * sqrt(L/g);

syms phi(t) y(t);
if first
    phi(t) = starting_phi * cos(w*t);
    x(t) = sin(phi(t));
    y(t) = -cos(phi(t));
    
    fplot(x(t*T)/pi, [0 5]); grid on;
    title("Гармоническое движение маятника");
    xlabel('t/T');
    ylabel('\theta/\pi');
    
    figure();
    hold on;
    fanimator(@fplot,x,y,'ko','MarkerFaceColor','k','AnimationRange',[0 5*T]);%отображение шара
    axis equal;
    fanimator(@(t) plot([0 x(t)],[0 y(t)],'k-'),'AnimationRange',[0 5*T]);%отображение нити
    fanimator(@(t) text(-0.3,0.3,"Время прошло: "+num2str(t,2)+" с"),'AnimationRange',[0 5*T]);%отображение таймера
    xlim([-1.5 1.5]);
    ylim([-1.5 1.5]);
    playAnimation;
    hold off;
else
    phi(t) =starting_phi *exp(-friction_coeff * t / (2 * mass)) * cos(w*t);
    x(t) = sin(phi(t));
    y(t) = -cos(phi(t));
    
    fplot(x(t*T)/pi, [0 5]); grid on;
    title("Затухающее движение маятника");
    xlabel('t/T');
    ylabel('\theta/\pi');
    
    figure();
    hold on;
    fanimator(@fplot,x,y,'ko','MarkerFaceColor','k','AnimationRange',[0 5*T]);%отображение шара
    axis equal;
    fanimator(@(t) plot([0 x(t)],[0 y(t)],'k-'),'AnimationRange',[0 5*T]);%отображение нити
    fanimator(@(t) text(-0.3,0.3,"Время прошло: "+num2str(t,2)+" с"),'AnimationRange',[0 5*T]);%отображение таймера
    xlim([-1.5 1.5]);
    ylim([-1.5 1.5]);
    playAnimation;
    hold off;
end

