clc;

year_start = 1945;
year_end = 1970;
y0 = 10; % start population

r = 1.1;
K = 4500;

t0 = 0;
t1 = year_end - year_start;

n = 40;
h = (t1 - t0) / n;
Y = zeros(1, n + 1);
Y(1) = y0;

delta = r / K;

YTOCHN = zeros(1, n + 1);
T = t0:h:t1;
C = (r - delta * y0) * exp(r * t0) / (y0 * r);

for i = 1: n + 1
    YTOCHN(i) = r / (delta + C * r * exp(-r * T(i)));
end

figure
plot(T, YTOCHN, 'bo-');
xlabel('Years');
ylabel('Популяция белощокой Казарки');
axis([0 t1 0 K+50]);
set(gca, 'XTickLabel', {year_start:2:year_end});

for i = 1 : n
    Y(i + 1) = Y(i) + h * Y(i) * (r - delta * Y(i));
end
abs1 = abs(YTOCHN - Y);
otn1=abs1 ./ YTOCHN;
points = 1 : (n + 1);

plot(T, Y, 'm*-');
hold on;
grid on; 
plot(T, YTOCHN, 'bo-');
xlabel('Years');
ylabel('Популяция белощекой казарки');
axis([0 t1 0 K+50]); 
set(gca,'XTickLabel',{year_start:2:year_end});

legend('Приближенное значение', 'Точное значение');
%2


YTOCHN = zeros(1, n + 1);
C1 = r * t0 + 1 / y0;
for i = 1: n + 1
    YTOCHN(i) = 1 / (C1 - r * T(i));
end

%YTOCHN = zeros(1, n + 1);
Y = zeros(1, n + 1);
Y(1) = y0;

%[T1, YTOCHN] = ode45(@(t, x)r * x^2, t0:h:t1, y0);

for i = 1 : n
    Y(i + 1) = Y(i) + h * r * Y(i)^2;
end



abs2 = abs(YTOCHN - Y)
otn2 = abs2 ./ YTOCHN


figure('NumberTitle', 'off', 'Name','Задание 2')
plot(T + year_start, YTOCHN, 'bo-', T + year_start, Y, 'g*-');
hold on;
grid on; 
xlabel('Года');
ylabel('Популяция белощекой казарки');
axis([year_start year_end 0 K * 1.1]);
legend('Точное значение', 'Приближенное значение');

% 3

y0 = 1;
tau = 3;
a = 1;
b = 3;

YTOCHN = zeros(1, n + 1);
Y = zeros(1, n + 1);
Y(1) = y0;

[T1, YTOCHN] = ode45(@(t, x)a * (b * x^2)/(b + tau * x), t0:h:t1, y0);

for i = 1:n
    Y(i + 1) = Y(i) + h * a * b * Y(i)^2 / (b + tau * Y(i));
end

abs3 = abs(YTOCHN - Y)
otn3 = abs3 ./ YTOCHN

figure('NumberTitle', 'off', 'Name','Задание 3')
plot(T + year_start, YTOCHN, 'bo-', T + year_start, Y, 'g*-');
hold on;
grid on; 
xlabel('Года');
ylabel('Популяция белощекой казарки');
axis([year_start year_end 0 K * 1.1]);
legend('Точное значение', 'Приближенное значение');

%4

y0 = 1;
tau = 1;
a = 0.6;
b = 1;
d = 0.1;

YTOCHN = zeros(1, n+1);
Y = zeros(1, n+1);
Y(1) = y0;

[T1, YTOCHN] = ode45(@(t, x)a * (b * x^2)/(b + tau * x) - d, t0:h:t1, y0);

for i = 1 : n
    Y(i+1) = Y(i) + h * a * b * Y(i)^2 / (b + tau * Y(i)) - d * h * Y(i);
end

abs4 = abs(YTOCHN - Y);
otn4 = abs4 ./ YTOCHN;

figure('NumberTitle', 'off', 'Name','Задание 4')
plot(T + year_start, YTOCHN, 'bo-', T + year_start, Y, 'g*-');
hold on;
grid on; 
xlabel('Года');
ylabel('Популяция белощекой казарки');
axis([year_start year_end 0 K * 1.1]);
legend('Точное значение', 'Приближенное значение');

% 5

y0 = 1;
tau = 1;
a = 2;
b = 1;
delta = r / K;
d = 1 / K; 

YTOCHN = zeros(1, n + 1);
Y = zeros(1, n + 1);
Y(1) = y0;

[T1, YTOCHN] = ode45(@(t, x)a * (b * x^2)/(b + tau * x) - d - delta * x^2, t0:h:t1, y0);

for i = 1 : n
    Y(i + 1) = Y(i) + h * a * b * Y(i)^2 / (b + tau * Y(i)) - d * Y(i) - delta * Y(i)^2;
end

abs5 = abs(YTOCHN - Y);
otn5 = abs1 ./ YTOCHN;

figure('NumberTitle', 'off', 'Name','Задание 5')
plot(T + year_start, YTOCHN, 'bo-', T + year_start, Y, 'g*-');
hold on;
grid on; 
xlabel('Года');
ylabel('Популяция белощекой казарки');
axis([year_start year_end 0 2 * K * 1.1]);
legend('Точное значение', 'Приближенное значение');

