clear,clc,format long, format compact % Nicole Hemming-Schroeder
% Project 1: Joos, et al, 2013, Carbon Response Function Analysis 
syms t, a0=0.2173; a=[0.2240, 0.2824, 0.2763];
time=[394.4, 36.54, 4.304];
f=a0;
for i=1:3
    f=f+a(i)*exp(-t/time(i));
end
% for 0<t<1000
percent=[0.25,0.6];

bisectionapprox=[500,1000;19,21];
newtonapprox=[1000;20];
secantapprox=[950,1000;19,21];
for k=1:2
    fstar=f-percent(k); dfstar=diff(fstar,t,1);
    
    % Bisection Method
    n=100; error=0.0001; a=zeros(1,n);b=zeros(1,n);p=zeros(1,n);
    a(1)=bisectionapprox(k,1); b(1)=bisectionapprox(k,2);
    for i=1:n
        p(i)=a(i)+((b(i)-a(i))/2);
        if subs(fstar,t,a(i))*subs(fstar,t,p(i))<0 % then it crossed
            a(i+1)=a(i); b(i+1)=p(i);
        elseif subs(fstar,t,a(i))*subs(fstar,t,p(i))>0 % then it didn't cross
            a(i+1)=p(i); b(i+1)=b(i);
        elseif subs(fstar,t,a(i))*subs(fstar,t,p(i))==0
            disp('Bisection root found'), bisectionroot=p(i), disp(' '), break
        end
        if abs(b(i+1)-a(i+1))<error
            disp('Bisection within error'),bisectionroot=p(i),disp(' '),
            i, break 
        end
    end
    % Newton's Method
    tn(1)=newtonapprox(k,1);
    for i=1:n
        tn(i+1)=tn(i)-subs(fstar,t,tn(i))/subs(dfstar,t,tn(i));
        if subs(fstar,t,tn(i+1))==0
            disp('Newton root found'), newtonroot=tn(i+1), disp(' '), break
        elseif abs(tn(i+1)-tn(i))<error
            disp('Newton within error margin'), newtonroot=tn(i+1), disp(' ')
            i, break 
        end
    end
    % Secant Method
    x(1)=secantapprox(k,1); x(2)=secantapprox(k,2);
    for i=1:n
        x(i+2)=x(i+1)-subs(fstar,t,x(i+1))*(x(i+1)-x(i))/...
            (subs(fstar,t,x(i+1))-subs(fstar,t,x(i)));
        if subs(fstar,t,x(i+2))==0
            disp('Secant root found'), secantroot=x(i+2), break
        elseif abs(x(i+2)-x(i+1))<error
            disp('Secant within error margin'), secantroot=x(i+2)
            i, break 
        end
    end
end
figure
hold on
set(0,'DefaultTextFontsize',18, ...
    'DefaultTextFontname','Times New Roman', ...
    'DefaultTextFontWeight','bold', ...
    'DefaultAxesFontsize',18, ...
    'DefaultAxesFontname','Times New Roman', ...
    'DefaultLineLineWidth', 4)
ax = gca; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
axis([0, 30, 0, 1])
xlabel('Time (years)')
ylabel('Proportion of Initial CO_2 Pulse')
fplot(f,[0,30],'r')
fplot(0.6,[0,30],'k--')
legend('CO_2 Pulse','60% of the Initial Carbon')
hold off

figure, hold on
set(0,'DefaultTextFontsize',18, ...
    'DefaultTextFontname','Times New Roman', ...
    'DefaultTextFontWeight','bold', ...
    'DefaultAxesFontsize',18, ...
    'DefaultAxesFontname','Times New Roman', ...
    'DefaultLineLineWidth', 4)
ax = gca; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
axis([19.2, 19.5, 0.6-0.0006,0.6+0.0006]) 
scatter(p(3),0.6,'b','filled'), scatter(tn(3),0.6,'k','filled'), 
scatter(x(3),0.6,'m','filled'), fplot(f,'r'),fplot(0.6,'k--')
xlabel('Time (years)')
ylabel('Proportion of Initial CO_2 Pulse')
legend('Bisection Method','Newton''s Method','Secant Method','f(t)',...
'g(t)')
hold off
