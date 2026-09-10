clc;
clear all;
close all;
dr=200;t=12;r=8;brlos=10^((-32.4-30*log10(dr)-20*log10(3.5))/10);
fc=3.5*10^9;B=60*10^6;P=10^-3*10^(24/10);sigma=10^-3*10^(-17.4)*B;
n=30;epsilon=10^-6;
a=log2(exp(1))*qfuncinv(epsilon)/sqrt(n);
Ns=min(t,r);
rho=1;etarho=2;
mont=5000;
D=5;
targetLen = 25;

repro_cell = cell(1, mont);
repe_cell  = cell(1, mont);

for w = 1:mont
    H = sqrt(brlos) * (randn(t, r) + 1i * randn(t, r)) / sqrt(2);
    lambda = svd(H * H');
    gamma = lambda.' ./ sigma;
    gamma = gamma(1:Ns);
    Q = rand(Ns, Ns);
    Q = Q * P / (sum(Q(:)));
    repro_cell{w} = runSTCCproposed(gamma, P, D, a, Ns, rho, etarho, Q);
    repe_cell{w}  = runSTCCpe(gamma, P, D, a, Ns, Q);
end

repro25 = zeros(targetLen, mont);
repe25  = zeros(targetLen, mont);

for w = 1:mont
    x = repro_cell{w}(:);
    y = repe_cell{w}(:);

    % 不足 25 次时，用最后一个值补齐
    if numel(x) < targetLen
        x(end+1:targetLen, 1) = x(end);
    else
        x = x(1:targetLen);
    end

    if numel(y) < targetLen
        y(end+1:targetLen, 1) = y(end);
    else
        y = y(1:targetLen);
    end

    repro25(:, w) = x;
    repe25(:, w)  = y;
end

% 对所有列求均值
repro_mean = mean(repro25, 2);
repe_mean  = mean(repe25, 2);

%% plot
X=1:25;
figure;
set(0,'DefaultAxesFontName','Times New Roman')
set(0,'DefaultTextFontName','Times New Roman')
plot(X,repro_mean,'-r*',X,repe_mean,'-bd','LineWidth',1.5);
xlim([1,25]);
xlabel('Index of iteration $\tau$', 'Interpreter', 'latex','FontSize', 12);
ylabel('Rate / bit/channel use','FontSize', 12)
legend('STCC-PACA','STCC-BMCA','FontSize', 12);
legend('location','northwest');
grid on;
