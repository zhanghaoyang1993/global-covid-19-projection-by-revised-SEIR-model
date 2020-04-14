%%%%修订后的SEIR模型%%%%
%Powered by Denny·Z
%A Junior Metaheuristic Learner
%此模型为预测湖北省疫情数据
%比较是否考虑潜伏者的传染情况
%初始化
clear all
close all
clc
%参数设置
%模型初值设定
%S1=59170000;%湖北省人口总数
S1=59170000;%湖北省人口总数
N0=S1;
S2=S1;
S3=S1;
E1=4007;% 1月23至1月29日新增确诊病例数
E2=E1;
E3=E1;
I1=786;%感染者
I2=786;
I3=786;
Sq1=2776;%尚在接受医学观察的人数
Sq2=Sq1;
Sq3=Sq1;
Eq1=400;%估计值，为正在被隔离的潜伏者
Eq2=Eq1;
Eq3=Eq1;
H1=1186;%正在住院的患者，为感染者和被隔离的潜伏者之和
H2=H1;
H3=H1;
R1=31;%官方公布的出院人数
R2=R1;
R3=R1;
%模型参数设定
c1=2;%接触率
c2=2.2;%接触率
c3=2.5;%接触率
deltaI=0.13;%感染者的隔离速度
deltaq=0.13;%隔离潜伏者向隔离感染者的转化速率
gammaI=0.007;%感染者的恢复率
gammaH=0.014;%隔离感染者的恢复速率
beta1=2.05*10^(-9);%传染概率
beta2=2.05*10^(-9);%传染概率
beta3=2.05*10^(-9);%传染概率
q1=1*10^(-6);%隔离比例
q2=1*10^(-6);%隔离比例
q3=1*10^(-6)/2;%隔离比例
alpha1=2.7*10^(-4);%病死率
alpha2=2.7*10^(-4);%病死率
alpha3=2.7*10^(-4);%病死率
rho1=1;%有效接触系数，参考取1
rho2=1;%有效接触系数，参考取1
rho3=1;%有效接触系数，参考取1
theta1=0.2;%潜伏者相对于感染者的传染能力比值
theta2=0.2;%潜伏者相对于感染者的传染能力比值
theta3=0.2;%潜伏者相对于感染者的传染能力比值
lambda=1/14;%隔离解除速度，为14天的倒数
sigma=1/7;%潜伏者向感染者的转化速度，平均潜伏期为7天，为7天的倒数
%差分迭代方程
newinfected1(1)=0;
newinfected2(1)=0;
newinfected3(1)=0;
E_sum1=E1;%总感染人数
E_sum2=E1;%总感染人数
E_sum3=E1;%总感染人数
T=1:300;
for i =1:length(T)-1
    %%
    %乐观情景
    dS1=(rho1*c1*beta1+rho1*c1*q1*(1-beta1))*S1(i)*(I1(i)+theta1*E1(i))-lambda*Sq1(i);
    newinfected1(i+1)=rho1*c1*beta1*S1(i)*(I1(i)+theta1*E1(i));
    S1(i+1)=S1(i)-dS1;%易感人数迭代
    E1(i+1)=E1(i)+rho1*c1*beta1*(1-q1)*S1(i)*(I1(i)+theta1*E1(i))-sigma*E1(i);%潜伏者人数迭代
    I1(i+1)=I1(i)+sigma*E1(i)-(deltaI+alpha1+gammaI)*I1(i);%感染者人数迭代
    Sq1(i+1)=Sq1(i)+rho1*c1*q1*(1-beta1)*S1(i)*(I1(i)+theta1*E1(i))-lambda*Sq1(i);%隔离易感染着人数迭代
    Eq1(i+1)=Eq1(i)+rho1*c1*beta1*q1*S1(i)*(I1(i)+theta1*E1(i))-deltaq*Eq1(i);%隔离潜伏者人数迭代
    H1(i+1)=H1(i)+deltaI*I1(i)+deltaq+Eq1(i)-(alpha1+gammaH)*H1(i);%住院患者人数迭代
    R1(i+1)=R1(i)+gammaI*I1(i)+gammaH*H1(i);%康复人数迭代 
    E_sum1=E_sum1+newinfected1(i+1);
    if E_sum1/N0 > 1/2000
        rho1 = rho1/1.01;  %rho1/1.1，加强隔离措施
        q1=1*10^(-6)*1.01; %q1*1.3，扩大隔离比例
        beta1=2.05*10^(-9)/1.01;%传染概率降低
    end
    
    %%
    %中情景
    dS2=(rho2*c2*beta2+rho2*c2*q2*(1-beta2))*S2(i)*(I2(i)+theta2*E2(i))-lambda*Sq2(i);
    newinfected2(i+1)=rho2*c2*beta2*S2(i)*(I2(i)+theta2*E2(i));
    S2(i+1)=S2(i)-dS2;%易感人数迭代
    E2(i+1)=E2(i)+rho2*c2*beta2*(1-q2)*S2(i)*(I2(i)+theta2*E2(i))-sigma*E2(i);%潜伏者人数迭代
    I2(i+1)=I2(i)+sigma*E2(i)-(deltaI+alpha2+gammaI)*I2(i);%感染者人数迭代
    Sq2(i+1)=Sq2(i)+rho2*c2*q2*(1-beta2)*S2(i)*(I2(i)+theta2*E2(i))-lambda*Sq2(i);%隔离易感染着人数迭代
    Eq2(i+1)=Eq2(i)+rho2*c2*beta2*q2*S2(i)*(I2(i)+theta2*E2(i))-deltaq*Eq2(i);%隔离潜伏者人数迭代
    H2(i+1)=H2(i)+deltaI*I2(i)+deltaq+Eq2(i)-(alpha2+gammaH)*H2(i);%住院患者人数迭代
    R2(i+1)=R2(i)+gammaI*I2(i)+gammaH*H2(i);%康复人数迭代 
    E_sum2=E_sum2+newinfected2(i+1);
    if E_sum2/N0 > 1/100
        rho2 = rho2/1.002;
        q2=1*10^(-6)*1.01; %扩大隔离比例
        beta2=2.05*10^(-9)/1.005;%传染概率降低
    end
    
   %%
    %悲观情景
    dS3=(rho3*c3*beta3+rho3*c3*q3*(1-beta3))*S3(i)*(I3(i)+theta3*E3(i))-lambda*Sq3(i);
    newinfected3(i+1)=rho3*c3*beta3*S3(i)*(I3(i)+theta3*E3(i));
    S3(i+1)=S3(i)-dS3;%易感人数迭代
    E3(i+1)=E3(i)+rho3*c3*beta3*(1-q3)*S3(i)*(I3(i)+theta3*E3(i))-sigma*E3(i);%潜伏者人数迭代
    I3(i+1)=I3(i)+sigma*E3(i)-(deltaI+alpha3+gammaI)*I3(i);%感染者人数迭代
    Sq3(i+1)=Sq3(i)+rho3*c3*q3*(1-beta3)*S3(i)*(I3(i)+theta3*E3(i))-lambda*Sq3(i);%隔离易感染着人数迭代
    Eq3(i+1)=Eq3(i)+rho3*c3*beta3*q3*S3(i)*(I3(i)+theta3*E3(i))-deltaq*Eq3(i);%隔离潜伏者人数迭代
    H3(i+1)=H3(i)+deltaI*I3(i)+deltaq+Eq3(i)-(alpha3+gammaH)*H3(i);%住院患者人数迭代
    R3(i+1)=R3(i)+gammaI*I3(i)+gammaH*H3(i);%康复人数迭代 
    E_sum3=E_sum3+newinfected3(i+1);
    %if E_sum3/N0 > 1/100
       % rho3 = rho3/1.002;
       % q3=1*10^(-6)*1.01; %扩大隔离比例
      %  beta3=2.05*10^(-9)/1.005;%传染概率降低
   % end
end
plot(T,I1,'LineWidth',2);
hold on;
plot(T,I2,'LineWidth',2);
hold on;
plot(T,I3,'LineWidth',2);
grid on;
legend('乐观情景','中情景','悲观情景');
xlabel('日期');
ylabel('感染人数');