%%
load fjarrvarme89.dat
load fjarrvarme90.dat

pow89 = fjarrvarme89(:,2);

air89 = fjarrvarme89(:,3);

wat89 = fjarrvarme89(:,4);

%%
model_data_length = 1350;
val = (80+1350+1):(80+1350+2*7*24);

y_mean = mean(pow89(80:80+model_data_length));
y = pow89(80:80+model_data_length) - y_mean;

a_mean = mean(air89(80:80+model_data_length));
a = air89(80:80+model_data_length) - a_mean;

w_mean = mean(wat89(80:80+model_data_length));
w = wat89(80:80+model_data_length) - w_mean;

val = (80+1350+1):(80+1350+2*7*24);
y_val = pow89(val) - y_mean;
a_val = air89(val) - a_mean;
w_val = wat89(val) - w_mean;

y_test1 = pow89(1767:1934) - y_mean; 
y_test2 = pow89(3440:3607) - y_mean;

a_test1 = air89(1767:1934) - a_mean;
a_test2 = air89(3440:3607) - a_mean;

w_test1 = wat89(1767:1934) - w_mean;
w_test2 = wat89(3440:3607) - w_mean;

%% New Start
model_data_length = 1350;
%a = air89(80:80+model_data_length) - mean(air89(50:50+model_data_length));
plotter(a, 80)

%% Trying outlier detection and replacement
a_initial = a;
%clear a
a  = filloutliers(a_initial, 'spline', 'quartile');
figure
plot(1:length(a), a_initial, 1:length(a), a)
plotter(a)
%% WITH OUTLIERS REPLACED
a_dat = iddata(a);
am_init = idpoly([1 zeros(1, 26)],[],[1 zeros(1, 24)]);
%am_init.Structure.a.Free = [0 1 1 1 zeros(1,19) 1 1 1 1];
%am_init.Structure.c.Free = [0 1 1 1 zeros(1, 19) 1 1];

am_init.Structure.a.Free = [0 1 1 1 1 zeros(1,17) 1 1 1 1 1];
am_init.Structure.c.Free = [0 1 1 1 zeros(1, 19) 1 1];


am = pem(a_dat, am_init);
a_pw = filter(am.a, am.c, a);
a_pw = a_pw(length(am.c)+1:end);
plotter(a_pw, 80, 0)
lbpTest(a_pw)
%%
%y = pow89(80:80+model_data_length) - mean(pow89(80:80+model_data_length));
y_pw = filter(am.a, am.c, y);
y_pw = y_pw(length(am.c):end);

plotter(y_pw)

figure
crosscorr(a_pw, y_pw, 'NumLags', 70)
hold on
M=30;
plot(-M:M, 2/sqrt(length(y_pw))*ones(1,  2*M+1),'--')
plot(-M:M, -2/sqrt(length(y_pw))*ones(1,  2*M+1),'--')
hold off

%%
A2 = [1]; % r
B = [0 1];  % d zeros + s
z_dat = iddata(y,a); 
M1_init = idpoly([],B,[],[],A2);
M1 = pem(z_dat, M1_init);

v_hat = resid(M1, z_dat);
plotter(v_hat.y)

%% NOT ADDRESSING POTENTIAL OUTLIERS
plotter(a, 80)
a_dat = iddata(a);
am_init = idpoly([1 zeros(1, 25)],[],[1 zeros(1, 23)]);
am_init.Structure.a.Free = [0 1 1 1 1 1 1 1 zeros(1,13) 1 1 0 1 1];
am_init.Structure.c.Free = [0 0 0 1 zeros(1, 19) 1];

am = pem(a_dat, am_init);
a_pw = filter(am.a, am.c, a);
a_pw = a_pw(length(am.c)+1:end);
plotter(a_pw, 80, 0)
lbpTest(a_pw, 80)
%%
figure
subplot(211)
plot(a)
subplot(212)
plot(y)
%%
%y = pow89(80:80+model_data_length) - mean(pow89(80:80+model_data_length));
y_pw = filter(am.a, am.c, y);
y_pw = y_pw(length(am.c):end);

plotter(y_pw)

figure
crosscorr(a_pw, y_pw, 'NumLags', 80)
hold on
M=30;
plot(-M:M, 2/sqrt(length(y_pw))*ones(1,  2*M+1),'--')
plot(-M:M, -2/sqrt(length(y_pw))*ones(1,  2*M+1),'--')
hold off

%% (d,r,s) = (1,0,0) or (1,1,0) or (1,2,0)
A2 = [1 ]; % r
B = [0 1];  % d zeros + s
z_dat = iddata(y,a); 
M1_init = idpoly([],B,[],[],A2);
M1 = pem(z_dat, M1_init);

v_hat = resid(M1, z_dat);
plotter(v_hat.y, 50)
%% v_hat modeled as AR(2)
A1 = [1 zeros(1,26)];
C1 = [1 zeros(1,24)];

v_hat = y - filter(M1.b, M1.f, a);
v_hat = v_hat(length(M1.f):end);
v_dat = iddata(v_hat);
%% diff
v_diff = filter([1 zeros(1, 23) -1], [1],v_hat);
v_diff = v_diff(24:end);
plotter(v_diff, 50)
%%
M2_init = idpoly(A1,[],C1);
M2_init.Structure.a.Free = [0 1 1 1 0 0 0 1 zeros(1,15) 1 1 1 1];
M2_init.Structure.c.Free = [0 1 1 zeros(1,20) 1 1];
M2 = pem(v_dat, M2_init);


res_v = resid(M2, v_dat);
plotter(res_v.y, 80, 0);
[h, p] = kstest((res_v.y-mean(res_v.y)/std(res_v.y)))
lbpTest(res_v.y)
%% everything estimated together


%bj_init = idpoly([], [0 1], [1 zeros(1,24)], [1 zeros(1,26)], [1]);
bj_init = idpoly([], [0 1], [1], [1 zeros(1,26)], [1]);

%bj_init.Structure.c.Free = [0 1 zeros(1,21) 1 1];
%bj_init.Structure.d.Free = [0 1 1 1 zeros(1,20) 1 1 1];
%bj_init.Structure.c.Free = [0 1 1 zeros(1,20) 1 1];
%bj_init.Structure.d.Free = [0 1 1 1 0 0 0 1 zeros(1,15) 1 1 1 1];

bj = pem(z_dat, bj_init);

res_bj = resid(bj, z_dat);
plotter(res_bj.y, 80, 0)
lbpTest(res_bj.y, 80)

%% Prediction
% u = exogenous input i.e. air temp. y = power
A = bj.d;
B = conv(bj.d, bj.b);
C = bj.c;

k = 6;
[Fy_1, Gy_1] = polydiv(C, A, k);

[Fu_1, Gu_1] = polydiv(conv(B, Fy_1), C, k);

uhat_1 = filter(Gu_1, C, a);

yhat0_1 = filter(Gy_1, C, y);

yhat_1 = uhat_1 + yhat0_1;
%yhat_3 = yhat0_3;
out = k+length(C);%max(k,length(C))+1;

errors_1 = yhat_1(out:end) - y(out:end);
var(errors_1)

figure
set(gcf, 'Position',  [500, 500, 1000, 400])

subplot(121)
plot(y(out:end))
hold on 
plot(yhat_1(out:end))
hold off

subplot(122)
plot(errors_1)

%% Prediction on validation dataset
a_val = air89(val) - mean(air89(80:80+model_data_length));
y_val = pow89(val) - mean(pow89(80:80+model_data_length));
%%
A = bj.d;
B = conv(bj.d, bj.b);
C = bj.c;

k = 6;
[Fy_k, Gy_k] = polydiv(C, A, k);
[Fu_k, Gu_k] = polydiv(conv(B, Fy_k), C, k);

uhat_k = filter(Gu_k, C, [a; a_val]);
yhat0_k = filter(Gy_k, C, [y; y_val]);

%uhat_k = uhat_k(end-length(a_val)+1:end);
uhat_k  =  uhat_k( length(y)-k+1 : length([y;y_val])-k );
yhat0_k = yhat0_k( length(y)-k+1 : length([y;y_val])-k );

yhat_k = uhat_k + yhat0_k;% + mean(pow89(80:80+model_data_length))+ mean(air89(50:50+model_data_length));


%errors_k = yhat_k(k+1:end) - y_val(1:end-k);
errors_k = yhat_k(k+2:end) - y_val(1:end-k-1);

var(errors_k)

figure
set(gcf, 'Position',  [500, 500, 1000, 400])

subplot(121)
plot(y_val(1:end-k-1))
%plot(y_val)
hold on 
plot(yhat_k(k+2:end))
%plot(yhat_k)

hold off

subplot(122)
plot(errors_k)
hold on
bound_k = 1.96*sqrt(var(res_bj.y)*sqrt(sum(Fy_k.^2)));
M = 335;
plot(1:M,  bound_k*ones(1,  M),'--')
plot(1:M, -bound_k*ones(1,  M),'--')
%%
m = mean(pow89(80:80+model_data_length));
figure
plot(80:(80+length(y)-1), y, (80+length(y)):(80+length(y)+length(y_val)-1), y_val, (80+length(y)):(80+length(y)+length(y_val)-1) , yhat_k)
hold on
plot(80:(80+length(pow89(80:80+1350+2*7*24))-1) , pow89(80:80+1350+2*7*24))
hold off
%% AIR
%% Prediction on test1 data set
% 1767:1934
y_temp = pow89(1350:1766) - y_mean;
a_temp = air89(1350:1766) - a_mean;
w_temp = wat89(1350:1766) - w_mean;

A = bj.d;
B = conv(bj.d, bj.b);
C = bj.c;

k = 6;
[Fy_k, Gy_k] = polydiv(C, A, k);
[Fu_k, Gu_k] = polydiv(conv(B, Fy_k), C, k);

uhat_k = filter(Gu_k, C, [a_temp; a_test1]);
yhat0_k = filter(Gy_k, C, [y_temp; y_test1]);

%uhat_k = uhat_k(end-length(a_val)+1:end);
uhat_k  =  uhat_k( length(a_temp)-k+1 : length([a_temp;a_test1])-k );
yhat0_k = yhat0_k( length(y_temp)-k+1 : length([y_temp;y_test1])-k );

yhat_k = uhat_k + yhat0_k;% + mean(pow89(80:80+model_data_length))+ mean(air89(50:50+model_data_length));


%errors_k = yhat_k(k+1:end) - y_val(1:end-k);
errors_k = yhat_k(k+2:end) - y_test1(1:end-k-1);

var(errors_k)

figure
set(gcf, 'Position',  [500, 500, 1000, 400])

subplot(121)
plot(y_test1(1:end-k-1))
%plot(y_val)
hold on 
plot(yhat_k(k+2:end))
%plot(yhat_k)

hold off

subplot(122)
plot(errors_k)
hold on
bound_k = 1.96*sqrt(var(res_bj.y)*sqrt(sum(Fy_k.^2)));
M = 335;
plot(1:M,  bound_k*ones(1,  M),'--')
plot(1:M, -bound_k*ones(1,  M),'--')

%% Prediction on test2 data set
% 1767:1934
y_temp = pow89(1350:1766) - y_mean;
a_temp = air89(1350:1766) - a_mean;
w_temp = wat89(1350:1766) - w_mean;

A = bj.d;
B = conv(bj.d, bj.b);
C = bj.c;

k = 6;
[Fy_k, Gy_k] = polydiv(C, A, k);
[Fu_k, Gu_k] = polydiv(conv(B, Fy_k), C, k);

uhat_k = filter(Gu_k, C, [a_temp; a_test2]);
yhat0_k = filter(Gy_k, C, [y_temp; y_test2]);

%uhat_k = uhat_k(end-length(a_val)+1:end);
uhat_k  =  uhat_k( length(a_temp)-k+1 : length([a_temp;a_test2])-k );
yhat0_k = yhat0_k( length(y_temp)-k+1 : length([y_temp;y_test2])-k );

yhat_k = uhat_k + yhat0_k;% + mean(pow89(80:80+model_data_length))+ mean(air89(50:50+model_data_length));


%errors_k = yhat_k(k+1:end) - y_val(1:end-k);
errors_k = yhat_k(k+2:end) - y_test2(1:end-k-1);

var(errors_k)

figure
set(gcf, 'Position',  [500, 500, 1000, 400])

subplot(121)
plot(y_test2(1:end-k-1))
%plot(y_val)
hold on 
plot(yhat_k(k+2:end))
%plot(yhat_k)

hold off

subplot(122)
plot(errors_k)
hold on
bound_k = 1.96*sqrt(var(res_bj.y)*sqrt(sum(Fy_k.^2)));
M = 335;
plot(1:M,  bound_k*ones(1,  M),'--')
plot(1:M, -bound_k*ones(1,  M),'--')
%% B 
w = fjarrvarme89(80:80+model_data_length,4) - mean(fjarrvarme89(80:80+model_data_length,4));

crosscorr(y, a, 'NumLags', 80)
figure
crosscorr(y, w, 'NumLags', 80)

%% Outlier stuff
plotter(w)
w = filloutliers(w, 'spline','quartiles');
plotter(w)
%%
plotter(a)
a = filloutliers(a, 'spline', 'quartiles');
plotter(a)
%%
plotter(y)
y = filloutliers(y, 'spline', 'quartiles');
plotter(y)
%% Model for water temp
plotter(w, 50)
w_data = iddata(w);
w_init = idpoly([1 0 0 0],[],[1 0 0]);
%w_init.Structure.a.Free = [0 1 0 1];
w_model = pem(w_data, w_init);

w_pw = filter(w_model.a, w_model.c, w);
w_pw = w_pw(length(w_model.c):end);
plotter(w_pw, 80, 0)


%% pw new y
%y_new = y - filter(bj.b, [1], a);
y_new = y - filter(M1.b, M1.f, a);
y_new_pw = filter(w_model.a, w_model.c, y_new);

figure
crosscorr(w_pw,y_new_pw, 'NumLags', 80);
%% (d,r,s) = (0,0,0)
z_wat_dat = iddata(y_new, w);
M2_init = idpoly([],[0 0 1],[],[],[1]);
M2 = pem(z_wat_dat, M2_init);
v_hat = resid(M2, z_wat_dat);

plotter(v_hat.y, 80, 1)
%%
v_hat = y_new - filter(M2.b, M2.f, w);
v_hat = v_hat(length(M2.f):end);

plotter(v_hat, 80)
%%
v_diff = filter([1 zeros(1, 23) -1], [1], v_hat);
plotter(v_diff, 80)

m = idpoly([])
%%
v_data = iddata(v_hat);
v_hat_init = idpoly([1 zeros(1,26)], [],[1 zeros(1,24)]);
%v_hat_init = idpoly([1 zeros(1,26)], [],[1 zeros(1,24)]);
v_hat_init.Structure.a.Free = [0 1 1 1 1 1 1 zeros(1, 14) 1 1 1 1 1 1];
v_hat_init.Structure.c.Free = [0 1 zeros(1, 21) 1 1 ];
%v_hat_init.Structure.a.Free = [0  1 1 1  1 1 1  1 1 1  1 1 1  1 1 1  1 1 1  1 1 1  1 1 1  1 1 1  1];
%v_hat_init.Structure.c.Free = [0  1 1 1  1 1 1  1 1 1  1 1 1  1 1 1  1 1 1  1 1 1  1 1 1  1 1 ];

v_hat_model = pem(v_data, v_hat_init);

v_res = resid(v_hat_model, v_data);
plotter(v_res.y, 80, 0)
lbpTest(v_res.y, 80)

%% 
bj_final_init = idpoly([1], {[0 1], [0 0 1]}, [1 zeros(1, 24)], [1 zeros(1, 26)], {[1], [1]});

%bj_final_init.Structure.d.Free = [0  1 1 zeros(1, 20) 0 1 1 1];
%bj_final_init.Structure.c.Free = [0 zeros(1, 23) 1 ];

%Current
%bj_final_init.Structure.d.Free = [0 1 1 1 1 1 1 zeros(1, 14) 1 1 1 1 1 1];
%bj_final_init.Structure.c.Free = [0 1 zeros(1, 21) 1 1 ];

bj_final_init.Structure.d.Free = [0 1 1 zeros(1, 21) 1 1 1];
bj_final_init.Structure.c.Free = [0 zeros(1, 23) 1 ];


%bj_init = idpoly([], [0 1], [1 zeros(1,23)], [1 zeros(1,25)], [1]);
%bj_init.Structure.c.Free = [0 1 zeros(1,21) 1];
%bj_init.Structure.d.Free = [0 1 1 1 zeros(1,19) 1 1 1];

final_data = iddata(y,[a w]);
bj_final = pem(final_data, bj_final_init);

plotter(resid(bj_final, final_data).y, 80, 0)

%% Prediction with air and water temps as inputs

A = bj_final.d;
B1 = conv(cell2mat(bj_final.b(:,1)), bj_final.d);
B2 = conv(cell2mat(bj_final.b(:,2)), bj_final.d);
C = bj_final.c;

% Set prediction step
k = 6;

[Fy_k, Gy_k] = polydiv(C, A, k);

[F_a_k, G_a_k] = polydiv(conv(B1, Fy_k), C, k);
[F_w_k, G_w_k] = polydiv(conv(B2, Fy_k), C, k);

x
a_hat_k = filter(G_a_k, C, a);
w_hat_k = filter(G_w_k, C, w);


yhat0_k = filter(Gy_k, C, y);

yhat_k = a_hat_k + w_hat_k + yhat0_k;


out = 30;

errors_k = yhat_k(out:end) - y(out:end);
var(errors_k)

figure
set(gcf, 'Position',  [500, 500, 1000, 400])

subplot(121)
plot(y(out:end))
hold on 
plot(yhat_k(out:end))
hold off

subplot(122)
plot(errors_k)

%% Prediction on validation dataset

% 3440:3607 

A = bj_final.d;
B1 = conv(cell2mat(bj_final.b(:,1)), bj_final.d);
B2 = conv(cell2mat(bj_final.b(:,2)), bj_final.d);
C = bj_final.c;

% Set prediction step
k = 6;

[Fy_k, Gy_k] = polydiv(C, A, k);

[F_a_k, G_a_k] = polydiv(conv(B1, Fy_k), C, k);
[F_w_k, G_w_k] = polydiv(conv(B2, Fy_k), C, k);
%

a_hat_k = filter(G_a_k, C, [a; a_val]);
w_hat_k = filter(G_w_k, C, [w; w_val]);

%a_hat_k = a_hat_k(end-length(a_val):end);
%w_hat_k = w_hat_k(end-length(w_val):end);

a_hat_k = a_hat_k( length(y)-k+1 : length([y;y_val])-k );
w_hat_k = w_hat_k( length(y)-k+1 : length([y;y_val])-k );


yhat0_k = filter(Gy_k, C, [y; y_val]);
yhat0_k = yhat0_k( length(y)-k+1 : length([y;y_val])-k );

yhat_k = a_hat_k + w_hat_k + yhat0_k;


%out = 30;

%errors_k = yhat_k(out:end) - y_val(out:end);
errors_k = yhat_k(k+2:end) - y_val(1:end-k-1);
var(errors_k)

figure
set(gcf, 'Position',  [200, 200, 1000, 400])

subplot(121)
plot(y_val(1:end-k-1))
hold on 
plot(yhat_k(k+2:end))
hold off

subplot(122)
plot(errors_k)
hold on

bound_k = 1.96*sqrt(var(resid(bj_final, final_data).y)*sqrt(sum(Fy_k.^2)));
M = 335;
plot(1:M,  bound_k*ones(1,  M),'--')
plot(1:M, -bound_k*ones(1,  M),'--')
hold off

%% Prediction on test1 data set
% 1767:1934
y_temp = pow89(1350:1766) - y_mean;
a_temp = air89(1350:1766) - a_mean;
w_temp = wat89(1350:1766) - w_mean;

A = bj_final.d;
B1 = conv(cell2mat(bj_final.b(:,1)), bj_final.d);
B2 = conv(cell2mat(bj_final.b(:,2)), bj_final.d);
C = bj_final.c;

% Set prediction step
k = 1;

[Fy_k, Gy_k] = polydiv(C, A, k);

[F_a_k, G_a_k] = polydiv(conv(B1, Fy_k), C, k);
[F_w_k, G_w_k] = polydiv(conv(B2, Fy_k), C, k);
%

a_hat_k = filter(G_a_k, C, [a_temp; a_test1]);
w_hat_k = filter(G_w_k, C, [w_temp; w_test1]);

%a_hat_k = a_hat_k(end-length(a_val):end);
%w_hat_k = w_hat_k(end-length(w_val):end);

a_hat_k = a_hat_k( length(a_temp)-k+1 : length([a_temp;a_test1])-k );
w_hat_k = w_hat_k( length(w_temp)-k+1 : length([w_temp;w_test1])-k );


yhat0_k = filter(Gy_k, C, [y_temp; y_test1]);
yhat0_k = yhat0_k( length(y_temp)-k+1 : length([y_temp;y_test1])-k );

yhat_k = a_hat_k + w_hat_k + yhat0_k;


%out = 30;

%errors_k = yhat_k(out:end) - y_val(out:end);
errors_k = yhat_k(k+2:end) - y_test1(1:end-k-1);
var(errors_k)

figure
set(gcf, 'Position',  [500, 500, 1000, 400])

subplot(121)
plot(y_test1(1:end-k-1))
hold on 
plot(yhat_k(k+2:end))
hold off

subplot(122)
plot(errors_k)
hold on

bound_k = 1.96*sqrt(var(res_bj.y)*sqrt(sum(Fy_k.^2)));
M = 335;
plot(1:M,  bound_k*ones(1,  M),'--')
plot(1:M, -bound_k*ones(1,  M),'--')
hold off

%% Prediction on test2 data set

y_temp = pow89(3000:3439) - y_mean;
a_temp = air89(3000:3439) - a_mean;
w_temp = wat89(3000:3439) - w_mean;

A = bj_final.d;
B1 = conv(cell2mat(bj_final.b(:,1)), bj_final.d);
B2 = conv(cell2mat(bj_final.b(:,2)), bj_final.d);
C = bj_final.c;

% Set prediction step
k = 1;

[Fy_k, Gy_k] = polydiv(C, A, k);

[F_a_k, G_a_k] = polydiv(conv(B1, Fy_k), C, k);
[F_w_k, G_w_k] = polydiv(conv(B2, Fy_k), C, k);
%

a_hat_k = filter(G_a_k, C, [a_temp; a_test2]);
w_hat_k = filter(G_w_k, C, [w_temp; w_test2]);

%a_hat_k = a_hat_k(end-length(a_val):end);
%w_hat_k = w_hat_k(end-length(w_val):end);

a_hat_k = a_hat_k( length(y_temp)-k+1 : length([y_temp;y_test2])-k );
w_hat_k = w_hat_k( length(y_temp)-k+1 : length([y_temp;y_test2])-k );


yhat0_k = filter(Gy_k, C, [y_temp; y_test2]);
yhat0_k = yhat0_k( length(y_temp)-k+1 : length([y_temp;y_test2])-k );

yhat_k = a_hat_k + w_hat_k + yhat0_k;


%out = 30;

%errors_k = yhat_k(out:end) - y_val(out:end);
errors_k = yhat_k(k+2:end) - y_test2(1:end-k-1);
var(errors_k)

figure
set(gcf, 'Position',  [500, 500, 1000, 400])

subplot(121)
plot(y_test2(1:end-k-1))
hold on 
plot(yhat_k(k+2:end))
hold off

subplot(122)
plot(errors_k)
hold on

bound_k = 1.96*sqrt(var(res_bj.y)*sqrt(sum(Fy_k.^2)));
M = 335;
plot(1:M,  bound_k*ones(1,  M),'--')
plot(1:M, -bound_k*ones(1,  M),'--')
hold off

%%
A = bj_final.d;
B1 = conv(cell2mat(bj_final.b(:,1)), bj_final.d);
B2 = conv(cell2mat(bj_final.b(:,2)), bj_final.d);
C = bj_final.c;
%% RECURSIVE ESTIMATION with 41

% state sp eqs
sigma_e = 10^(-5);
A = eye(41);
Re = sigma_e*eye(41);
%Re = [10^(-4) 0; 0 0]; 
Rw = 1.25;


% initial vals
Rxx_1 = 0.001*eye(41);
a_coeff_init = [bj_final.d(2:7), bj_final.d(22:27)];
c_coeff_init = [bj_final.c(2), bj_final.c(24:25)];
xtt_1 = [a_coeff_init, bj_final.b{1}(2)*[1, a_coeff_init], bj_final.b{2}(3)*[1, a_coeff_init], c_coeff_init]';
 
% where we store vals
N = length(y);
e = randn(N,1);
xsave = zeros(41,N-28);
% kalman filter

for k=29:N
    % bj_final_init.Structure.d.Free = [0 1 1 1 1 1 1 zeros(1, 14) 1 1 1 1 1 1];
    % bj_final_init.Structure.c.Free = [0 1 zeros(1, 21) 1 1 ]
    y_part = -flip([y(k-26:k-21)', y(k-6:k-1)']);
    a_part =  flip([a(k-27:k-22)', a(k-7:k-1)']);
    w_part =  flip([w(k-28:k-23)', w(k-8:k-2)']);
    e_part =  flip([e(k-24:k-23)', e(k-1)']);
    
    C = [y_part, a_part, w_part, e_part];
        
    % update
    Ryy = C * Rxx_1 * C' + Rw;
    Kt = Rxx_1*C'*inv(Ryy);
    xtt = xtt_1 + Kt*(y(k) - C*xtt_1);
    Rxx = Rxx_1 - Kt * Ryy * Kt';
    
    % save
    xsave(:,k-28) = xtt;
    
    % predict
    Rxx_1 = A * Rxx * A' + Re;
    xtt_1 = A*xtt;
end;

figure
plot(xsave')
%% RECURSIVE ESTIMATION with 18

% state sp eqs
sigma_e = 10^(-5);
A = eye(18);
Re = sigma_e*eye(18);
%Re = [10^(-4) 0; 0 0]; 
Rw = 3.7;


% initial vals
Rxx_1 = 0.001*eye(18);
a_coeff_init = [bj_final.d(2:3), bj_final.d(25:27)];
c_coeff_init = bj_final.c(25);
xtt_1 = [a_coeff_init, bj_final.b{1}(2)*[1, a_coeff_init], bj_final.b{2}(3)*[1, a_coeff_init], c_coeff_init]';
 
% where we store vals
N = length(y);
e = sqrt(Rw)*randn(N,1);
xsave = zeros(18,N-28);
% kalman filter

for k=29:N
    %bj_final_init.Structure.d.Free = [0 1 1 zeros(1, 21) 1 1 1];
    %bj_final_init.Structure.c.Free = [0 zeros(1, 23) 1 ];

    y_part = -flip([y(k-26:k-24)', y(k-2:k-1)']);
    a_part =  flip([a(k-27:k-25)', a(k-3:k-1)']);
    w_part =  flip([w(k-28:k-26)', w(k-4:k-2)']);
    e_part =  e(k-24);
    
    C = [y_part, a_part, w_part, e_part];
        
    % update
    Ryy = C * Rxx_1 * C' + Rw;
    Kt = Rxx_1*C'*inv(Ryy);
    xtt = xtt_1 + Kt*(y(k) - C*xtt_1);
    Rxx = Rxx_1 - Kt * Ryy * Kt';
    
    % save
    xsave(:,k-28) = xtt;
    
    % predict
    Rxx_1 = A * Rxx * A' + Re;
    xtt_1 = A*xtt;
end;



figure()
set(gcf, 'Position',  [200, 200, 1600, 500])
subplot(131)
plot(xsave(1,:)')
hold on
plot(xsave(2,:)')
hold on
plot(xsave(3,:)')
hold on
plot(xsave(4,:)')
hold on
plot(xsave(5,:)')
hold on
plot(xsave(18,:)')
legend('a1','a2', 'a24', 'a25', 'a26', 'c24')

subplot(132)
plot(xsave(6,:)')
hold on
plot(xsave(7,:)')
hold on
plot(xsave(8,:)')
hold on
plot(xsave(9,:)')
hold on
plot(xsave(10,:)')
hold on
plot(xsave(11,:)')

legend('b1,0','b1,1', 'b1,24', 'b1,25', 'b1,26')

subplot(133)
plot(xsave(12,:)') 
hold on
plot(xsave(13,:)') 
hold on
plot(xsave(14,:)')
hold on
plot(xsave(15,:)')
hold on
plot(xsave(16,:)')
hold on
plot(xsave(17,:)')
hold on

legend('b2,0','b2,1', 'b2,24', 'b2,25', 'b2,26')

%% RECURSIVE ESTIMATION with 18 + prediction 6 steps

y_kal = pow89;
a_kal = air89;
w_kal = wat89;


% state sp eqs
sigma_e = 10^(-5);
A = eye(18);
Re = sigma_e*eye(18);
Rw = 3.7;

e_kal = sqrt(Rw)*randn(length(y_kal),1);


% initial vals
Rxx_1 = 0.001*eye(18);
a_coeff_init = [bj_final.d(2:3), bj_final.d(25:27)];
c_coeff_init = bj_final.c(25);
xtt_1 = [a_coeff_init, bj_final.b{1}(2)*[1, a_coeff_init], bj_final.b{2}(3)*[1, a_coeff_init], c_coeff_init]';
 
% where we store vals
N = length(pow89);
xsave = zeros(18,length(y_kal)-28-6);
ypredict = zeros(1, length(y_kal)-28-6);
% kalman filter

for k=29:N-6
    %bj_final_init.Structure.d.Free = [0 1 1 zeros(1, 21) 1 1 1];
    %bj_final_init.Structure.c.Free = [0 zeros(1, 23) 1 ];

    y_part = -flip([y_kal(k-26:k-24)', y_kal(k-2:k-1)']);
    a_part =  flip([a_kal(k-27:k-25)', a_kal(k-3:k-1)']);
    w_part =  flip([w_kal(k-28:k-26)', w_kal(k-4:k-2)']);
    e_part =  e_kal(k-24);
    
    C = [y_part, a_part, w_part, e_part];
        
    % update
    Ryy = C * Rxx_1 * C' + Rw;
    Kt = Rxx_1*C'*inv(Ryy);
    xtt = xtt_1 + Kt*(y_kal(k) - C*xtt_1);
    Rxx = Rxx_1 - Kt * Ryy * Kt';
    
    % save
    xsave(:,k-28) = xtt;
    
    % predict
    Rxx_1 = A * Rxx * A' + Re;
    xtt_1 = A*xtt;
    
    C1 = [-flip([y_kal(k-25:k-23)', y_kal(k-1:k)']), flip([a_kal(k-26:k-24)', a_kal(k-2:k)']), flip([w_kal(k-27:k-25)', w_kal(k-3:k-1)']), e_kal(k-23)];
    y_1 = C1*xtt_1;
    
    C2 = [y_1, -flip([y_kal(k-24:k-22)', y_kal(k)]), flip([a_kal(k-25:k-23)', a_kal(k-1:k+1)']), flip([w_kal(k-26:k-24)', w_kal(k-2:k)']), e_kal(k-22)];
    y_2 = C2*xtt_1;
    
    C3 = [y_2, y_1, -flip(y_kal(k-23:k-21)'), flip([a_kal(k-24:k-22)', a_kal(k:k+2)']), flip([w_kal(k-25:k-23)', w_kal(k-1:k+1)']), e_kal(k-21)];
    y_3 = C3*xtt_1;
    
    C4 = [y_3, y_2, -flip(y_kal(k-22:k-20)'), flip([a_kal(k-23:k-21)', a_kal(k+1:k+3)']), flip([w_kal(k-24:k-22)', w_kal(k:k+2)']), e_kal(k-20)];
    y_4 = C4*xtt_1;
    
    C5 = [y_4, y_3, -flip(y_kal(k-21:k-19)'), flip([a_kal(k-22:k-20)', a_kal(k+2:k+4)']), flip([w_kal(k-23:k-21)', w_kal(k+1:k+3)']), e_kal(k-19)];
    y_5 = C5*xtt_1;
    
    C6 = [y_5, y_4, -flip(y_kal(k-20:k-18)'), flip([a_kal(k-21:k-19)', a_kal(k+3:k+5)']), flip([w_kal(k-22:k-20)', w_kal(k+2:k+4)']), e_kal(k-18)];
    y_6 = C6*xtt_1;
    
    y_predict(1, k-28) = y_6;
        
end;
%y_predict = y_predict(N+1:end);

%% Plotting

figure()
set(gcf, 'Position',  [50, 50, 1500, 500])
subplot(131)
plot(xsave(1,:)')
hold on
plot(xsave(2,:)')
hold on
plot(xsave(3,:)')
hold on
plot(xsave(4,:)')
hold on
plot(xsave(5,:)')
hold on
plot(xsave(18,:)')
legend('a1','a2', 'a24', 'a25', 'a26', 'c24')

subplot(132)
plot(xsave(6,:)')
hold on
plot(xsave(7,:)')
hold on
plot(xsave(8,:)')
hold on
plot(xsave(9,:)')
hold on
plot(xsave(10,:)')
hold on
plot(xsave(11,:)')

legend('b1,0','b1,1', 'b1,24', 'b1,25', 'b1,26')

subplot(133)
plot(xsave(12,:)') 
hold on
plot(xsave(13,:)') 
hold on
plot(xsave(14,:)')
hold on
plot(xsave(15,:)')
hold on
plot(xsave(16,:)')
hold on
plot(xsave(17,:)')
hold on

legend('b2,0','b2,1', 'b2,24', 'b2,25', 'b2,26')
figure
plot(1:length(y_predict), y_val(1:end-28-6), 1:length(y_predict), y_predict)
figure
plot(y_val)
hold on
plot(y_predict)