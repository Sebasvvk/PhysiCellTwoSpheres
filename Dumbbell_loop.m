clear all
hold off

%load data
outputfiles = dir('output0*_cells_physicell.mat'); 
iwant = cell((N),1);N = length(outputfiles);

for i = 1:N
    alloutput = load(outputfiles(i).name);
    iwant{i} = alloutput;
end

angleresult = zeros(size(iwant));
opts = optimset('Display','off');

for j = 1:(N)
    
    %make data
    allcells_xy =  transpose(iwant{j,:}.cells([2,3,6],:));
    xall = allcells_xy(:,1);
    yall = allcells_xy(:,2);

    %split data into two sets for the two phenotypes
    type02mask = or(allcells_xy(:,3) == 0, allcells_xy(:,3) == 2);
    type13mask = or(allcells_xy(:,3) == 1, allcells_xy(:,3) == 3);
    
    tf02 = any(type02mask,2);
    tf13 = any(type13mask,2);
    
    pheno02 = allcells_xy(tf02,:);
    pheno13 = allcells_xy(tf13,:);
    
    x02 = pheno02(:,1);
    y02 = pheno02(:,2);
    x13 = pheno13(:,1);
    y13 = pheno13(:,2);

    %boundary function
    boundaryall = boundary(xall,yall,0.8);
    xboundary = transpose(xall(boundaryall));
    yboundary = transpose(yall(boundaryall));

    %centers of mass
    x0 = sum(x02)./length(x02);
    y0 = sum(y02)./length(y02);
    x1 = sum(x13)./length(x13);
    y1 = sum(y13)./length(y13);
    
    com02 = [x0, y0];
    com13 = [x1, y1];

    %radius guess
    x02max = abs(max(x02)-x0);
    y02max = abs(max(y02)-y0);
    x13max = abs(max(x13)-x1);
    y13max = abs(max(y13)-y1);
    
    r0 = mean([x02max y02max]);
    r1 = mean([x13max y13max]);

    %fitting
    inputarguments = [x0,y0,r0,x1,y1,r1];
    fitcoeffs = lsqcurvefit(@dumbbell,inputarguments,xboundary,abs(yboundary),[],[],opts);
    
    angleresult(j) = angle(fitcoeffs,max(xboundary));
    
end

%construct data for plot
time = 0:1:(length(angleresult)-1);
cosa = transpose(cos(angleresult));
data = [time(1,:);cosa(1,:)];

%remove NaN values
datanonan = data(:,all(~isnan(data)));
angleresultnonan = datanonan(2,:);

cosfinalangle = mean(angleresultnonan((length(angleresultnonan)-20):length(angleresultnonan)));
xlin = linspace(min(datanonan(1,:)), max(datanonan(1,:)));
%plot the angle as a function of time
scatter(datanonan(1,:),datanonan(2,:),80, '.')
hold on

ft = fittype('exp(-(x-a)/b)',...
    'dependent',{'y'},'independent',{'x'},...
    'coefficients',{'a','b'});
f = fit(transpose(datanonan(1,:)),transpose(datanonan(2,:)),ft,'StartPoint',[500,50]);
coeffs = coeffvalues(f);
fitcurve = exp(-(xlin-coeffs(1))./coeffs(2));

xlim([min(datanonan(1,:))-20 max(datanonan(1,:))+20])
ylim([0 1])
xlabel('time (hours)','Interpreter','latex')
ylabel('cos $\theta$','Interpreter','latex')
plot(xlin,fitcurve,'black','LineWidth',2)
legend('Measurement points', 'Exponential fit','Interpreter','latex')
currentdirectory = pwd;
adhesion = currentdirectory(32:34);
speed = currentdirectory(38:40);
box on

%title(sprintf("Evolution of contact angle \\theta for adhesion "  + adhesion + " and speed " + speed))

%add line at final contact angle
%yline(cosfinalangle,'--r','\theta_f_i_n_a_l','LabelVerticalAlignment','bottom','LabelHorizontalAlignment','left','LineWidth',1.5)

%dumbbell function
function y = dumbbell(x,xdata)
    x0 = x(1);
    y0 = x(2);
    r0 = x(3);
    x1 = x(4);
    y1 = x(5);
    r1 = x(6);

    d = sqrt((x1-x0)^2+(y1-y0)^2);
    phi = asin((y1-y0)/d); 
    u = (r0^2-r1^2+d^2)/(2*d);
    
    angle0 = acos(u/r0);
    angle1 = acos((d-u)/r1);

    y = zeros(size(xdata));
    
    conditionarray = (xdata-x0)/r0;
    indexestoconsider1 = (-1 < conditionarray & conditionarray < 1);
    thetaall = acos(conditionarray);
    index1 = (thetaall-phi) > angle0;
    allindex1 = and(indexestoconsider1,index1);
    y(allindex1) = y0+sin(thetaall(allindex1))*r0;

    conditionarray = (xdata-x1)/r1;
    indexestoconsider2 = (-1 < conditionarray & conditionarray < 1);
    thetaall = acos(conditionarray);
    index2 = (pi-(thetaall-phi)) > angle1;
    allindex2 = and(indexestoconsider2,index2);
    y(allindex2) = y1+sin(thetaall(allindex2))*r1;
    
    y(y==0) = NaN;
    y = fillmissing(y,'linear');
end

%function to calculate angle from dumbbell
function angle0 = angle(x,width)
    x0 = x(1);
    y0 = x(2);
    r0 = x(3);
    x1 = x(4);
    y1 = x(5);
    r1 = x(6);
   
    d=sqrt((x1-x0)^2+(y1-y0)^2);
    u=(r0^2-r1^2+d^2)/(2*d);
    
    angle0 = NaN;
    
    if 0 <= u && u <= r0 
        angle0 = acos(u/r0);
    end 
    
    if u < 0 %this happens when the centers of mass move past each other
        angle0 = NaN;
    end 
    
    if ~isreal(angle0)
        angle0 = NaN;
    end 
    
    if (((x1+r1) <  (x0+r0)) && ((x0-r0) < (x1-r1))) || (((x0+r0) <  (x1+r1)) && ((x1-r1) < (x0-r0))) %gives a 0 when the two circles overlap completely
        angle0 = pi/2;
    end 
    
    if (((x0+u) > 0.5*width) || ((x0+u) < -0.5*width)) 
        angle0 = pi/2;
    end
end 
