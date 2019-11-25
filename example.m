% Before running this example, you need to unpack the package robot-10.2.zip.
% And put all the files in the search path to make sure they can be found.

clear all
mdl_puma560

%set the numerical test object
a1 = 72;
a2 = 42;
a3 = -44;
a4 = 52;
a5 = 2;
a6 = 106;
T0 = p560.fkine([a1 a2 a3 a4 a5 a6]/180*pi);

%sovle all solutions of the inverse kinematics
k = 40;
Q = zeros(k,6);
Q1 = zeros(k,6);
tic%timing
 for j = 1:k
     for n = 1:6
         Q(j,n) = Q(j,n) + normrnd(0,pi/2);
     end
     Q1(j,:) = iksovle(p560,T0,Q(j,:))/pi*180;
 end
Q2 = roundn(Q1(:,1:6),-4);
Q3 = rmmissing(Q2);
unique(Q3,'rows') 
toc
