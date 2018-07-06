function dydt = PSIIODES (t,y,k) % 'ODE' for 'Ordinary Differential Equation'
 
dydt = zeros(31,1); 

if rem(t/50,2) <1
    h = 1;
%     dummy = 0.01;
else
    h = 0;
%     dummy = -0.01;
end
  
r(1) = h*k(1)*y(1);
r(2) = k(2)*y(2)*y(3);
r(3) = k(3)*y(4)*y(5);
r(4) = k(4)*y(6)*y(7);
r(5) = k(5)*y(4)*y(8); %0
r(6) = k(6)*y(8)*y(9);
r(7) = k(7)*y(8)*y(10)*y(11);
r(8) = k(8)*y(11)*y(12);
r(9) = k(9)*y(13);
r(10) = k(10)*y(31)*(y(21))^2;
r(11) = h*k(11)*y(33);
r(12) = k(12)*y(15)*y(16);
r(13) = k(13)*y(17)*y(18);
r(14) = k(14)*y(19)*y(20);
r(15) = k(15)*y(20)*y(22);
r(16) = k(16)*y(22)*y(23);
r(17) = k(17)*y(25)*y(11)*(y(24))^2;
r(18) = k(18)*y(15);
r(19) = k(19)*y(2);
r(20) = k(20)*y(26);
r(21) = k(21)*y(24); %y(11) 
r(22) = k(22)*y(27)*y(30);
r(23) = k(23)*y(29)*(y(28))^3*(y(26))^2;
r(24) = k(24)*y(32);
r(25) = k(25)*y(28);
r(26) = k(26)*y(30);
r(27) = k(27)*y(14);
r(28) = k(28)*y(31);
r(29) = k(29)*y(26);
r(30) = 0; %r(21)b k(30)*y(24);
r(31) = h*k(31);


dydt(1) = -r(1) + r(3) + r(5) + r(19); 
dydt(2) = r(1) - r(2) - r(19);
dydt(3) = -r(2) + r(4);
dydt(4) = r(2) - r(3) - r(5); 
dydt(5) = -r(3) + r(31); 
dydt(6) = r(2) - r(4); 
dydt(7) = -r(4) +r(5) + r(6) +r(7); 
dydt(8) = + r(4) - r(5) - r(6) -r(7); 
dydt(9) = -r(6) + r(24); 
dydt(10)= r(6) - r(7); 
dydt(11)= -r(7) - r(8) - r(17) + r(20) - r(21) + r(22) + r(23) + r(26); 
dydt(12)= r(7) - r(8); 
dydt(13)= r(8) - r(9); 
dydt(14)= 0;%-r(27); 
dydt(15)= r(11) - r(12) - r(18);
dydt(16)= -r(12) + r(13);
dydt(17)= r(12) - r(13);
dydt(18)= -r(13) + r(15) + r(16); 
dydt(19)= 2*r(10) - r(14);
dydt(20)= r(12) - r(14) - r(15);
dydt(21)= -2*r(10) + r(14);
dydt(22)= r(13) - r(15) - r(16);
dydt(23)= -r(16) + 2*r(17) + r(21) - r(30);
dydt(24)= r(16) - 2*r(17) - r(21) + r(30);
dydt(25)= -r(17) + r(20) + 2*r(23) - r(29);
dydt(26)= r(17) - r(20) - 2*r(23) + r(29); 
dydt(27)= -r(22) + 3*r(23) + r(25); 
dydt(28)= r(22) - 3*r(23) - r(25);
dydt(29)= r(27) - r(23);
dydt(30)= 2*r(10) - r(22) - r(26)+ r(31);
dydt(31)= r(9) - r(10) - r(28);
dydt(32)= -r(24) + r(28) +r(10);
dydt(33)= -r(11) + r(14) + r(15) + r(18); 


end 