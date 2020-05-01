clear
R=0
for(e=1)              //e stands for energy level, use to select energies that need to be plotted
    //create grid
    maxX = 10;                      //bounds of the problem
    delta = 0.2;
    x = -maxX:delta:maxX;
    n = length(x);                         //define number of points
    x = x';
    E = e;      //Energy
    
    //potential energy
    l=2;
    A=0;
    B=10;
    k=sqrt(2*E);
    for(i=1:n)
        xi(i)=-exp(2*%pi*x(i)/l);               //xi value for Eckart potential
        V(i)=-A*xi(i)/(1-xi(i))-B*xi(i)/((1-xi(i))^2);     //Potential energy function
        //V(i)=(0.5-exp(x(i)))^2;
    end
    //V=zeros(n,1);
    
    //determine Q matrix
    for(j=1:n)
        w(j) = 2*(V(j)-E);
    end
    Q(n) =exp(%i*k*delta);                     //First input Q 
    for(l=n:-1:2)
        Q(l-1) = (w(l)*delta^2+2-Q(l))^(-1);              //Next Q point
    end
    
    //Reflectance and Transmission
    R(e) = (exp(%i*k*x(2))-Q(1)*exp(%i*k*x(1)))/(Q(1)*exp(-%i*k*x(1))-exp(-%i*k*x(2)))
    disp(R)
    Refl(e) = abs(R(e))^2
    Trans(e) = 1-Refl(e);
    
    //wavefunction and propagation
    u(1)=exp(-%i*k*x(1)) + R(e)*exp(%i*k*x(1));
    for(s=1:n-1)
        u(s+1) = Q(s)*u(s);                             //wavefunction
    end
    ur = real(u);
    ui = imag(u);
    scf(0);clf();
    plot(x,real(u));     //Real part in blue
    plot(x,imag(u),'k'); //imag part in black
    plot(x,V,'r');  //potential in red
    title('Tunneling through Eckart barrier','fontsize',7)
    xlabel('X direction','fontsize',6)
    ylabel('Amplitude','fontsize',6)
    hl=legend("Real part","Imaginary part","Potential")
    hl.font_size=4
    //plot(x,real(Q),'c');
    //plot(x,imag(Q),'g');
end
