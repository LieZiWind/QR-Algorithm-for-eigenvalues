N=100;

A  =rand(N);
tic
A = Doublestep_QR(A,N,0.00001);
t = toc


%% First we upper Hessenberg a matrix A.  
%% To do so, we need householder transformation.


%% Next, we use implicit QR algorithm to compute U.
%% To do so, we need to pop the bulge
%% For that, a givens rotation is required.
%% Here we use single-step iteration and test on real eigs
%% Later we will try the double-step iteration and test on complex eigs
%% In fact we need to decide whether it is hessenberg or not!

%time = [];
%for i = 5:N
%    t=0;
%    for j = 1:2
%        A = rand(i);
%        tic
%        A = Doublestep_QR(A,i,0.0000001);
%        v = toc;
%        t = t + v;
%    end
%    t = t/2;
%    time = [time t];
%end
%plot(time)









function A = Upperhessenberg(A,n) % A the matrix to be dealt with, n the order of A
    for i = 1:n-1
        flag=0; % 0 means all zeros, 1 means find non-zero element
        %In fact, we find the first non-zero element other than the diag
        for j = i+1:n
            if A(j,i)~=0
                flag = 1;
                temp = Householder(A(j:n,i),n-j+1,n);
                break;
            end
        end
        if flag == 1
            A = temp*A*temp';
        end
    end
end

function H = Householder(v,l,M) % v the vector that need to be changed, l the length, M the order of the matrix
    v = v / norm(v,"inf"); %normalize
    
    h=zeros(1,l)';
   
    h(2:l) = v(2:l);   %we want the first element to be zero,so start from 2
    h(1) = v(1) + v(1)/abs(v(1))*norm(v);
    h=h/norm(h);
    H=zeros(M,M);


    if M+1-l>=1
        H(1:M-l,1:M-l)=eye(M-l);
    end
    H(M+1-l:M,M+1-l:M) = eye(l)-  2*(h *h');
end

function G = Givens(v,pos,M)% pos the place where zero is hoped, v & M the same as Householder
     G = eye(M);
     G(pos+2,pos+2)=0;
     G(pos+1,pos+1)=0;
    %We know that tan(theta) = v(2,1)/v(1,1)
    if v(1,1)==0
       G(pos+1,pos+1)=0;
       G(pos+1,pos+2)=1;
       G(pos+2,pos+1)=-1;
       G(pos+2,pos+2)=0;
    else 
        t = v(2,1)/v(1,1);
        G(pos+2,pos+2)= 1/sqrt(1+t*t);
        G(pos+1,pos+1) = G(pos+2,pos+2);
        %Make sure tan(theta) = v(2,1)/v(1,1) instead of -v(2,1)/v(1,1)
        if sign(t) ==1
            G(pos+1,pos+2)=sqrt(1-1/(t*t+1));
            G(pos+2,pos+1)=-G(pos+1,pos+2);
        else
            G(pos+1,pos+2)=-sqrt(1-1/(t*t+1));
            G(pos+2,pos+1)=-G(pos+1,pos+2);
        end
    end                                                
end

function A = Singlestep_Iteration(A,n,m)% A the matrix to be dealt with, n the order of A, m the order of a rather smaller A(tilde where outside is done)
    %The first givens rotation, introducing the bulge
 
    G = eye(n);
    G(n-m+1,n-m+1)=0;
    G(n-m+2,n-m+2)=0;
    norm = sqrt(A(n-m+1,n-m+1)^2 + A(n-m+2,n-m+1)^2);
    G(n-m+1,n-m+1) = A(n-m+1,n-m+1)/norm;
    G(n-m+2,n-m+2) =G(n-m+1,n-m+1);
    G(n-m+1,n-m+2) = A(n-m+2,n-m+1)/norm;
    G(n-m+2,n-m+1) = -G(n-m+1,n-m+2);
    

    A = G*A*G';
  
    %Pop the bulge using Givens rotation
    for i=n-m+1:n-2
        temp = Givens(A(i+1:i+2,i),i,n);
        A = temp*A*temp';
    end
end

function A = Singlestep_QR(A,n,tol)%A the matrix to be dealt with, n the order of A, tol the tolerence of accurary 
    A = Upperhessenberg(A,n);
    
    for i = 1:n-1
        while abs(A(i+1,i)) >= tol
            A = Singlestep_Iteration(A,n,n-i+1);
        end
    end
end

function A = Doublestep_Iteration(A,n,m,r,tol) %n the focusing matrix's size; B the whole matrix, m its size, r the upperleft corner  
    if n-r <=2 
        l = eye(m);
        [l(r:n,r:n) ,t] = schur(A(r:n,r:n));
        A = l'*A*l;

    else


    flag = 0 ;
    % Below we introduce the first bubble(Using Francis displacement):
    a = A(n,n);
    b = A(n-1,n-1);
    c = A(n-1,n);
    d = A(n,n-1);
    %sigma = solve('(x^2-(a+b)x-c*d+a*b=0','x');
    sigma = ((a+b)+sqrt((a+b)^2-4*(a*b-c*d)))/2;
    re = real(sigma);
    no = norm(sigma);

    H_begin = zeros(3,1);
    H_begin(1,1) = A(r,r)^2+A(r,r+1)*A(r+1,r)-2*re*A(r,r)+no^2;
    H_begin(2,1) = A(r+1,r)*(A(r,r)+A(r+1,r+1)-2*re);
    H_begin(3,1) = A(r+1,r) *A(r+2,r+1);



    Hnorm = norm(H_begin);

    H_begin(1,1) = H_begin(1,1)/Hnorm;
    H_begin(2,1) = H_begin(2,1)/Hnorm;
    H_begin(3,1) = H_begin(3,1)/Hnorm;

    H_rest = null(H_begin');
    H_small = [H_begin H_rest];

    H = zeros(m);
    H(r:r+2,r:r+2) = H_small;
    for i =r+3:m
        H(i,i)=1;
    end
    for i =1:r-1
        H(i,i)=1;
    end
    

    A = H'*A*H;

    % Next we chase the bulge

    for i = r+1:n-2
        if abs(norm(A(i+1:i+2,i-1))) >= tol 
            h_small = Householder(A(i:i+2,i-1),3,3);
            house = eye(m);
            house(i:i+2,i:i+2) = h_small;
            A = house'*A*house;
        else
            flag = 1;
            break
        end
    end
    
    if flag == 0
        % Finally a givens rotation
        g_small = Givens(A(n-1:n,n-2),n-2,n);
        g = eye(m);
        g(1:n,1:n) = g_small;
        A = g*A*g';
    end

    if flag == 1
        A = Doublestep_Iteration(A,n,m,i,tol);
    end
    end
end

function A = Doublestep_QR(A,n,tol) %A the matrix to be dealt with, n the order of A, tol the tolerence of accurary
    A = Upperhessenberg(A,n);
    i = n;
    r = 1;


    while i >= r+3
        flag =0;
        while abs(A(i-1,i-2)) > tol
            if abs(A(i,i-1)) > tol
                if abs(A(r+1,r)) >= tol
                    display(i)
                    display(r)
                    
                    A = Doublestep_Iteration(A,i,n,r,tol);

                else
                    r= r+1;
                    display(i)
                    display(r)
                    A = Doublestep_Iteration(A,i,n,r,tol);

                end
            else
                flag = 1;
                i = i-1;
                break
            end
        end
        if flag == 0
            i = i-2;
        end
    end


    U = eye(n);
    [U(i+1:i+2,i+1:i+2) ,t] = schur(A(i+1:i+2,i+1:i+2));
    A = U'*A*U;

    U = eye(n);
    [U(r+1:r+2,r+1:r+2) ,t] = schur(A(r+1:r+2,r+1:r+2));
    A = U'*A*U;

    U = eye(n);
    [U(1:3,1:3) ,t] = schur(A(1:3,1:3));
    A = U'*A*U;
    
    A = Clean(A,n,tol);
    

end

function A = Clean(A,n,tol)
    cleanup = 1;
    
    while cleanup <=n-1
        if abs(A(cleanup+1,cleanup)) <= tol
            cleanup = cleanup +1;
            continue
        end
        u = eye(n);
        [u(cleanup:cleanup+1,cleanup:cleanup+1) t] = schur(A(cleanup:cleanup+1,cleanup:cleanup+1));
       
        A = u'*A*u;
        
        cleanup = cleanup+2;
    end
end



