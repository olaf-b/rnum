% mesuring performance of Matlab
n = 500;
m = linspace(0,n*n-1,n*n);
m = rem(m, n+1) + 1.0;
m = reshape(m,n,n);
x = reshape(linspace(1,n*n,n*n),n,n)+0.0;

disp(sprintf('solving %ix%i matrix...' ,n ,n));
t1 = cputime ();

y = m \ x;

t2 = cputime ();
disp(sprintf('Time: %5.4f sec.', t2-t1))
