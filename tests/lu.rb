require 'mybench'
n = 1000
m = RNum::zeros([n,n]).indgen!.rem!(n+1).add!(1)
m.t!

printf "executing %ix%i LU...\n",n,n
bench_time(1){m.getrf}
