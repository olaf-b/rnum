require 'mybench'

n = 500

a = RNum::zeros([n,n]).indgen!(1).rem(n+1).add!(1)
b = RNum::zeros([n,n]).indgen!(1).rem(n-1).add!(1)
c = 0
printf "executing %ix%i Matrix product...\n",n,n
puts 'c=a*b'

bench_time(1) { c=a*b }

# time: 4.21 sec

