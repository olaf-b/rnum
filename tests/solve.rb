require 'rnum'
n = 500
m = RNum::zeros([n,n]).indgen!.rem!(n+1).add!(1).t
x = RNum::zeros([n,n]).indgen!.t
printf "solving  y=x/m ...\n",n,n
#y = x/m
t1 = Process.times.utime
x.solve!(m)
#y = x/m
t2 = Process.times.utime
puts "Time: %.4f sec\n\n" % [t2-t1]
