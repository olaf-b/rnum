require 'rnum'

REPEAT = 100
ARRSZ = 1_000_000

def bench_array
	[ RNum::Vector::zeros(ARRSZ).indgen!, 
	  RNum::Vector::zeros(ARRSZ).indgen! ]
end

def bench_float
	 bench_array
end

def bench_time(n=REPEAT)
	t1 = Process.times.utime
	for i in 1..n
		yield
	end
	t2 = Process.times.utime
	puts " Time: %.2f sec\n\n" % [t2-t1]
end
