require 'mybench'

a, b = bench_float
c = 0

print "calculating c = a+b ...\n"

#bench_time{ a.add!(b) }
bench_time{ c = a + b }

