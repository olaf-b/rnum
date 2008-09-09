require 'mybench'

a, b = bench_float
c = 0
puts "calculating c = a*b ..."

bench_time{ c = a * b }

