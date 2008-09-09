require 'mybench'

a = RNum::Vector::zeros(1000).indgen!
b = a.clone

puts "calculating dyad"
bench_time{ a.dyad(b) }

