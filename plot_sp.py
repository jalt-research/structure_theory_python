from machine import *

"""
    SP lattice sizes and enumeration times
    for random small finite state machines
"""

symbols = ['a','b','c','d','e','f','g','h','i','j','k','l']
# You don't need the rest of the alphabet where we're going

(start, end) = (2, 9)
I = 200
avg_times_naive = []
avg_times_algorithm = []
avg_counts_naive = []
avg_counts_algorithm = []
max_a = []
max_t = []

for s in range(start,end):
    print(symbols[:s])
    times_naive = []
    times_algorithm = []
    count_algorithm = []
    count_naive = []
    for i in range(I):
        M3 = Machine.random(symbols[:s],['0','1'])
        sp = 0

        s1 = time.time()
        P = [p for p in Partition.list(M3.S) if M3.is_SP(p)]
        s2 = time.time()
        count_naive += [len(P)]
        
        s3 = time.time()
        P1 = M3.enumerate_SP(True)
        s4 = time.time()
        count_algorithm += [len(P1)]
        
        times_naive += [s2 - s1 + 0.0]
        times_algorithm += [s4 - s3 + 0.0]
    avg_times_naive += [sum(times_naive)/I]
    avg_times_algorithm += [sum(times_algorithm)/I]
    avg_counts_naive += [sum(count_naive)/I]
    avg_counts_algorithm += [sum(count_algorithm)/I]
    max_a += [max(count_algorithm)]
    max_t += [max(times_algorithm)]
print('time, naive method \n\t', avg_times_naive)
print('time, efficient method \n\t', avg_times_algorithm)
print('SP partitions \n\t', avg_counts_naive, '\n (check) \n', avg_counts_algorithm, ')')
print('largest SP lattice \n\t', max_a)
print('time on largest \n\t', max_t)
    
plot(avg_times_algorithm, label='Average running time (Hartmanis and Stearns)')
plot(avg_times_naive, label='Average running time (Naive)')
plot(max_t, label='Longest running time (H&S)')
legend()
show()
plot(avg_counts_algorithm, label='Average SP lattice size')
plot(max_a, label='Largest SP lattice seen')
legend()
show()

