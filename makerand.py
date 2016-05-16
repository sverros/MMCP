import numpy as np

g1 = open('rand1.txt', 'w')
g2 = open('rand2.txt', 'w')
g3 = open('rand3.txt', 'w')
g4 = open('rand4.txt', 'w')

numpts = 86000
randnums = np.random.randn(10000*numpts)

for i in range(0, 2500*numpts):
    s = str(randnums[i])
    g1.write(s + '\n')

for i in range(2500*numpts, 5000*numpts):
    s = str(randnums[i])
    g2.write(s + '\n')

for i in range(5000*numpts, 7500*numpts):
    s = str(randnums[i])
    g3.write(s + '\n')

for i in range(7500*numpts, 10000*numpts):
    s = str(randnums[i])
    g4.write(s + '\n')




g1.close()
g2.close()
g3.close()
g4.close()
