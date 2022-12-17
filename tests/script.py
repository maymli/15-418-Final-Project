import random

'''
f1 = open("complete-5000.txt", "w")
f1.write("5000\n")

for i in range(5000):
    for j in range(i + 1, 5000):
        f1.write(str(i) + " " + str(j) + "\n")
f1.close()
'''

'''
f1 = open("random-5000.txt", "w")
f1.write("5000\n")
#pairs_count = 0
#while pairs_count < 2500000:
for _ in range(2500000):
    first = random.randrange(0, 5000)
    second = random.randrange(0, 5000)
    if first != second:
        f1.write(str(first) + " " + str(second) + "\n")
        # pairs_count += 1
f1.close()
'''
'''
f1 = open("random-50000.txt", "w")
f1.write("50000\n")
#pairs_count = 0
#while pairs_count < 2500000:
for _ in range(250000000):
    first = random.randrange(0, 50000)
    second = random.randrange(0, 50000)
    if first != second:
        f1.write(str(first) + " " + str(second) + "\n")
        # pairs_count += 1
f1.close()
'''
'''
f1 = open("sparse-50000.txt", "w")
f1.write("50000\n")
for _ in range(1000000):
    first = random.randrange(0, 5000)
    second = random.randrange(0, 5000)
    if first != second:
        f1.write(str(first) + " " + str(second) + "\n")
f1.close()
'''
'''
f1 = open("sparse-50000.txt", "w")
f1.write("50000\n")
for _ in range(2500000):
    first = random.randrange(0, 5000)
    second = random.randrange(0, 5000)
    if first != second:
        f1.write(str(first) + " " + str(second) + "\n")
f1.close()
'''
'''
f1 = open("very-sparse-50000.txt", "w")
f1.write("50000\n")
for _ in range(25000):
    first = random.randrange(0, 5000)
    second = random.randrange(0, 5000)
    if first != second:
        f1.write(str(first) + " " + str(second) + "\n")
f1.close()
'''

'''
f2 = open("complete-50000.txt", "w")
f2.write("50000")
for i in range(8000):
    for j in range(i + 1, 8000):
        f2.write(str(i) + " " + str(j) + "\n")
f2.close()
'''

'''
f1 = open("corner-50000.txt", "w")
f1.write("50000\n")

# vertices 0 to 1000 are in a complete graph
for i in range(1000):
    for j in range(i + 1, 1000):
        f1.write(str(i) + " " + str(j) + "\n")

# remaining vertices are sparse
for _ in range(2500000):
    first = random.randrange(0, 50000)
    second = random.randrange(0, 50000)
    if first != second:
        f1.write(str(first) + " " + str(second) + "\n")
f1.close()
'''

'''
# n-cycle graph (also sparse)
f1 = open("n-cycle-50000.txt", "w")
f1.write("50000\n")

for i in range(50000 - 1):
    f1.write(str(i) + " " + str(i + 1) + "\n")
f1.write(str(50000 - 1) + " " + str(0) + "\n")
f1.close()
'''

f1 = open("components-5000.txt", "w")
f1.write("5000\n")

# vertices 0 to 1000 are in a complete graph
for k in range(5):
    start = 1000 * k
    end = 1000 * (k + 1)
    for _ in range(0, 200000):
        first = random.randrange(start, end)
        second = random.randrange(start, end)
        if first != second:
            f1.write(str(first) + " " + str(second) + "\n")
f1.close()

