import random

# Test 1: complete-5000
# |V| = 5,000, |E| = 12,497,500
f1 = open("complete-5000.txt", "w")
f1.write("5000\n")
for i in range(5000):
    for j in range(i + 1, 5000):
        f1.write(str(i) + " " + str(j) + "\n")
f1.close()


# Test 2: random-5000
# |V| = 5,000, |E| = 2,500,000
f1 = open("random-5000.txt", "w")
f1.write("5000\n")
pairs_count = 0
while pairs_count < 2500000:
    first = random.randrange(0, 5000)
    second = random.randrange(0, 5000)
    if first != second:
        f1.write(str(first) + " " + str(second) + "\n")
        pairs_count += 1
f1.close()

# Test 3: random-50000
# |V| = 50,000, |E| = 250,000,000
f1 = open("random-50000.txt", "w")
f1.write("50000\n")
pairs_count = 0
while pairs_count < 250000000:
    first = random.randrange(0, 50000)
    second = random.randrange(0, 50000)
    if first != second:
        f1.write(str(first) + " " + str(second) + "\n")
        pairs_count += 1
f1.close()

# Test 4: sparse-50000
# |V| = 50,000, |E| = 1,000,000
f1 = open("sparse-50000.txt", "w")
f1.write("50000\n")
pairs_count = 0
while pairs_count < 1000000:
    first = random.randrange(0, 50000)
    second = random.randrange(0, 50000)
    if first != second:
        f1.write(str(first) + " " + str(second) + "\n")
        pairs_count += 1
f1.close()

# Test 5: sparse-500000
# |V| = 500,000, |E| = 10,000,000
f1 = open("sparse-500000.txt", "w")
f1.write("500000\n")
pairs_count = 0
while pairs_count < 10000000:
    first = random.randrange(0, 500000)
    second = random.randrange(0, 500000)
    if first != second:
        f1.write(str(first) + " " + str(second) + "\n")
        pairs_count += 1
f1.close()

# Test 6: very-sparse-50000
# |V| = 50,000, |E| = 100,000
f1 = open("very-sparse-50000.txt", "w")
f1.write("50000\n")
pairs_count = 0
while pairs_count < 100000:
    first = random.randrange(0, 50000)
    second = random.randrange(0, 50000)
    if first != second:
        f1.write(str(first) + " " + str(second) + "\n")
        pairs_count += 1
f1.close()

# Test 7: corner-50000
# |V| = 50,000, |E| = 2,500,000
f1 = open("corner-50000.txt", "w")
f1.write("50000\n")
# vertices 0 to 1000 are in a complete graph
for i in range(1000):
    for j in range(i + 1, 1000):
        f1.write(str(i) + " " + str(j) + "\n")
# remaining vertices are sparse
pairs_count = 499500
while pairs_count < 2500000:
    first = random.randrange(0, 50000)
    second = random.randrange(0, 50000)
    if first != second:
        f1.write(str(first) + " " + str(second) + "\n")
        pairs_count += 1
f1.close()

# Test 8: components-5000
# |V| = 5,000, |E| = 1,000,000 (5 highly-connected components of 200,000 edges)
f1 = open("components-5000.txt", "w")
f1.write("5000\n")
for k in range(5):
    start = 1000 * k
    end = 1000 * (k + 1)
    pairs_count = 0
    while pairs_count < 200000:
        first = random.randrange(start, end)
        second = random.randrange(start, end)
        if first != second:
            f1.write(str(first) + " " + str(second) + "\n")
            pairs_count += 1
f1.close()

# Test 9: components-50000
# |V| = 50,000, |E| = 10,000,000 (5 highly-connected components of 2,000,000 edges)
f1 = open("components-50000.txt", "w")
f1.write("50000\n")
for k in range(5):
    start = 10000 * k
    end = 10000 * (k + 1)
    pairs_count = 0
    while pairs_count < 2000000:
        first = random.randrange(start, end)
        second = random.randrange(start, end)
        if first != second:
            f1.write(str(first) + " " + str(second) + "\n")
            pairs_count += 1
f1.close()
