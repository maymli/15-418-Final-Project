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
for _ in range(2500000):
    first = random.randrange(0, 5000)
    second = random.randrange(0, 5000)
    if first != second:
        f1.write(str(first) + " " + str(second) + "\n")
f1.close()
'''

f1 = open("sparse-50000.txt", "w")
f1.write("50000\n")
for _ in range(5000000):
    first = random.randrange(0, 50000)
    second = random.randrange(0, 50000)
    if first != second:
        f1.write(str(first) + " " + str(second) + "\n")
f1.close()

'''
f2 = open("complete-8000.txt", "w")
f2.write("8000")
for i in range(8000):
    for j in range(i + 1, 8000):
        f2.write(str(i) + " " + str(j) + "\n")
f2.close()
'''
