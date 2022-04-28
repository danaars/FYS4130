import numpy as np

L = 20
j = 12

counter = 0
for _ in range(j+1, L):
    print(_)
    counter += 1

print(f"L = {L}, j = {j}\nCounter = {counter}")
print("søren klybe")
'''
np.random.seed(1337+420+69)

L = 16
beta = 0
J = 1
#p = 1 - np.exp( - 2*beta * J)
p = 0.5

chain = np.random.randint(0, 2, L)

def flip(arr,site, f, l, flips):
    if arr[site] == 0:
        arr[site] = 1
    elif arr[site] == 1:
        arr[site] = 0

    flips.append(site)

    right = 1
    left = -1

    if arr[(site + right) % l] != arr[site]:# and site != f:
        nextsite = (site + right) % l
        if np.random.rand() > p:
            flip(arr, nextsite, site, l, flips)

    if arr[(site + left) % l] != arr[site]:# and site != f:
        nextsite = (site + left) % l
        if np.random.rand() > p:
            flip(arr, nextsite, site, l, flips)

print("-1 mod 5")
print(-1%5)
'''
'''
print(chain)

for i in range(10):
    flips = []
    s = np.random.randint(0, len(chain))
    flip(chain, s, None, len(chain), flips)

    print(flips)
    print(chain)
'''
