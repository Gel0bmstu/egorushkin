import os
import re
import matplotlib.pyplot as plt

# X = [[], [], []]
X = []
Z = [[], [], []]


with open('X_values.txt', 'r') as f:
    content = f.read().split('\t')
    try:
        for v in content:
            X.append(float(v))
    except Exception as e:
        pass

print(X[4800:4810])

plt.suptitle('Z', fontsize=20)
plt.plot(X)
# plt.plot(Z[1])
plt.grid()
plt.show()
