import csv
import json
import math
from collections import OrderedDict
from math import ceil

import matplotlib.pyplot as plt
from scipy.stats import norm

round_ = lambda x: round(x, 2)
csv_out = csv.writer(open('out_main.csv', 'w'))
data = json.load(open('input1.txt'))
n = len(data)
x_max, x_min = max(data), min(data)
m = 12
h = (x_max - x_min) / m
inters = [0] * m

for x in data:
    i = ceil((x - x_min) / (h)) or 1
    inters[i - 1] += 1

results = OrderedDict([
    ('i', ['Интервалы']),
    ('bs', ['Границы интервалов']),
    ('x', ['xᵢ']),
    ('n', ['nᵢ']),
    ('p', ['pᵢ']),
    ('m', ['xᵢ*pᵢ']),
    ('d', ['xᵢ²*pᵢ'])
])
En = 0
Ep = 0
X_ = 0
m2 = 0

for i in range(m):
    results['i'].append(i + 1)
    results['bs'].append(
        '[{}, {}{}'.format(round_(x_min + i * h), round_(x_min + (i + 1) * h), ')' if i != n - 1 else ']'))
    x_i = x_min + (i + 0.5) * h
    results['x'].append(round_(x_i))
    n_i = inters[i]
    En += n_i
    results['n'].append(n_i)
    p = n_i / n
    Ep += p
    results['p'].append(round_(p))
    m_i = x_i * p
    X_ += m_i
    results['m'].append(round_(m_i))
    d_i = x_i ** 2 * p
    m2 += d_i
    results['d'].append(round_(d_i))

results['m'].append('X̅ = {}'.format(round_(X_)))
results['d'].append('m₂ = {}'.format(round_(m2)))

plt.bar(results['i'][1:], results['p'][1:])

plt.show()

x_cur = x_min
F = 0
for i, p in zip(results['i'][1:], results['p'][1:]):
    plt.plot([x_cur, x_cur], [F, F + p], '--k')
    plt.plot([x_cur], [F + p], 'ok', mfc='none')
    F += p
    plt.plot([x_cur, x_min + i * h], [F, F], '-k')
    plt.plot([x_min + i * h], [F], 'ok')
    x_cur += h

plt.show()

S_2 = m2 - X_ ** 2
S = S_2 ** 0.5
print('S = {}'.format(round_(S)))

t = 1.95
e = t * S / n ** 0.5
print('{} < m < {}'.format(round_(X_ - e), round_(X_ + e)))

# 5
new_iters = results['n'][1:]

iters_left = 0
while new_iters[0] < 8:
    new_iters = [new_iters.pop(0) + new_iters.pop(0)] + new_iters
    iters_left += 1

iters_right = 0
while new_iters[-1] < 8:
    new_iters = new_iters + [new_iters.pop() + new_iters.pop()]
    iters_right += 1

new_xis = [x_min + ((iters_left + 1)) * h] + [xi + h / 2 for xi in results['x'][iters_left + 2:-iters_right - 1]] + [
    math.inf]

r = len(new_iters) - 3
xi_1 = None


def ksii(ni, xi):
    global xi_1
    zi = (xi - X_) / S
    if xi_1 is not None:
        zi_1 = (xi_1 - X_) / S
    pi = norm.cdf(zi) - (norm.cdf(zi_1) if xi_1 is not None else 0)
    xi_1 = xi
    return ni ** 2 / (n * pi)


ksi2 = sum([ksii(ni, xi) for ni, xi in zip(new_iters, new_xis)]) - n

ksirp = float(input('Введите χᵣₚ² для r={} и p=0.05: '.format(r)))
print('χ² = {} < χᵣₚ² - {}'.format(round(ksi2, 2), ksi2 < ksirp))

results['bs'].append('')
results['x'].append('')
results['i'].append('Доп. резы')
results['p'].append('Σpᵢ = {}'.format(round_(Ep)))
results['n'].append('Σnᵢ = {}'.format(round_(En)))

for val in results.values():
    csv_out.writerow(val)

print('Таблицу с интервалами вы можете увидеть в файле out_main.csv')
