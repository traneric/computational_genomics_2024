#!/usr/bin/env python3

from scipy.stats import binom_test

x = 2406
n = 268761
p = 0.010451

p_value = binom_test(x, n, p, alternative='two-sided')
p_value = format(p_value, ".2e")

print(p_value)