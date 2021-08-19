---
title: "First post"
author: Matteo Murdaca
last_modified_at: 2021-07-19T11:08:00+1
categories:
  - Blog
tags:
  - Post Formats
  - readability
  - standard
---

This post has been updated and should show a modified date if used in a layout.

All children, except one, grow up. They soon know that they will grow up, and the way Wendy knew was this. One day when she was two years old she was playing in a garden, and she plucked another flower and ran with it to her mother. I suppose she must have looked rather delightful, for Mrs. Darling put her hand to her heart and cried, "Oh, why can't you remain like this for ever!" This was all that passed between them on the subject, but henceforth Wendy knew that she must grow up. You always know after you are two. Two is the beginning of the end.

```python
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# parameter values
u = 0
R0 = X
t_incubation = 5.1
t_infective = 3.3

# various initial numbers
e_initial = 0.00015
i_initial = 0.0015
r_initial = 0
v_initial = Y
s_initial = 1 - e_initial - i_initial - r_initial - (v_initial * 0.71)
d_initial = 0
o_initial = 0
icu_initial = 0

alpha = 1/t_incubation
gamma = 1/t_infective
beta = R0*gamma

# SEIR model differential equations.
def deriv(x, t, alpha, beta, gamma):
    s, e, i, r, d, o, v, icu = x
    dsdt = -(1-u-v_initial*0.71)*beta * s * i 
    dedt =  (1-u-v_initial*0.71)*beta * s * i - alpha * e
    didt = alpha * e - gamma * i 
    drdt =  gamma * i
    dddt = (gamma * i) * 0.005
    dodt = (alpha * e - gamma * i) * 0.1
    dvdt = v
    dicudt = (alpha * e - gamma * i) * 0.07
    return [dsdt, dedt, didt, drdt, dddt, dodt, dvdt, dicudt]

t = np.linspace(0, 160, 160)
x_initial = s_initial, e_initial, i_initial, r_initial, d_initial, o_initial, v_initial, icu_initial
soln = odeint(deriv, x_initial, t, args=(alpha, beta, gamma))
s, e, i, r, d, o, v, icu = soln.T

def plotdata(t, s, i, e, d, o, v, icu=None):
    # plot the data
    fig = plt.figure(figsize=(12,6))
    ax = [fig.add_subplot(221, axisbelow=True),
          fig.add_subplot(222, axisbelow=True), 
          fig.add_subplot(223, axisbelow=True),
          fig.add_subplot(224, axisbelow=True)]

    ax[0].plot(t, s, lw=3, label='Fraction Susceptible')
    ax[0].plot(t, i, lw=3, label='Fraction Infective')
    ax[0].plot(t, r, lw=3, label='Recovered')
    ax[0].plot(t, icu, lw=3, label='Intensive Care')
    ax[0].set_title('Susceptible and Recovered Populations')
    ax[0].set_xlabel('Time /days')
    ax[0].set_ylabel('Fraction')

    ax[3].plot(t, o, lw=3, label='Hospitalized')
    ax[3].plot(t, d, lw=3, label='Dead')
    ax[3].plot(t, icu, lw=3, label='Intensive Care')
    ax[3].set_title('Susceptible and Recovered Populations')
    ax[3].set_xlabel('Time /days')
    ax[3].set_ylabel('Fraction')

    ax[1].plot(t, i, lw=3, label='Infective')
    ax[1].set_title('Infectious Population')
    if e is not None: ax[1].plot(t, e, lw=3, label='Exposed')
    ax[1].set_ylim(0, 0.07)
    ax[1].set_xlabel('Time /days')
    ax[1].set_ylabel('Fraction')

    ax[2].plot(t, d, lw=3, label='Fraction Dead')
    ax[2].plot(t, o, lw=3, label='Fraction Hospitalized')
    ax[2].plot(t, i, lw=3, label='Fraction Infective')
    ax[2].set_title('Dead and Hospitalized population')
    ax[2].set_xlabel('Time /days')
    ax[2].set_ylabel('Fraction')


    for a in ax: 
        a.grid(True)
        a.legend()

    plt.tight_layout()

plotdata(t, s, e, i, d, o, v, icu)
```