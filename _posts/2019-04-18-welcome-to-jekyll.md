---
title: "Welcome to Jekyll!"
date: 2019-04-18T15:34:30-04:00
categories:
  - blog
tags:
  - Jekyll
  - update
---

You'll find this post in your `_posts` directory. Go ahead and edit it and re-build the site to see your changes. You can rebuild the site in many different ways, but the most common way is to run `jekyll serve`, which launches a web server and auto-regenerates your site when a file is updated.

To add new posts, simply add a file in the `_posts` directory that follows the convention `YYYY-MM-DD-name-of-post.ext` and includes the necessary front matter. Take a look at the source for this post to get an idea about how it works.

Jekyll also offers powerful support for code snippets:

```ruby
def print_hi(name)
  puts "Hi, #{name}"
end
print_hi('Tom')
#=> prints 'Hi, Tom' to STDOUT.
```
Check out the [Jekyll docs][jekyll-docs] for more info on how to get the most out of Jekyll. File all bugs/feature requests at [Jekyll’s GitHub repo][jekyll-gh]. If you have questions, you can ask them on [Jekyll Talk][jekyll-talk].

[jekyll-docs]: https://jekyllrb.com/docs/home
[jekyll-gh]:   https://github.com/jekyll/jekyll
[jekyll-talk]: https://talk.jekyllrb.com/

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