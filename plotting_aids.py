import numpy as np
import matplotlib.pyplot as plt
import scipy.stats

def plot_linfit(x, y, color=[1, 0.4, 0.4], alpha=0.2):
    # fit a curve to the data using a least squares 1st order polynomial fit
    z = np.polyfit(x,y,1)
    p = np.poly1d(z)
    fit = p(x)
    # get the coordinates for the fit curve
    c_y = [np.min(fit),np.max(fit)]
    c_x = [np.min(x),np.max(x)]
    
    # predict y values of origional data using the fit
    p_y = z[0] * x + z[1]

    # calculate the y-error (residuals)
    y_err = y -p_y

    # create series of new test x-values to predict for
    p_x = np.linspace(np.min(x), np.max(x), num=50)
    
    # now calculate confidence intervals for new test x-series
    mean_x = np.mean(x)  # mean of x
    n = len(x)  # number of samples in origional fit
    t = scipy.stats.t.isf(0.05 / 2, len(x)-1)  # t value (where n=9, two tailed 95%)
    s_err = np.sum(y_err ** 2)  # sum of the squares of the residuals

    confs = t * np.sqrt(
        (s_err / (n - 2)) * (1.0 / n + (p_x - x.mean())**2 / (((x**2).sum()) - n*(x.mean()**2)))
    )

    # now predict y based on test x-values
    p_y = z[0]*p_x+z[1]

    # get lower and upper confidence limits based on predicted y and confidence intervals
    lower = p_y - abs(confs)
    upper = p_y + abs(confs)

    # plot line of best fit
    plt.plot(p_x, p_y, '-', color=color,
             label='Regression line')
    # fill the confidence intervals
    plt.fill_between(p_x, lower, upper, color=color, alpha=alpha)