import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button

"""
This function is used to plot and compare the different orders from an image.
The plot starts by plotting the 0th order, other orders can be shown by changing 
the slider. The gives an interactive idea of how the orders change staring from the 
first order (left on the flatfield) and ending at the last order (right on the flatfield).
"""
def plotOrdersWithSlider(orders, xValues=None, xMax=10560, yMax=30000, xMin=0):
    # The parametrized function to be plotted
    def f(t, m):
        m = int(m)
        mth_order = orders[m]
        return mth_order

    def t(m):
        m = int(m)
        if xValues is None:
            mth_order = orders[m]
            return np.arange(len(mth_order))
        else:
            return xValues[m]


    # Define initial parameters
    init_order = 0

    # Create the figure and the line that we will manipulate
    fig, ax = plt.subplots()
    line, = plt.plot(t(init_order), f(t(init_order), init_order), lw=2)
    ax.set_xlabel('x Position [pxl]')
    plt.xlim([xMin, xMax])
    plt.ylim([0, yMax])
    
    # adjust the main plot to make room for the sliders
    plt.subplots_adjust(left=0.25)

    # Make a vertically oriented slider to control the amplitude
    axamp = plt.axes([0.1, 0.25, 0.0225, 0.63])
    amp_slider = Slider(
        ax=axamp,
        label="Amplitude",
        valmin=1,
        valmax=len(orders),
        valinit=init_order,
        orientation="vertical"
    )


    # The function to be called anytime a slider's value changes
    def update(val):
        line.set_ydata(f(t(amp_slider.val-1), amp_slider.val-1))
        line.set_xdata(t(amp_slider.val-1))
        fig.canvas.draw_idle()        


    # register the update function with each slider
    amp_slider.on_changed(update)

    plt.show()


def plotGIF(figures, background):

    for i in np.arange(5):
        plt.imshow(background, origin='lower')
        plt.imshow(figures[i], alpha=0.1, origin='lower')
        plt.show()

