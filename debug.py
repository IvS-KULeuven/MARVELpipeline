import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons

"""
This function is used to plot and compare the different orders from an image.
The plot starts by plotting the 0th order, other orders can be shown by changing
the slider. The gives an interactive idea of how the orders change staring from the
first order (left on the flatfield) and ending at the last order (right on the flatfield).
"""



def plotOrdersWithSlider(orders, xValues=None, xMax=10560, yMax=30000, xMin=0, nOrders=5, title=None):

    # Define initial parameters
    init_order = 33
    init_fiber = 1

    # The function to be called anytime a slider's value changes
    fiber = init_fiber
    order = init_order
    fiber_value = {"fiber 1": 1, "fiber 2": 2, "fiber 3": 3, "fiber 4": 4, "fiber 5": 5}

    # The parametrized function to be plotted
    def f(o, f):
        idx = (int(o)-33)*nOrders + f - 1
        spectrum = orders[idx]
        return spectrum

    def t(o, f):
        idx = (int(o)-33)*nOrders + f - 1
        if xValues is None:
            spectrum = orders[idx]
            return np.arange(len(spectrum))
        else:
            return xValues[idx]




    # Create the figure and the line that we will manipulate
    fig, ax = plt.subplots()
    line, = plt.plot(t(init_order, init_fiber), f(init_order, init_fiber), lw=2)
    ax.set_xlabel('x Position [pxl]')
    plt.xlim([xMin, xMax])
    plt.ylim([0, yMax])

    if not title is None:
            plt.title(title)

    # adjust the main plot to make room for the sliders
    plt.subplots_adjust(left=0.25)

    # Make radio boxes
    rax = fig.add_axes([0.05, 0.7, 0.15, 0.15])
    radio = RadioButtons(rax, ("fiber 1", "fiber 2", "fiber 3", "fiber 4", "fiber 5"))

    # Make a vertically oriented slider to control the amplitude
    axamp = fig.add_axes([0.1, 0.13, 0.0225, 0.5])
    amp_slider = Slider(
        ax=axamp,
        label="Order",
        valmin=33,
        valmax=int(len(orders)/nOrders)+33,
        valstep=np.arange(33, int(len(orders)/nOrders+33)),
        valinit=init_order,
        orientation="vertical"
    )




    def update(val):
        draw_image()

    def update_fiber(label):
        draw_image()

    def draw_image():
        fiber = fiber_value[radio.value_selected]
        order = amp_slider.val

        line.set_ydata(f(order, fiber))
        line.set_xdata(t(order, fiber))

        fig.canvas.draw_idle()



    # register the update function with each slider
    amp_slider.on_changed(update)
    radio.on_clicked(update_fiber)

    plt.show()


def plotGIF(figures, background):

    for i in np.arange(5):
        plt.imshow(background, origin='lower')
        plt.imshow(figures[i], alpha=0.1, origin='lower')
        plt.show()

