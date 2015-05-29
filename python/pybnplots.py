
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_pdf import PdfPages 

def plot3d(x, y, z, 
           xlabel="x", ylabel="y", zlabel="z", 
           plotlabel="", 
           export=None):

    fsize=16
    mpl.rc('text', usetex=True)
    mpl.rc('font', size=fsize, family="times")

    # fig = plt.figure(figsize=plt.figaspect(1))
    fig = plt.figure(figsize=(5,4))
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    ax.plot_surface(x, y, z, 
                    rstride=1, cstride=1, 
                    shade=True,
                    cmap=mpl.cm.coolwarm,
                    edgecolor='none', 
                    antialiased=False)
    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)

    ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(nbins=5))
    ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(nbins=5))
    ax.zaxis.set_major_locator(mpl.ticker.MaxNLocator(nbins=5))

    # This. Is. The. Most. Perverse. Latex. Hack. I. Have. EVER. Written:
    ax.zaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda x, p: "\\rule[-10cm]{{7mm}}{{1pt}}${:5.1f}$".format(x)))
    # mplot3d (in the present version) does not allow controlling the
    # padding between the tick marks and the tick labels. The default
    # formatter placed the labels on top of the actual axes, resulting
    # in an unsightly layout of the figure. Now, the FuncFormatter
    # allows controlling the label output through a custom function
    # (the lambda above), so I figured I just add some spaces in front
    # of the number to effect some additional distance between the
    # axis and the number. But NOOO, that wouldn't work, since the
    # formatter strips ALL whitespace from the label. And I mean all!
    # I spent at least an hour -- AN HOUR -- googling, trying stuff,
    # to no avail. Finally, I came up with the above: Place a rule of
    # an appropriate length, AND DROP IT VERTICALLY SO FAR DOWN THAT
    # IT IS NOT VISIBLE IN THE PLOT.
    
    
    ax.xaxis._axinfo['label']['space_factor'] = 2
    ax.yaxis._axinfo['label']['space_factor'] = 2.4
    ax.zaxis._axinfo['label']['space_factor'] = 3
    
    ax.text2D(0.1, 0.9, plotlabel, 
              fontdict={'fontsize': fsize}, 
              transform=ax.transAxes)

    plt.subplots_adjust(left=-0.05, bottom=0.02, 
                        right=0.90, top=1.05)

    if export:
        pp = PdfPages(export)
        pp.savefig(fig)
        pp.close()
    else:
        plt.show()

    plt.close()


