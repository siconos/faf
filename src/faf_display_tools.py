def setAxLinesBW(ax):
    """
    Take each Line2D in the axes, ax, and convert the line style to be
    suitable for black and white viewing.
    """
    MARKERSIZE = 5

    COLORMAP = {
        'b': {'marker': None, 'dash': (None,None)},
        'g': {'marker': None, 'dash': [5,5]},
        'r': {'marker': None, 'dash': [5,3,1,3]},
        'c': {'marker': None, 'dash': [1,3]},
        'm': {'marker': None, 'dash': [5,2,5,2,5,10]},
        'y': {'marker': None, 'dash': [5,3,1,2,1,10]},
        'k': {'marker': 'o',  'dash': [1,2,1,10]}
        }

    MARKER = [None, 'o', '+', '*']

    lines_to_adjust = ax.get_lines()
    print(lines_to_adjust)
    count =1
    for line in lines_to_adjust:
        print(len(lines_to_adjust))
        origColor = line.get_color()
        line.set_color('black')
        line.set_dashes(COLORMAP[origColor]['dash'])
        print(count/7)
        line.set_marker(MARKER[count/7])
        count = count +1
        line.set_markersize(MARKERSIZE)

    try:
        lines_legend_to_adjust = ax.get_legend().get_lines()
        print(lines_legend_to_adjust)
        count=1

        for line in lines_legend_to_adjust:
            print(len(lines_to_adjust))
            origColor = line.get_color()
            line.set_color('black')
            line.set_dashes(COLORMAP[origColor]['dash'])
            print(count/7)
            line.set_marker(MARKER[count/7])
            count = count +1
            line.set_markersize(MARKERSIZE)




    except AttributeError:
        pass


def setFigLinesBW(fig):
    """
    Take each axes in the figure, and for each line in the axes, make the
    line viewable in black and white.
    """
    for ax in fig.get_axes():
        setAxLinesBW(ax)
