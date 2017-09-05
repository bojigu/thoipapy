
def create_dict_organising_subplots(n_plots_per_fig,n_rows,n_cols):
    '''
    Function to help organise the creation of figures that contain multiple plots.
    For example, 15 histograms printed in figures with 8 histograms per figure/page.
    Returns a dict that gives a tuple for each plot/graph.
    newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr
    row_nr and col_nr are used to index pyplot subplots as follows
    fig, axarr = plt.subplots(2,2)
    _ = axarr[row_nr,col_nr].plot(x, y)
    '''
    dict_organising_subplots = {}
    #figure number
    fig_nr = 0
    #plot number in figure
    plot_nr_in_fig = 0
    #row number in figure
    row_nr = 0
    #column number in figure
    col_nr = 0
    #whether the figure needs to be saved
    savefig = False
    #whether a new figure needs to be created
    newfig = True

    for plotnr in range(1, 500):
        #add current counters to dict
        dict_organising_subplots[plotnr] = (newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr)
        plot_nr_in_fig += 1
        row_nr += 1
        newfig = False
        savefig = False
        #if plot_nr_in_fig is the last one before the new figure, then savefig = True
        if plot_nr_in_fig % (n_plots_per_fig - 1) == 0 and plot_nr_in_fig != 0:
            savefig = True
        #if plot_nr_in_fig is in a multiple of n_rows, then the plot goes to the second column
        if plot_nr_in_fig % n_rows == 0 and plot_nr_in_fig != 0:
            col_nr += 1
            row_nr = 0
        #if the plotnr is in a multple of n_plots_per_fig, then a new figure needs to created, and everything else reset
        if plotnr % n_plots_per_fig == 0 and plotnr != 0:
            #go to second figure
            fig_nr += 1
            #reset values
            plot_nr_in_fig = 0
            row_nr = 0
            col_nr = 0
            newfig = True
    return dict_organising_subplots