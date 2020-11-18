from matplotlib import pyplot as plt
import pandas as pd
import numpy as np

# helper function to get average stats table
def average_stats(stats_table):
    avg_stat = pd.DataFrame()
    for stat,_ in stats_table:
        avg_stat = avg_stat.add(stat,fill_value=0)
    # plot average result
    return avg_stat/len(stats_table)


# plot cases vs time
def plot_cases(ax, stats_table, title, save_name = None):
    for stat,_ in stats_table:
        ax.plot(stat['E'],color ='gray', alpha = 0.1)
        ax.plot(stat['symptomatic_cases'],color ='gray', alpha = 0.1)
        ax.plot(stat['accumulated_cases'].diff(),color ='gray', alpha = 0.1)
    # plot average result
    avg_stat = average_stats(stats_table)
    ax.plot(avg_stat['E'],'--',label = 'Asymptomatic Cases')
    ax.plot(avg_stat['symptomatic_cases'],label = 'Symptomatic Cases')
    ax.plot(avg_stat['accumulated_cases'].diff(),'--',label = 'Daily New Cases') 
    ax.set_xlabel('Days')
    ax.set_ylabel('Fraction of Population')
    ax.set_xlim(0,len(avg_stat)-1)
    ax.set_ylim(bottom=0)
    ax.legend()
    ax.set_title(title)
    return avg_stat

# Plot the population being quaratined and isolated
def plot_quarantine(ax, stats_table,title, save_name = None):
    avg_stat = average_stats(stats_table)
    ax.stackplot(avg_stat.index, avg_stat['J'],avg_stat['Q'],
                 labels = ['Isolated','Quarantined'],alpha=0.5)
    ax.set_xlabel('Days')
    ax.set_ylabel('Fraction of Population')
    ax.set_xlim(0,len(avg_stat)-1)
    ax.set_ylim(bottom=0)
    ax.legend()
    ax.set_title(title)
    return avg_stat
            
def plot_rs(ax, stats_table,title):
    axr = ax.twinx()
    for stat, _ in stats_table:
        ax.plot(stat['R'],color ='gray', alpha = 0.1)
        axr.plot(stat['accumulated_cases'],color = 'gray',alpha = 0.1)
    # plot average result
    avg_stat = average_stats(stats_table)
    ax.plot(avg_stat['R'],label = 'Reproduction Number')
    ax.annotate(f' R0: {avg_stat["R"].max():.2f}', 
                (np.argmax(avg_stat["R"]),avg_stat["R"].max()))
    axr.plot(avg_stat['accumulated_cases'],'r--',label = 'accumulated Cases')
    ax.set_xlabel('Days')
    ax.set_xlim(0,len(avg_stat)-1)
    ax.set_ylabel('Reproduction Number',color='tab:blue')
    axr.set_ylabel('Accumulated Cases (Fraction)',color='tab:red')
    ax.set_ylim(bottom=0)
    axr.set_ylim(bottom=0)
    ax.set_title(title)
    return avg_stat
    
def plot_compare(ax, stat_tables, column, labels = None, save_name = None):
    if labels is None:
        labels = [None for _ in stat_tables]
    # plot lines
    for i, stat in enumerate(stat_tables):
        ax.plot(stat[column],label = labels[i])
        ax.set_xlabel('Days')
        ax.set_ylabel(column)
        ax.set_xlim([0,len(stat)])
        ax.set_ylim(bottom=0)
        ax.legend()

def plot_max_compare(ax, vals, stat_tables, column, title = None,
                     x_name = None, y_name = None, save_name = None):
    ax.plot(vals, [res[column].max() for res in stat_tables],'.-')
    ax.set_title(title)
    ax.set_xlabel(x_name)
    ax.set_ylabel(y_name)
    ax.set_xlim([min(vals),max(vals)])
    ax.set_ylim(bottom = 0)
