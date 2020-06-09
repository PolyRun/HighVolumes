import matplotlib.pyplot as plt
import numpy as np
import re
import pprint
import os


SAVEEPS = True
SAVEPNG = True

if 'ERRBAR' in os.environ and os.environ['ERRBAR'] == 'On': 
    PLOT_ERRORBARS = True
else:
    PLOT_ERRORBARS = False

PEAK_PERFORMANCE = 16
MEMORY_BANDWIDTH = 96
STREAM_BANDWIDTH = 12

#		 [Runt , Perf , I/O  , Roofl
ADD_PERF_ROOFS = [False, True, False, True]
ADD_MEM_ROOFS  = [False, False, False, True]

# Determines if PEAK_PERFORMANCE and MEMORY_BANDWIDTH are included into plots
# ADD_X_ROOFS has to be set to true to show any roof
MACHINE_ROOFS  = [False, False, False, False]

ROOFLINE_LOG = False


def finalize_and_save(plt, name):
    
    x1, x2 = plt.xlim()
    y1, y2 = plt.ylim()
    if "XLOG" in os.environ and os.environ["XLOG"] == 'On':
        plt.xscale('log')
        #assert(x1 > 0 and "cannot set aspect ratio of log axis with left border <= 0")
        print(str(x1), str(x2))
        x1, x2 = map(np.log10, plt.xlim())

    plt.ylim(bottom=0)
	
    plt.legend(bbox_to_anchor=(1, 1), loc="upper left")

    # show at most 12 xlabels
    ax = plt.gca()
    labels = ax.get_xaxis().get_ticklabels()
    keepith = (len(labels) + 11) // 12
    for i,label in enumerate(ax.get_xaxis().get_ticklabels()):
        if i % keepith != 0:
            label.set_visible(False)
            
    #set aspect ratio to 2/3
    ratio = (3.0*(x2 - x1))/(5.0*(y2 - y1))
    #print(str(x2-x1), str(y2-y1), str(ratio))
    ax.set_aspect(ratio)

        
    if SAVEEPS:
        plt.savefig(name+".eps", bbox_inches = "tight", format = 'eps', dpi = 1200)
    if SAVEPNG:
        plt.savefig(name+".png", bbox_inches = "tight")

    plt.clf()


def plot(path, plot_name, dict_list, x_option, title, x_label, y_label, perf_roofs, mem_roofs):
    if plot_name == None:
        plot_name = "plot"
    plot_name = plot_name + "_"
    #perf_roofs += [PEAK_PERFORMANCE,]
    mem_roofs += [STREAM_BANDWIDTH, ] #MEMORY_BANDWIDTH
   
    x_values = []
    x_flops = {}
    x_bytes = {}
    time_function_names = []
    time_function_heights = {}
    time_function_std_devs = {}
    time_function_ci_high = {}
    time_function_ci_low = {}
    for dict in dict_list:
        # cut out x_label option
        pattern = '{}=[^,]*'.format(x_option)
        name = ''.join(re.compile(pattern).split(dict[0]))
        time_dict = (dict[1])['time']
        ival = dict[2]
        time_fun_name = time_dict['name_t']+name
        time_fun_height = float(time_dict['mean'])
        time_fun_std_dev = float(time_dict['std_dev'])
        time_fun_ci_high = float(time_dict['ci_high'])  
        time_fun_ci_low = float(time_dict['ci_low'])   
        performance_counters = (dict[1])['performance_counter']
        flops = performance_counters['flops']
        bytes = performance_counters['bytes']
        if not ival in x_values:
            x_values.append(ival)
        if time_fun_name not in time_function_names:
            time_function_names.append(time_fun_name)
            time_function_heights[time_fun_name] = [time_fun_height]
            time_function_std_devs[time_fun_name] = [time_fun_std_dev]
            time_function_ci_high[time_fun_name] = [time_fun_ci_high-time_fun_height]
            time_function_ci_low[time_fun_name] = [time_fun_height-time_fun_ci_low]
            x_flops[time_fun_name] = [flops]
            x_bytes[time_fun_name] = [bytes]
        else:
            time_function_heights[time_fun_name].append(time_fun_height)
            time_function_std_devs[time_fun_name].append(time_fun_std_dev)
            time_function_ci_high[time_fun_name].append(time_fun_ci_high-time_fun_height)
            time_function_ci_low[time_fun_name].append(time_fun_height-time_fun_ci_low)
            x_flops[time_fun_name].append(flops)
            x_bytes[time_fun_name].append(bytes)

    # Runtime Plot
    plt.title(title[0] + '\n' + y_label[0] + " vs. " + x_label[0], loc='left', fontsize=12, pad=6)
	        
    x_ticks = [int(val) for val in x_values]
    plt.xticks(x_ticks, fontsize=12)
    plt.yticks(None, fontsize=12)

    for i, name in enumerate(time_function_names):
        #plt.plot(x_ticks, time_function_heights[name], label=name)
        pprint.pprint(x_ticks)
        pprint.pprint(time_function_heights[name])
        x_ticks_tmp = [tick for tick,val in zip(x_ticks, time_function_heights[name])]
        pprint.pprint(x_ticks_tmp)
        assert(len(x_ticks_tmp) == len(time_function_heights[name]) and "maybe you used 'generator' in some parameter? or some configs identical?")
        plt.errorbar(
            x_ticks_tmp,
            time_function_heights[name],
            label=name,
            yerr= [time_function_ci_low[name], time_function_ci_high[name]] if PLOT_ERRORBARS else [
                [0 for i in time_function_ci_low[name]],[0 for i in time_function_ci_high[name]]],
            capsize=4
        )

    finalize_and_save(plt, path+"/plots/"+plot_name+"runtime_mean")


    # Performance plot
    plt.title(title[1] + '\n' + y_label[1] + " vs. " + x_label[1], loc='left', fontsize=12, pad=6)

    plt.xticks(x_ticks, fontsize=12)
    plt.yticks(None, fontsize=12)

    for i, name in enumerate(time_function_names):
        x_ticks_tmp = [tick for tick,val in zip(x_ticks, time_function_heights[name])]
        time_function_performance = []
        time_function_performance_ci_low = []
        time_function_performance_ci_high = []
        for index, item in enumerate(time_function_heights[name]):
            time_function_performance.append(x_flops[name][index]/item)
            
            time_function_performance_ci_low.append(
                x_flops[name][index]/item - x_flops[name][index]/(time_function_heights[name][index] - time_function_ci_low[name][index])
                if PLOT_ERRORBARS
                else 0
            )
            time_function_performance_ci_high.append(
                x_flops[name][index] / (time_function_heights[name][index] + time_function_ci_high[name][index]) - x_flops[name][index]/item
                if PLOT_ERRORBARS
                else 0
            )
        #plt.plot(x_ticks, time_function_performance, label=name)
        plt.errorbar(x_ticks_tmp, time_function_performance, label=name, yerr=[time_function_performance_ci_low, time_function_performance_ci_high], capsize=4)

    if ADD_PERF_ROOFS[1]:
        if MACHINE_ROOFS[1]:
            perf_roofs.append(PEAK_PERFORMANCE)
        for r in perf_roofs:
            plt.axhline(r, linestyle='--', color='#808080')
        if MACHINE_ROOFS[1]:
            perf_roofs.pop()

    finalize_and_save(plt, path+"/plots/"+plot_name+"performance_mean")


    # I/O plot
    plt.title(title[2] + '\n' + y_label[2] + " vs. " + x_label[2], loc='left', fontsize=12, pad=6)

    plt.xticks(x_ticks, fontsize=12)
    plt.yticks(None, fontsize=12)

    for i, name in enumerate(time_function_names):
        x_ticks_tmp = [tick for tick,val in zip(x_ticks, time_function_heights[name])]
        time_function_bytes = []
        time_function_bytes_ci_low = []
        time_function_bytes_ci_high = []
        for index, item in enumerate(time_function_heights[name]):
            time_function_bytes.append(x_bytes[name][index]/item)            
            time_function_bytes_ci_low.append(
                x_bytes[name][index]/item - x_bytes[name][index]/(time_function_heights[name][index] - time_function_ci_low[name][index])
                if PLOT_ERRORBARS
                else 0
            )
            time_function_bytes_ci_high.append(
                x_bytes[name][index]/(time_function_heights[name][index] + time_function_ci_high[name][index]) - x_bytes[name][index]/item
                if PLOT_ERRORBARS
                else 0
            )
        #plt.plot(x_ticks, time_function_bytes, label=name)
        plt.errorbar(x_ticks_tmp, time_function_bytes, label=name, yerr=[time_function_bytes_ci_low, time_function_bytes_ci_high], capsize=4)

    if ADD_MEM_ROOFS[2]:
        if MACHINE_ROOFS[2]:
            mem_roofs.append(MEMORY_BANDWIDTH)
        for r in mem_roofs:
            plt.axhline(r, linestyle='--', color='#808080')
        if MACHINE_ROOFS[2]:
            mem_roofs.pop()

    finalize_and_save(plt, path+"/plots/"+plot_name+"io_mean")

    
    # Roofline plot
    plt.title("Roofline measurements \nPerformance [Flops/Cycle] vs. Operational Intensity [Flops/Byte]", loc='left', fontsize=12, pad=6)
	
    max_intensity = 0
    min_intensity = 1024
    min_performance = PEAK_PERFORMANCE+1

    i = 0
    f = None
    if ROOFLINE_LOG:
        f = open(path+"/plots/"+plot_name+"roofline.out", "w")

    for name in time_function_names:
        time_function_performance = []
        time_function_intensity = []
        if x_flops[name][0] == 0 or x_bytes[name][0] == 0:
            continue
        for index, item in enumerate(time_function_heights[name]):
            time_function_performance.append(x_flops[name][index]/item)
            time_function_intensity.append(x_flops[name][index]/x_bytes[name][index])
            max_intensity = max(max_intensity, x_flops[name][index]/x_bytes[name][index])
            min_intensity = min(min_intensity, x_flops[name][index]/x_bytes[name][index])
            min_performance = min(min_performance, x_flops[name][index]/item)
        if ROOFLINE_LOG:
            print("{}:".format(name))
            print("Performance: ", time_function_performance)
            print("O-Intensity: ", time_function_intensity)
            f.write("{}:\n".format(name))
            f.write("Performance: "+ str(time_function_performance) + "\n")
            f.write("O-Intensity: "+ str(time_function_intensity) + "\n")
        plt.plot(time_function_intensity, time_function_performance, label=name, marker=".")
        i += 1

    if ROOFLINE_LOG:
        f.close()

    right_add = 0.5*max_intensity
    left_sub = 0.5*min_intensity
    bottom_sub = 0.5*min_performance
    right_bound= max_intensity+right_add
    left_bound = min_intensity-left_sub
    bottom_bound = min_performance-bottom_sub

    if ADD_PERF_ROOFS[3]:
        if MACHINE_ROOFS[3]:
            perf_roofs.append(PEAK_PERFORMANCE)
            mem_roofs.append(MEMORY_BANDWIDTH)
        for r in perf_roofs:
            plt.axhline(r, linestyle='--', color='#808080')
        for r in mem_roofs:
            ridge_point_intensity = PEAK_PERFORMANCE / r
            mem_bound_left_x = left_bound
            mem_bound_left_y = PEAK_PERFORMANCE*(mem_bound_left_x/ridge_point_intensity)
            mem_bound_right_x = max_intensity + right_add
            mem_bound_right_y = PEAK_PERFORMANCE*(mem_bound_right_x/ridge_point_intensity)
            plt.plot([mem_bound_left_x, mem_bound_right_x], [mem_bound_left_y, mem_bound_right_y], linestyle='--', color='#808080')
        if MACHINE_ROOFS[3]:
            perf_roofs.pop()
            mem_roofs.pop()


    plt.xscale('log')
    plt.yscale('log')
    plt.xticks(None, fontsize=12)
    plt.yticks(None, fontsize=12)

    
    plt.xlim((left_bound, right_bound))
    
    plt.ylim((bottom_bound, 1.5*max(perf_roofs + [plt.ylim()[1]])))

    
    #set aspect ratio to 2/3
    if (left_bound <= 0 or bottom_bound <= 0):
        print("cannot set aspect ratio because left_bound or bottom_bound <= 0")
    else:
        ax = plt.gca()
        x1, x2 = map(np.log10, plt.xlim())
        y1, y2 = map(np.log10, plt.ylim())
        ratio = (3.0*(x2 - x1))/(5.0*(y2 - y1))
        ax.set_aspect(ratio)
    
    plt.legend(bbox_to_anchor=(1, 1), loc="upper left")

    if SAVEEPS:
        plt.savefig(path+"/plots/"+plot_name+"roofline.eps", bbox_inches = "tight", format = 'eps', dpi=1200)
    if SAVEPNG:
        plt.savefig(path+"/plots/"+plot_name+"roofline.png", bbox_inches = "tight")

    plt.clf()


