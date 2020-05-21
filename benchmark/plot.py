import matplotlib.pyplot as plt
import numpy as np
import re
import pprint

SAVEEPS = False
PEAK_PERFORMANCE = 20
MEMORY_BANDWIDTH = 32

def plot(path, plot_name, dict_list, x_option, title, x_label, y_label):
    if plot_name == None:
        plot_name = "plot"
    plot_name = plot_name + "_"
   
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
    plt.title(title[0])
    plt.xlabel(x_label[0])
    plt.ylabel(y_label[0])
	        
    x_ticks = []
    for val in x_values:
        x_ticks.append(int(val))
    plt.xticks(x_ticks)
	

    i = 0
    for name in time_function_names:
        #plt.plot(x_ticks, time_function_heights[name], label=name)
        pprint.pprint(x_ticks)
        pprint.pprint(time_function_heights[name])
        assert(len(x_ticks) == len(time_function_heights[name]) and "maybe you used 'generator' in some parameter? or some configs identical?")
        plt.errorbar(x_ticks, time_function_heights[name], label=name, yerr=[time_function_ci_low[name], time_function_ci_high[name]], capsize=4)
        i += 1
	
    plt.ylim(bottom=0)
	
    plt.legend(bbox_to_anchor=(1, 1), loc="upper left")

    if SAVEEPS:
        plt.savefig(path+"/plots/"+plot_name+"runtime_mean.eps", bbox_inches = "tight", format = 'eps', dpi = 1200)
    else:
        plt.savefig(path+"/plots/"+plot_name+"runtime_mean.png", bbox_inches = "tight")

    plt.clf()

    # Performance plot
    plt.title(title[1])
    plt.xlabel(x_label[1])
    plt.ylabel(y_label[1])

    plt.xticks(x_ticks)

    i = 0
    for name in time_function_names:
        time_function_performance = []
        time_function_performance_ci_low = []
        time_function_performance_ci_high = []
        for index, item in enumerate(time_function_heights[name]):
            time_function_performance.append(x_flops[name][index]/item)
            time_function_performance_ci_low.append(x_flops[name][index]/item - x_flops[name][index]/(time_function_heights[name][index] - time_function_ci_low[name][index]))
            time_function_performance_ci_high.append(x_flops[name][index] / (time_function_heights[name][index] + time_function_ci_high[name][index]) - x_flops[name][index]/item)
        #plt.plot(x_ticks, time_function_performance, label=name)
        plt.errorbar(x_ticks, time_function_performance, label=name, yerr=[time_function_performance_ci_low, time_function_performance_ci_high], capsize=4)
        i += 1
	
    plt.ylim(bottom=0)
	
    plt.legend(bbox_to_anchor=(1, 1), loc="upper left")

    if SAVEEPS:
        plt.savefig(path+"/plots/"+plot_name+"performance_mean.eps", bbox_inches = "tight", format = 'eps', dpi=1200)
    else:
        plt.savefig(path+"/plots/"+plot_name+"performance_mean.png", bbox_inches = "tight")

    plt.clf()

    # I/O plot
    plt.title(title[2])
    plt.xlabel(x_label[2])
    plt.ylabel(y_label[2])

    plt.xticks(x_ticks)

    i = 0
    for name in time_function_names:
        time_function_bytes = []
        time_function_bytes_ci_low = []
        time_function_bytes_ci_high = []
        for index, item in enumerate(time_function_heights[name]):
            time_function_bytes.append(x_bytes[name][index]/item)            
            time_function_bytes_ci_low.append(x_bytes[name][index]/item - x_bytes[name][index]/(time_function_heights[name][index] - time_function_ci_low[name][index]))
            time_function_bytes_ci_high.append(x_bytes[name][index]/(time_function_heights[name][index] + time_function_ci_high[name][index]) - x_bytes[name][index]/item)
        #plt.plot(x_ticks, time_function_bytes, label=name)
        plt.errorbar(x_ticks, time_function_bytes, label=name, yerr=[time_function_bytes_ci_low, time_function_bytes_ci_high], capsize=4)
        i += 1
	
    plt.ylim(bottom=0)
	
    plt.legend(bbox_to_anchor=(1, 1), loc="upper left")

    if SAVEEPS:
        plt.savefig(path+"/plots/"+plot_name+"io_mean.eps", bbox_inches = "tight", format = 'eps', dpi=1200)
    else:
        plt.savefig(path+"/plots/"+plot_name+"io_mean.png", bbox_inches = "tight")

    plt.clf()

    # Roofline plot
    plt.title(title[3])
    plt.xlabel(x_label[3])
    plt.ylabel(y_label[3])
    plt.title("Roofline measurements")
    plt.xlabel("Operational Intensity [Flops/Byte]")
    plt.ylabel("Performance [Flops/Cycle]")
	
    max_intensity = 0;

    i = 0
    for name in time_function_names:
        time_function_performance = []
        time_function_intensity = []
        for index, item in enumerate(time_function_heights[name]):
            time_function_performance.append(x_flops[name][index]/item)
            time_function_intensity.append(x_flops[name][index]/x_bytes[name][index])
            max_intensity = max(max_intensity, x_flops[name][index]/x_bytes[name][index])
        #print("Performance:", time_function_performance)
        #print("O-Intensity:", time_function_intensity)
        plt.plot(time_function_intensity, time_function_performance, label=name)
        i += 1

    LEFT_BOUND = 0.05
    LOWER_BOUND = 0.05
    right_add = 0.1

    ridge_point_intensity = PEAK_PERFORMANCE / MEMORY_BANDWIDTH

    mem_bound_left_x = LEFT_BOUND
    mem_bound_left_y = PEAK_PERFORMANCE*(mem_bound_left_x/ridge_point_intensity)
    mem_bound_right_x = max_intensity + right_add
    mem_bound_right_y = PEAK_PERFORMANCE*(mem_bound_right_x/ridge_point_intensity)
	
    # Draw Memory bound until max_intensity (+right_add)
    plt.plot([mem_bound_left_x, mem_bound_right_x], [mem_bound_left_y, mem_bound_right_y], linestyle='--', alpha=0.5, color='black')
    # Draw Memory bound until ridge point
    #plt.plot([mem_bound_left_x, ridge_point_intensity], [mem_bound_left_y, PEAK_PERFORMANCE], linestyle='--', alpha=0.5, color='black')
		
    plt.axhline(PEAK_PERFORMANCE, linestyle='--', alpha=0.5, color='black')

    plt.xscale('log')
    plt.yscale('log')
	
    plt.xlim((LEFT_BOUND, max_intensity + right_add))
    plt.ylim((LOWER_BOUND, 2*PEAK_PERFORMANCE))
	
    plt.legend(bbox_to_anchor=(1, 1), loc="upper left")

    if SAVEEPS:
        plt.savefig(path+"/plots/"+plot_name+"roofline.eps", bbox_inches = "tight", format = 'eps', dpi=1200)
    else:
        plt.savefig(path+"/plots/"+plot_name+"roofline.png", bbox_inches = "tight")

    plt.clf()


