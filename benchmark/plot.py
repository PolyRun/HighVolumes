import matplotlib.pyplot as plt
import numpy as np

def plot(path, plot_name, dict_list):
    if plot_name == None:
        plot_name = "plot"
    plot_name = plot_name + "_"
	
    convergence = False
   
    time_labels = []
    time_heights = []
    time_std_devs = []
    convergence_labels = []
    convergence_heights = []
    convergence_std_devs = []
    for dict in dict_list:
        name = dict[0]
        time_dict = (dict[1])['time']
        time_labels.append(time_dict['name_t']+name)
        time_heights.append(float(time_dict['mean']))
        time_std_devs.append(float(time_dict['std_dev']))
        convergence_dict = (dict[1])['convergence']
        if convergence_dict:
            convergence = True
            convergence_labels.append(convergence_dict['name_c']+name)
            convergence_heights.append(float(convergence_dict['mean']))
            convergence_std_devs.append(float(convergence_dict['std_dev']))
    
    
    time_y_pos = np.arange(len(time_labels))
    convergence_y_pos = np.arange(len(convergence_labels))

    plt.title('Runtime comparison (mean)')
    plt.xlabel('Functions')
    plt.ylabel('Mean time')
    
    plt.xticks(time_y_pos, time_labels, rotation=90)
    bar_width = 0.1
	
    plt.bar(time_y_pos, time_heights, yerr=time_std_devs, ecolor='red', width=bar_width)
	
    plt.ylim(bottom=0)

    plt.savefig(path+"/plots/"+plot_name+"runtime_mean.png", bbox_inches = "tight")

    plt.clf()
	
	
    # Plot of convergence if available
    if convergence:
        plt.title('Convergence comparison (mean)')
        plt.xlabel('Functions')
        plt.ylabel('Mean convergence')
    
        plt.xticks(convergence_y_pos, convergence_labels, rotation=90)
	
        plt.bar(convergence_y_pos, convergence_heights, yerr=convergence_std_devs, ecolor='red', width=bar_width)
	
        plt.ylim(bottom=0)

        plt.savefig(path+"/plots/"+plot_name+"convergence_mean.png", bbox_inches = "tight")
    plt.clf()



def plot_input(path, plot_name, dict_list, x_label):
    if plot_name == None:
        plot_name = "plot"
    plot_name = plot_name + "_"
	
    convergence = False
   
    x_values = []
    x_flops = {}
    x_bytes = {}
    time_function_names = []
    time_function_heights = {}
    time_function_std_devs = {}
    convergence_function_names = []
    convergence_function_heights = {}
    convergence_function_std_devs = {}
    for dict in dict_list:
        name = "_"+dict[0].split("=")[1]
        time_dict = (dict[1])['time']
        ival = dict[2]
        time_fun_name = time_dict['name_t']+name
        time_fun_height = float(time_dict['mean'])
        time_fun_std_dev = float(time_dict['std_dev'])        
        performance_counters = (dict[1])['performance_counter']
        flops = performance_counters['flops']
        bytes = performance_counters['bytes']
        if not ival in x_values:
            x_values.append(ival)
        if time_fun_name not in time_function_names:
            time_function_names.append(time_fun_name)
            time_function_heights[time_fun_name] = [time_fun_height]
            time_function_std_devs[time_fun_name] = [time_fun_std_dev]
            x_flops[time_fun_name] = [flops]
        else:
            time_function_heights[time_fun_name].append(time_fun_height)
            time_function_std_devs[time_fun_name].append(time_fun_std_dev)
            x_flops[time_fun_name].append(flops)

            
        convergence_dict = (dict[1])['convergence']
        if convergence_dict:
            convergence = True
            convergence_fun_name = convergence_dict['name_c']+name
            convergence_fun_height = float(convergence_dict['mean'])
            convergence_fun_std_dev = float(convergence_dict['std_dev'])
            if convergence_fun_name not in convergence_function_names:
                convergence_function_names.append(convergence_fun_name)
                convergence_function_heights[convergence_fun_name] = [convergence_fun_height]
                convergence_function_std_devs[convergence_fun_name] = [convergence_fun_std_dev]
            else:
                convergence_function_heights[convergence_fun_name].append(convergence_fun_height)
                convergence_function_std_devs[convergence_fun_name].append(convergence_fun_std_dev)
    

    # Runtime Plot
    plt.title('Runtime comparison (mean)')
    plt.xlabel(x_label)
    plt.ylabel('Mean time')
	
    bar_width = 0.1
    
    y_pos = [ [] for j in range(len(time_function_names)) ]
    y_pos[0] = np.arange(len(x_values))

    for i in range(1, len(time_function_names)):
        for j in range(len(y_pos[i-1])):            
            y_pos[i].append(y_pos[0][j]+i*bar_width)
        i += 1
        
    x_ticks = [r + bar_width for r in range(len(x_values))]
    for index, item in enumerate(x_ticks):
        if len(time_function_names) % 2 == 1:
            x_ticks[index] += bar_width
        else:
            x_ticks[index] += 0.5*bar_width
    plt.xticks(x_ticks, x_values)
	

    i = 0
    for name in time_function_names:
        plt.bar(y_pos[i], time_function_heights[name], yerr=time_function_std_devs[name], ecolor='red', width=bar_width, label=name)
        i += 1
	
    plt.ylim(bottom=0)
	
    plt.legend(bbox_to_anchor=(1, 1), loc="upper left")

    plt.savefig(path+"/plots/"+plot_name+"runtime_mean.png", bbox_inches = "tight")

    plt.clf()

    # Performance plot
    plt.title('Performance comparison (mean)')
    plt.xlabel(x_label)
    plt.ylabel('Flops/cycles(mean)')
	
    x_ticks = []
    for val in x_values:
        x_ticks.append(int(val))
	
    plt.xticks(x_ticks)	

    i = 0
    for name in time_function_names:
        time_function_performance = []
        for index, item in enumerate(time_function_heights[name]):
            time_function_performance.append(x_flops[name][index]/item)
        plt.plot(x_ticks, time_function_performance, label=name)
        i += 1
	
    plt.ylim(bottom=0)
	
    plt.legend(bbox_to_anchor=(1, 1), loc="upper left")

    plt.savefig(path+"/plots/"+plot_name+"performance_mean.png", bbox_inches = "tight")

    plt.clf()
	
    # Plot of convergence if available
    if convergence:
        plt.title('Convergence comparison (mean)')
        plt.xlabel(x_label)
        plt.ylabel('Mean time')
	
        bar_width = 0.1
    
        y_pos = [ [] for j in range(len(convergence_function_names)) ]
        y_pos[0] = np.arange(len(x_values))

        for i in range(1, len(convergence_function_names)):
            for j in range(len(y_pos[i-1])):            
                y_pos[i].append(y_pos[0][j]+i*bar_width)
            i += 1
        
        x_ticks = [r + bar_width for r in range(len(x_values))]
        for index, item in enumerate(x_ticks):
            if len(convergence_function_names) % 2 == 1:
                x_ticks[index] += bar_width
            else:
                x_ticks[index] += 0.5*bar_width
        plt.xticks(x_ticks, x_values)
	

        i = 0
        for name in convergence_function_names:
            plt.bar(y_pos[i], convergence_function_heights[name], yerr=convergence_function_std_devs[name], ecolor='red', width=bar_width, label=name)
            i += 1
	
        plt.ylim(bottom=0)

        plt.savefig(path+"/plots/"+plot_name+"convergence_mean.png", bbox_inches = "tight")
		
        plt.clf()

	