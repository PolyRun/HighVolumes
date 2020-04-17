import matplotlib.pyplot as plt

def plot(path, plot_name, dict_list):
    print("Start plotting")
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
    
    
    time_y_pos = list(range(1, len(time_labels)+1))
    convergence_y_pos = list(range(1, len(convergence_labels)+1))

    plt.title('Runtime comparison (mean)')
    plt.xlabel('Functions')
    plt.ylabel('Mean time')
    
    plt.xticks(time_y_pos, time_labels, rotation=90)
	
    #plt.bar(time_y_pos, heights)
    plt.errorbar(time_y_pos, time_heights, time_std_devs, ecolor='red', linestyle='None', marker='^')
	
    plt.ylim(bottom=0)

    plt.savefig(path+"/plots/"+plot_name+"runtime_mean.png", bbox_inches = "tight")

    plt.clf()
	
    # Plot of convergence if available
    if convergence:
        plt.title('Convergence comparison (mean)')
        plt.xlabel('Functions')
        plt.ylabel('Mean convergence')
    
        plt.xticks(convergence_y_pos, convergence_labels, rotation=90)
	
        #plt.bar(convergence_y_pos, convergence_heights)
        plt.errorbar(convergence_y_pos, convergence_heights, convergence_std_devs, ecolor='red', linestyle='None', marker='^')
	
        plt.ylim(bottom=0)

        plt.savefig(path+"/plots/"+plot_name+"convergence_mean.png", bbox_inches = "tight")

    print ("Finished plotting")