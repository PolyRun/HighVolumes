import matplotlib.pyplot as plt
import numpy as np
import math

def plot(path, dict_list):
    print("Start plotting")
   
    y_pos = list(range(1, len(dict_list)+1))
    heights = []
    labels = []
    for dict in dict_list:
        name = dict[0]
        time_dict = (dict[1])['time']
        convergence_dict = (dict[1])['convergence']
        labels.append(time_dict['name_t']+name)
        heights.append(float(time_dict['mean']))
    
    


    plt.title('Runtime comparison (mean)')
    plt.xlabel('Functions')
    plt.ylabel('Mean time')
    
    plt.xticks(y_pos, labels, rotation=90)
	
    plt.bar(y_pos, heights)
	
    plt.ylim(bottom=0)

    plt.savefig(path+"/plots/runtime_mean.png", bbox_inches = "tight")

    print ("Finished plotting")