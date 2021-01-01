import datetime as d
import os
def today():
	return d.datetime.today().strftime('%d.%m.%Y_%H:%M:%S')

def make_file_name(n, num_of_meas):
	return f"{today()}_{n}x{n}_{num_of_meas}"

"""def save_data(file_name, magn_list, magn_sig_list, sus_list, sus_sig_list):
	try:
	    with open(f'Results/data_' + file_name + '.csv', 'a') as file:
	        for ind in range(len(magn_list)):
	            file.write(f"beta = {betas[ind]}, <M> = {magn_list[ind]}, s<M> = {magn_sig_list[ind]}, chi = {sus_list[ind]}, s_chi = {sus_sig_list[ind]}\n")
	        print("wrote to file!")
	except:
	    print("write to file failed, what the hell?")"""


def save_data(file_name, betas, magn_list, magn_sig_list, sus_list, sus_sig_list):
	with open(f'Results/data_' + file_name + '.csv', 'a') as file:
		file.write("beta,M,sM,chi,schi\n")
		for ind in range(len(magn_list)):
			file.write(f"{betas[ind]},{magn_list[ind]},{magn_sig_list[ind]},{sus_list[ind]},{sus_sig_list[ind]}\n")
	print("wrote to file!")

def open_saved(n, num_of_meas):
	data_dict = {'beta':[], 'M':[], 'sM':[],'chi':[],'schi':[]}
	keys = list(data_dict.keys())
	#print(keys)
	for filename in os.listdir('Results'):
		if "data" in filename and f"{n}x{n}" in filename and str(num_of_meas) in filename:
			with open('Results/'+filename,'r') as file:
				lines = file.readlines()
				for line in lines[1:]:
					for index,val in enumerate(line.split(',')):
						#print(data_dict)
						#print(keys[index])
						data_dict[keys[index]] += [float(val)]
	return data_dict

if __name__=="__main__":
	print(open_saved(20,4000)['beta'])

