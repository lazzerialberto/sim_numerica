import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
import os


data=np.loadtxt("./NSL_SIMULATOR/OUTPUT/parameters.dat",skiprows=1)

in_file=pd.read_fwf("./NSL_SIMULATOR/INPUT/input.dat",names=["word", "value"])
data_dict = {key: in_file["value"][i].tolist() for i, key in enumerate(in_file['word'])}

with open("./NSL_SIMULATOR/INPUT/input.dat", 'r') as file:
    lines = file.readlines()


model=LinearRegression()

if data_dict['R_CUT']==2.2 and data_dict['RHO']==1.1:
    model.fit(data[14:21,-1].reshape(-1,1),data[14:21,-2])
    y_pred=model.predict(data[14:21,-1].reshape(-1,1))
    x_pred=1./model.coef_[0]*(0.8-model.intercept_)
    insert_position = 9 
    lines.insert(insert_position,   f"INIT_TEMP\t\t\t   "+str(round(float(x_pred),4)))
    #plt.scatter(data[14:21,-2], data[14:21,-1], color='navy', label='Real data')
    plt.errorbar(data[14:21,-1], data[14:21,-2],yerr=np.sqrt(sum(((y_pred-data[14:21,-2])**2)/len(y_pred))), fmt='o',elinewidth=1.75, capsize=4, capthick=1.75, color='navy',ecolor='blue', label='Real data')
    plt.plot(data[14:21,-1], y_pred, color='red', label='Linear Fit')
    plt.xlabel('temp_init')
    plt.ylabel('temp_fin')
    plt.legend()
    plt.show()

elif data_dict['R_CUT']==2.5 and data_dict['RHO']==0.8:
    model.fit(data[7:14,-1].reshape(-1,1),data[7:14,-2])
    y_pred=model.predict(data[7:14,-1].reshape(-1,1))
    x_pred=1./model.coef_[0]*(1.1-model.intercept_)
    insert_position = 9 
    lines.insert(insert_position,   f"INIT_TEMP\t\t\t   "+str(round(float(x_pred),4)))
    #plt.scatter(data[7:14,-2], data[7:14,-1], color='navy', label='Real data')
    plt.errorbar(x=data[7:14,-1], y=data[7:14,-2],yerr=np.sqrt(sum(((y_pred-data[7:14,-2])**2)/len(y_pred))), fmt='o',elinewidth=1.75, capsize=4, capthick=1.75, color='navy',ecolor='blue', label='Real data')
    plt.plot(data[7:14,-1], y_pred, color='red', label='Linear Fit')
    plt.xlabel('temp_init')
    plt.ylabel('temp_fin')
    plt.legend()
    plt.show()

elif data_dict['R_CUT']==5.0 and data_dict['RHO']==0.05:
    model.fit(data[0:7,-1].reshape(-1,1),data[0:7,-2])
    y_pred=model.predict(data[0:7,-1].reshape(-1,1))
    x_pred=1./model.coef_[0]*(1.2-model.intercept_)
    insert_position = 9 
    lines.insert(insert_position,   f"INIT_TEMP\t\t\t   "+str(round(float(x_pred),4)))
    #plt.scatter(data[0:7,-2], data[0:7,-1], color='navy', label='Real data')
    plt.errorbar(data[0:7,-1], data[0:7,-2],yerr=np.sqrt(sum(((y_pred-data[0:7,-2])**2)/len(y_pred))), fmt='o',elinewidth=1.75, capsize=4, capthick=1.75, color='navy',ecolor='blue', label='Real data')
    plt.plot(data[0:7,-1], y_pred, color='red', label='Linear Fit')
    plt.xlabel('temp_init')
    plt.ylabel('temp_fin')
    plt.legend()
    plt.show()

    '''elif data_dict['R_CUT']==5.0 and data_dict['RHO']==0.05:
    model.fit(data[0:7,-2].reshape(-1,1),data[0:7,-1])
    y_pred=model.predict(data[0:7,-2].reshape(-1,1))
    insert_position = 9 
    lines.insert(insert_position,   f"INIT_TEMP\t\t\t   "+str(round(float(y_pred[4]),4)))
    #plt.scatter(data[0:7,-2], data[0:7,-1], color='navy', label='Real data')
    plt.errorbar(data[0:7,-2], data[0:7,-1],yerr=np.sqrt(sum(((y_pred-data[0:7,-1])**2)/len(y_pred))), fmt='o',elinewidth=1.75, capsize=4, capthick=1.75, color='navy',ecolor='blue', label='Real data')
    plt.plot(data[0:7,-2], y_pred, color='red', label='Linear Fit')
    plt.xlabel('temp_init')
    plt.ylabel('temp_fin')
    plt.legend()
    plt.show()'''



else:
    print("Error: cutoff radius doesn't match density")

insert_position = 10 
lines.insert(insert_position, "\n")

with open("./NSL_SIMULATOR/INPUT/input.dat", 'w') as file:
    file.writelines(lines)