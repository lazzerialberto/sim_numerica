import numpy as np
import torch
from torch import nn
from tqdm import tqdm
import os
import pandas as pd
import argparse

'''parser = argparse.ArgumentParser()
parser.add_argument("--r_cut", help="cutoff radius", type=float, required=True)
parser.add_argument("--rho", help="density", type=float, required=True)
parser.add_argument("--fin_temp", help="final temperature deisdered", type=float, required=True)
args = parser.parse_args()'''

# Neural Network definition
class NeuralNetwork(nn.Module):
    def __init__(self):
        super(NeuralNetwork, self).__init__()
        self.linear_elu_stack = nn.Sequential(
            nn.Linear(3,100,True),
            nn.ELU(),
            nn.Linear(100,100,True),
            nn.ELU(),
            nn.Linear(100,50,True),
            nn.ELU(),
            nn.Linear(50,20,True),
            nn.ELU(),
            nn.Linear(20,1,True)
        )

    def forward(self, x):
        output = self.linear_elu_stack(x)
        return output

# Get cpu, gpu device for training.
device = (
    "cuda"
    if torch.cuda.is_available()
    else "cpu"
)
print(f"Using {device} device")

in_file=pd.read_fwf("./NSL_SIMULATOR/INPUT/input.dat",names=["word", "value"])

data_dict = {key: in_file["value"][i].tolist() for i, key in enumerate(in_file['word'])}

norms=pd.read_csv("./IA_equilib/norm.csv",index_col=0)

data= torch.tensor([float((data_dict['R_CUT']-norms['Min']['r_cut'])/(norms['Max']['r_cut']-norms['Min']['r_cut'])),float((data_dict['RHO']-norms['Min']['rho'])/(norms['Max']['rho']-norms['Min']['rho'])),float((data_dict['TEMP']-norms['Min']['temp_fin'])/(norms['Max']['temp_fin']-norms['Min']['temp_fin']))], device=device)

model=NeuralNetwork().to(device)
model.load_state_dict(torch.load("./IA_equilib/model_val.pt", map_location=torch.device(device)))
model.eval()

with torch.no_grad():
    y=model(data)

print(round(float(y.cpu().numpy())*(norms['Max']['temp_init']-norms['Min']['temp_init'])+norms['Min']['temp_init'],4))

with open("./NSL_SIMULATOR/INPUT/input.dat", 'r') as file:
    lines = file.readlines()

insert_position = 9 
lines.insert(insert_position,   f"INIT_TEMP\t\t\t   "+str(round(float(y.cpu().numpy())*(norms['Max']['temp_init']-norms['Min']['temp_init'])+norms['Min']['temp_init'],4)))

insert_position = 10 
lines.insert(insert_position, "\n")

with open("./NSL_SIMULATOR/INPUT/input.dat", 'w') as file:
    file.writelines(lines)




