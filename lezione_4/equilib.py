import numpy as np
import torch
from torch import nn
from tqdm import tqdm
import os
import pandas as pd

# Neural Network definition
class NeuralNetwork(nn.Module):
    def __init__(self):
        super(NeuralNetwork, self).__init__()
        self.linear_elu_stack = nn.Sequential(
            nn.Linear(3,50,True),
            nn.ReLU(),
            nn.Linear(50,50,True),
            nn.ReLU(),
            nn.Linear(50,20,True),
            nn.ReLU(),
            nn.Linear(20,1,True)
        )

    def forward(self, x):
        output = self.linear_elu_stack(x)
        return output

#function to find the initial temperature
def Find_initial_temperature(r_cut,rho,fin_temp):
    # Get cpu, gpu device for training.
    device = (
        "cuda"
        if torch.cuda.is_available()
        else "cpu"
    )
    print(f"Using {device} device")

    norms=pd.read_csv("./IA_equilib/norm.csv",index_col=0)

    data= torch.tensor([float((r_cut-norms['Min']['r_cut'])/(norms['Max']['r_cut']-norms['Min']['r_cut'])),float((rho-norms['Min']['rho'])/(norms['Max']['rho']-norms['Min']['rho'])),float((fin_temp-norms['Min']['temp_fin'])/(norms['Max']['temp_fin']-norms['Min']['temp_fin']))], device=device)

    model=NeuralNetwork().to(device)
    model.load_state_dict(torch.load("./IA_equilib/model_val.pt", map_location=torch.device(device)))
    model.eval()

    with torch.no_grad():
        y=model(data)

    return round(float(y.cpu().numpy())*(norms['Max']['temp_init']-norms['Min']['temp_init'])+norms['Min']['temp_init'],4)





