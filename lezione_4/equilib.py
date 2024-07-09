import numpy as np
import torch
from torch import nn
from torch.utils.data import TensorDataset, DataLoader
from tqdm import tqdm
import os
from equilib_nn import NeuralNetwork
import pandas as pd

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

    return y.cpu().numpy()


print(Find_initial_temperature(2.0,1.2,0.3))



