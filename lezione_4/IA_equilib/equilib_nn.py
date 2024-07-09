import numpy as np
import matplotlib.pyplot as plt
import torch
from torch import nn
from torch.utils.data import TensorDataset, DataLoader
from tqdm import tqdm
import os
import csv


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


data=np.loadtxt("../NSL_SIMULATOR/OUTPUT/parameters.dat",skiprows=1)

keys = ["r_cut", "rho", "temp_fin", "temp_init"]
data_dict = {key: data[:, i].tolist() for i, key in enumerate(keys)}

def normalize(data):
    min_val = min(data)
    max_val = max(data)
    if max_val!=min_val:
        return [round((x - min_val) / (max_val - min_val),6) for x in data]
    else:
        return [round(x,6) for x in data]

normalized_dict = {key: normalize(values) for key, values in data_dict.items()}

min_max_dict = {key: {'min': None, 'max': None} for key in keys}

for key in keys:
    min_max_dict[key]['min'] = min(data_dict[key])
    min_max_dict[key]['max'] = max(data_dict[key])

print(min_max_dict)

with open('./norm.csv', 'w', newline='') as csv_file:
    writer = csv.writer(csv_file)
    # Scrivere l'intestazione
    writer.writerow(['Key', 'Min', 'Max'])
    # Scrivere i valori
    for key, value in min_max_dict.items():
        writer.writerow([key, value['min'], value['max']])


r_cut=torch.tensor(normalized_dict["r_cut"])
rho=torch.tensor(normalized_dict["rho"])
temp_fin=torch.tensor(normalized_dict["temp_fin"])
temp_init =torch.tensor(normalized_dict["temp_init"])

n_samples=len(temp_fin)
len_train_set=int(n_samples*0.8)

train_load = TensorDataset(torch.column_stack((r_cut[0:len_train_set-1],rho[0:len_train_set-1],temp_fin[0:len_train_set-1])),temp_init[0:len_train_set-1]) # create your dataset
train_set = DataLoader(train_load, batch_size=1)

val_load = TensorDataset(torch.column_stack((r_cut[len_train_set-1:n_samples],rho[len_train_set-1:n_samples],temp_fin[len_train_set-1:n_samples])),temp_init[len_train_set-1:n_samples]) # create your dataset
val_set = DataLoader(val_load, batch_size=1)
    
model = NeuralNetwork().to(device)
print(model)

optimizer = torch.optim.Adam(model.parameters(), lr=0.0001)
scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, 'min', factor=0.8, patience=1, verbose=True)
n_epochs=100
loss_train_epochs=[]
loss_val_epochs=[]

print("\nStart training")

min_train_loss=1000
min_val_loss=1000

for epoch in tqdm(range(n_epochs)):

    cumulative_loss=0
    size = len(train_set.dataset)
    model.train()

    for batch, (X , Y) in enumerate(train_set):
        X, Y = X.to(device), Y.to(device)

        pred = model(X)
        loss = 0
        for i in range(len(Y)):
            loss+=(Y[i]-pred[i])**2

        loss /= len(Y)
        cumulative_loss+=loss.detach().cpu().numpy()
        loss.backward()
        optimizer.step()
        optimizer.zero_grad()

    loss_train_epochs.append(cumulative_loss/len(train_set))

    if epoch==10:
        min_train_loss=cumulative_loss
        torch.save(model.state_dict(),"./model_train.pt")
    if min_train_loss>cumulative_loss:
        min_train_loss=cumulative_loss
        torch.save(model.state_dict(),"./model_train.pt")


    cumulative_val_loss=0
    num_batches = len(val_set)
    model.eval()
    test_loss = 0
    with torch.no_grad():
        for X, Y in val_set:
            X, Y = X.to(device), Y.to(device)
            pred = model(X)
            loss = 0
            for i in range(len(Y)):
                loss+=(Y[i]-pred[i])**2

            loss /= len(Y)
            test_loss+=loss
    test_loss /= num_batches
    loss_val_epochs.append(test_loss.cpu().numpy())

    if epoch==10:
        min_val_loss=cumulative_val_loss
        torch.save(model.state_dict(),"./model_val.pt")
    if min_val_loss>cumulative_val_loss:
        min_val_loss=cumulative_val_loss
        torch.save(model.state_dict(),"./model_val.pt")

    scheduler.step(cumulative_val_loss)


np.save("./train_loss_epochs.npy", loss_train_epochs)
np.save("./val_loss_epochs.npy", loss_val_epochs)


print("Finished train\n")
print("printing images\n")

fig=plt.figure(figsize=(6,5))
plt.plot(loss_train_epochs, c='navy', label='Train Loss')
plt.plot(loss_val_epochs, c='darkorange', label='Validation Loss')
plt.grid()
plt.legend()
plt.xlabel('Epochs')
plt.ylabel('MSE Loss')
plt.title('Losses')
fig.savefig("./loss.png")
