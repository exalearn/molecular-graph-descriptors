import os
import numpy as np
import csv
import pandas as pd
from ase.db import connect

import torch
#torch.nn.Module.dump_patches = True

def evaluate(
    args,
    model,
    train_loader,
    val_loader,
    test_loader,
    device,
    metrics,
    custom_header=None,
    to_kcal=False,
):

    header = []
    results = []
    #preds, inits =[], []
    if "train" in args.split:
        header += ["Train MAE", "Train RMSE"]
        results += evaluate_dataset(metrics, model, train_loader, device)
        df = predict_dataset(args, model, train_loader, device)
    #
    if "val" in args.split:
        header += ["Val MAE", "Val RMSE"]
        results += evaluate_dataset(metrics, model, val_loader, device)
        df = predict_dataset(args, model, val_loader, device)
    #
    if "test" in args.split:
        header += ["Test MAE", "Test RMSE"]
        results += evaluate_dataset(metrics, model, test_loader, device)
        df = predict_dataset(args, model, test_loader, device)
    #
    if custom_header:
        header = custom_header

    # unit conversion
    if to_kcal:
        results = [r * 23.06054 for r in results]

    eval_file = os.path.join(args.modelpath, "evaluation.txt")
    with open(eval_file, "w") as file:
        wr = csv.writer(file)
        wr.writerow(header)
        wr.writerow(results)

    df_file =  os.path.join(args.modelpath, "results_dataframe.csv")
    df.to_csv(df_file)

    '''pred_file = os.path.join(args.modelpath, "predictions.csv")
    with open(pred_file, "w") as file:
        wr = csv.writer(file)
        wr.writerows(preds)#wr.writerow(pred_head)

    inits_file = os.path.join(args.modelpath, "actuals.csv")
    with open(inits_file, "w") as file:
        wr = csv.writer(file)
        wr.writerows(inits)'''

def evaluate_dataset(metrics, model, loader, device):
    for metric in metrics:
        metric.reset()

    for batch in loader:
        #print(len(batch))
        batch = {k: v.to(device) for k, v in batch.items()}
        #batch = {k: v.to('cuda:0') for k, v in batch.items()}
        #print(batch)
        model = torch.nn.DataParallel(model.module)
        result = model(batch)

        for metric in metrics:
            metric.add_batch(batch, result)

    results = [metric.aggregate() for metric in metrics]
    return results

def predict_dataset(args, model, loader, device):
    results=[]
    initials=[]
    pred_file = os.path.join(args.modelpath, "predicted_coordinates.xyz")
    with open(pred_file,'a') as file:
        wr = csv.writer(file)
        for batch in loader:
            batch = {k: v.to(device) for k, v in batch.items()}
            model = torch.nn.DataParallel(model.module)
            outputs = model(batch)

            wr.writerow([str(len((batch['_atomic_numbers'].cpu().numpy()[0])))])
            wr.writerow(['Actual: '+str(batch['energy'].cpu().numpy()[0][0])+'  Predicted: '+str(outputs['energy'].cpu().numpy()[0][0])])

            for i in range(len(batch['_atomic_numbers'].cpu().numpy()[0])):
                if batch['_atomic_numbers'].cpu().numpy()[0][i]==8:
                    atom_type='O   '
                else:
                    atom_type='H   '
                line = atom_type+str(batch['_positions'].cpu().numpy()[0][i][0])+'  '+str(batch['_positions'].cpu().numpy()[0][i][1])+'  '+str(batch['_positions'].cpu().numpy()[0][i][2])
                wr.writerow([line])

            initials.append(batch['energy'].cpu().numpy())
            results.append(outputs['energy'].cpu().numpy())
    d = {'Actual':initials, 'Prediction':results}#coords,'Atomic Numbers':atomic_numbers}
    df = pd.DataFrame(data=d)
    return df
