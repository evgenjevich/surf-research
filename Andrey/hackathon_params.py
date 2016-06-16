import json


def hackParams(N, stp, dx, sweeps):
    
    params = {'N' : N,  'steps' : stp, 'dx' : dx, 'sweeps' : sweeps}

    with open('params.json', 'w') as fp:
        json.dump(params, fp)
        
