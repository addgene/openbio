import pandas as pd
import numpy as np
from sklearn import linear_model
import sklearn as sk


# Import csv, pick out useful columns, rename columns
def virus_model(csv_path):

    virus_data = pd.read_csv(csv_path)
    virus_data_rough = virus_data[['Unnamed: 0', 'Plasmid ID', 'Serotype', 'Plasmid Size', 'Team', 'STATUS ',
                                   'Final Titer GC/mL', '# SU', 'Total Yield (GC)', 'Vessel']]
    virus_data_rough = virus_data_rough.rename(columns={"Unnamed: 0": "plasmid_name", "Plasmid ID": "plasmid_id",
                                                        "Serotype": "serotype", "Plasmid Size": "plasmid_size",
                                                        "Team": "team", "STATUS ": "status",
                                                        "Final Titer GC/mL": "final_titer", "# SU": "su",
                                                        "Total Yield (GC)": "total_yield", "Vessel": "vessel"})

    # Remove rejected or on-hold preps
    virus_data_rough = virus_data_rough[virus_data_rough.status == "Online"]

    # Identify every unique plasmid feature from the plasmid names using set()  lower case before split
    all_feats = set()
    for name in virus_data_rough['plasmid_name']:
        name_split = name.split('-')[1:]
        for x in name_split:
            all_feats.add(x)

    # Convert to list
    all_feats = list(all_feats)

    # For each prep, identify which plasmid features are present
    for feat in all_feats:
        plas_feat = lambda row: 1 if feat in list(row.plasmid_name.split("-")[1:]) else 0
        virus_data_rough[feat] = virus_data_rough.apply(plas_feat, axis=1)

    # All serotypes that are currently being used for AAV Prep
    all_sero = ['AAV5', 'RepCap8', 'RepCapRh10', 'RepCap2', 'RepCap9', 'AAVrg', 'RepCap5', 'Repcap1', 'AAV1']

    # For each prep, identify which serotype was used
    for sero in all_sero:
        sero_feat = lambda row: 1 if sero in row.serotype else 0
        virus_data_rough[sero] = virus_data_rough.apply(sero_feat, axis=1)


    # Identify all sci-tech teams
    all_teams = set()
    for team in virus_data_rough['team']:
        all_teams.add(team)

    # Convert to list
    all_teams = list(all_teams)

    # Label each prep with team tag
    for team in all_teams:
        team_apply = lambda row: 1 if team == row.team else 0
        virus_data_rough[team] = virus_data_rough.apply(team_apply, axis=1)

    # Identify all vessels
    all_vessels = set()
    for vess in virus_data_rough['vessel']:
        all_vessels.add(vess)

    # Convert to list
    all_vessels = list(all_vessels)

    # Label each prep with vessel
    for this_vess in all_vessels:
        vess_apply = lambda row: 1 if this_vess == row.vessel else 0
        virus_data_rough[this_vess] = virus_data_rough.apply(vess_apply, axis=1)

    # Add features for square and cube of plasmid size
    virus_data_rough['size_squared'] = virus_data_rough['plasmid_size']**2
    virus_data_rough['size_cubed'] = virus_data_rough['plasmid_size']**3

    # Remove rows in which Total Yield is NaN or zero
    virus_data_rough = virus_data_rough.dropna()
    virus_data_rough = virus_data_rough[(virus_data_rough.total_yield != 0)]

    # Select columns to be used for X_data, reindex, and convert to an array
    X_data = virus_data_rough.drop(['plasmid_name', 'plasmid_id', 'serotype', 'team', 'status', 'final_titer', 'su',
                                    'total_yield', 'vessel'], 1)
    X_data = X_data.reset_index(drop=True)
    X_data = sk.preprocessing.normalize(X_data.values)

    # Select labels, reindex, reformat to correct label format for sklearn
    Y_data = virus_data_rough.total_yield
    Y_data = Y_data.reset_index(drop=True)
    Y_data = pd.DataFrame(Y_data, columns=['total_yield'])
    Y_data = np.ravel(Y_data.values)

    # Print shape to verify
    print(X_data.shape)
    print(Y_data.shape)

    # Linear Regression model using X_data and Y_data
    model = linear_model.LinearRegression()
    model.fit(X_data, Y_data)
    print("Training accuracy is " + str(model.score(X_data, Y_data)))

    return model, all_feats, all_sero, all_teams, all_vessels

# Function to predict yield from a given prep


def predict_yield(model, plasmid_name_inst, serotype_name, plasmid_size, sci_team, vess_name='CS10'):

    # Initiate feature vector and add plasmid size
    temp_x = [plasmid_size]

    # Add plasmid features
    for feature in all_feats:
        feature_apply = (feature in plasmid_name_inst.split("-")[1:])
        temp_x.append(feature_apply)

    # Add serotypes
    for sero_inst in all_sero:
        sero_apply = sero_inst in serotype_name
        temp_x.append(sero_apply)

    # Add team
    for curr_team in all_teams:
        team_add = (curr_team == sci_team)
        temp_x.append(team_add)

    # Add vessel
    for curr_vess in all_vessels:
        vess_add = (curr_vess == vess_name)
        temp_x.append(vess_add)

    # Add square and cube plasmid_size
    temp_x.append(plasmid_size ** 2)
    temp_x.append(plasmid_size ** 3)

    # Convert from 1D array to 2D array
    temp_x = [temp_x]

    # Normalize inputs
    temp_x = sk.preprocessing.normalize(temp_x)

    return model.predict(temp_x)


# Define file path for training data
virus_csv = "~/Desktop/AAV_2CS5mod.csv"

model_virus, all_feats, all_sero, all_teams, all_vessels = virus_model(virus_csv)

print("Michelle's predicted yield is " + str(predict_yield(model_virus, 'pAAV-EF1a-Nuc-flox-mCherry-EGFP', 'AAVrg', 6927, 'Kate', 'CS10')))
print("Kate's predicted yield is " + str(predict_yield(model_virus, 'pAAV-EF1a-Flpo', 'AAVrg', 6208, 'Kate', 'CS10')))
print("Erin's predicted yield is " + str(predict_yield(model_virus, 'pAAV-EF1a-Flpo', 'AAVrg', 6208, 'Erin', 'CS10')))


# Used to confirm object type
# if isinstance(Y_data, pd.DataFrame):
# print(True)
# else:
# print(False)