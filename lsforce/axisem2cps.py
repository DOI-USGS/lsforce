#!/usr/bin/env python

import os
import sys

import numpy as np
from pandas import read_table

# Limit model to this many layers; CPS supports up to 200
NUM_LAYERS = 195

# Model comment "can be no more than 80 characters long"
MAX_CHARACTERS = 80

# Four-space tab for help message and output model file
TAB = '    '

# From https://ds.iris.edu/ds/products/emc-syngine-retired/
MODELS = dict(
    ak135f='https://ds.iris.edu/media/product/emc-syngine/files/1dmodel_ak135f.txt',
    iasp91='https://ds.iris.edu/media/product/emc-syngine/files/1dmodel_iasp91.txt',
    prem_crust20_ocean='https://ds.iris.edu/media/product/emc-syngine/files/1dmodel_PREM_ocean.txt',
    prem_ani='https://ds.iris.edu/media/product/emc-syngine/files/1dmodel_PREMani.txt',
    prem_iso='https://ds.iris.edu/media/product/emc-syngine/files/1dmodel_PREMiso.txt',
)

help_message = (
    'Convert model files used for IRIS Syngine AxiSEM runs to\nComputer Programs in '
    f'Seismology (CPS) format.\n\nUsage:\n{TAB}axisem2cps <model>\n\n'
    f'Where <model> is one of:\n{TAB}' + f'\n{TAB}'.join(MODELS.keys())
)


def axisem2cps(model_nickname):
    """Convert model files used for IRIS Syngine AxiSEM runs to CPS format.

    CPS = Computer Programs in Seismology. The currently supported AxiSEM models are:

        * ak135f
        * iasp91
        * prem_crust20_ocean
        * prem_ani
        * prem_iso

    Args:
        model_nickname (str): Name of a supported AxiSEM model
    """

    # Error out if model_nickname is invalid
    if model_nickname not in MODELS.keys():
        model_key_strings = [f"'{m}'" for m in MODELS.keys()]
        raise ValueError(
            f'`model_nickname` must be one of: ' + ', '.join(model_key_strings)
        )

    # Read in model cleanly
    model_url = MODELS[model_nickname]
    tmp = read_table(model_url, sep=r'\s+', header=4, comment='#')
    model = tmp[tmp.columns[:-1]]
    model.columns = tmp.columns[1:]

    # Only retain the first 6 columns (this discards vph, vsh, eta for models with those)
    model = model[model.columns[:6]]

    # Unit conversions
    model.radius = model.radius / 1000.0  # To km
    model.rho = model.rho / 1000.0  # To g/cm^3
    model.vpv = model.vpv / 1000.0  # To km
    model.vsv = model.vsv / 1000.0  # To km

    # Convert to layers
    model['H(KM)'] = np.hstack([np.diff(np.abs(model.radius - model.radius[0])), 0])
    model = model[model['H(KM)'] != 0].reset_index()
    model = model.drop('radius', axis=1)

    # Add columns
    model['ETAP'] = np.zeros(model.shape[0])
    model['ETAS'] = np.zeros(model.shape[0])
    model['FREFP'] = np.ones(model.shape[0])
    model['FREFS'] = np.ones(model.shape[0])

    # Rename remaining columns
    model = model.rename(
        columns=dict(
            rho='RHO(GM/CC)', vpv='VP(KM/S)', vsv='VS(KM/S)', qka='QP', qmu='QS'
        )
    )

    # Put columns in correct order
    COLUMNS = [
        'H(KM)',
        'VP(KM/S)',
        'VS(KM/S)',
        'RHO(GM/CC)',
        'QP',
        'QS',
        'ETAP',
        'ETAS',
        'FREFP',
        'FREFS',
    ]
    model = model[COLUMNS]

    # Truncate to NUM_LAYERS and find maximum depth
    model = model[:NUM_LAYERS]
    max_depth = np.sum(model['H(KM)'])  # [km]

    # Form model file description and ensure it's not too long
    description = (
        f'AxiSEM input file "{os.path.basename(model_url)}" truncated to {max_depth:g} km '
        f'({NUM_LAYERS} layers)'
    )
    assert (
        len(description) <= MAX_CHARACTERS
    ), f'Description is more than {MAX_CHARACTERS} long!'

    # Write to the file
    output_filename = f'{model_nickname}.mod'
    with open(os.path.join(os.getcwd(), output_filename), 'w') as f:

        def line(string):
            return f.write(string + '\n')

        line('MODEL.01')
        line(description)
        line('ISOTROPIC')
        line('KGS')
        line('SPHERICAL EARTH')
        line('1-D')
        line('CONSTANT VELOCITY')
        line('LINE08')
        line('LINE09')
        line('LINE10')
        line('LINE11')
        line(TAB.join(COLUMNS))
        for _, row in model.iterrows():
            line(
                TAB.join(
                    [
                        '{: =4.1f}'.format(row['H(KM)']),
                        '{: =7.4f}'.format(row['VP(KM/S)']),
                        '{: =7.4f}'.format(row['VS(KM/S)']),
                        '{: =7.4f}'.format(row['RHO(GM/CC)']),
                        '{: =7.1f}'.format(row['QP']),
                        '{: =7.1f}'.format(row['QS']),
                        '{: =3.1f}'.format(row['ETAP']),
                        '{: =3.1f}'.format(row['ETAS']),
                        '{: =3.1f}'.format(row['FREFP']),
                        '{: =3.1f}'.format(row['FREFS']),
                    ]
                )
            )

    print(
        f'{output_filename} created ({NUM_LAYERS} layers, maximum depth = {max_depth:g} km)'
    )


def main():
    """Wrapper to faciliate ``axisem2cps`` command."""

    # Print help message if the user messed up
    if len(sys.argv) != 2 or sys.argv[1] not in MODELS.keys():
        print(help_message)
        sys.exit(1)

    # Call function
    model_nickname = sys.argv[1]
    axisem2cps(model_nickname)


if __name__ == '__main__':
    main()
