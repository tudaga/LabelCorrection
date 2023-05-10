#!/bin/bash

python 0_generate_features.py
python 1_generate_predictions.py
python 2_evaluate_predictions.py
