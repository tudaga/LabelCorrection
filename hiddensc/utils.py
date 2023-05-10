import numpy as np

RANDOM_SEED = 42


def set_random_seed(seed: int = RANDOM_SEED):
    np.random.seed(seed)
    print(f'Random seed set to {seed}')


def print_module_versions(module_list):
    """Print module versions"""
    for module in module_list:
        print(f'{module.__name__:<20s}: {module.__version__}')
