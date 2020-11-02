import matplotlib.pyplot as plt
import seaborn as sns
import subprocess
import pandas as pd 

# Excecute command
def execute_command(cmd):
    try:
        subprocess.check_call([cmd], shell=True)
    except subprocess.CalledProcessError:
        raise RuntimeError('Command "%s" does not work' % cmd)
    except OSError:
        raise Exception('Command "%s" does not exist' % cmd)

def plot_hist(data, variable_name, save_path):
    sns.set_context("paper")
    sns.set_style("ticks")
    plt.figure(dpi=300)
    df = pd.Series(data, name=variable_name)
    sns.distplot(df, kde=False)
    plt.savefig(save_path)

def plot_box(df, col_name_dict, save_path):
    sns.set_context("paper")
    sns.set_style("ticks")
    plt.figure(dpi=300)
    df = df[list(col_name_dict)].rename(columns=col_name_dict)
    print(pd.melt(df).sample())
    sns.boxplot(x="variable", y="value", data=pd.melt(df), palette='viridis')
    plt.savefig(save_path)
