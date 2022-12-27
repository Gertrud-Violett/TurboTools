import sys
import papermill
import argparse

parser=argparse.ArgumentParser()
parser.add_argument('input data filename', nargs='*', help='input data filename (without .xlsx)')
parser.add_argument('new notebook filename', nargs='*', help='new ipynb notebook filename (without .ipynb')
args=parser.parse_args()

#Read syntax and files
if __name__ == '__main__':
    args = sys.argv
    targetfile = args[1]
    notebookfile = args[2]

papermill.execute_notebook('GeneralPlot_MD.ipynb',notebookfile+'.ipynb',parameters=dict(filename=targetfile))