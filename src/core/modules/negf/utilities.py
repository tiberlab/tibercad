import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
    
    
def read_tibercad_file_as_df(name):
    
    file = open(name, 'r')
    commentedLine = file.readline()

    while True:
        if commentedLine[0]=='#':
            titleLine = commentedLine
            commentedLine = file.readline()
        else:
            break

    titleLine = titleLine.split()[1:]

    numberOfColumns = len(commentedLine.split())

    columns = [[] for i in range(numberOfColumns)]

    for i in range(numberOfColumns):
        columns[i].append(float(commentedLine.split()[i]))

    for numericLine in file:
        for i in range(numberOfColumns):
            columns[i].append(float(numericLine.split()[i]))
    file.close()
    
    df = pd.DataFrame(data = np.array(columns).T, columns = titleLine)
    return df

def read_tibercad_file_as_array(name):
    
    file = open(name, 'r')
    commentedLine = file.readline()

    while 0==0:
        if commentedLine[0]=='#':
            titleLine = commentedLine
            commentedLine = file.readline()
        else:
            break

    titleLine = titleLine.split()[1:]

    numberOfColumns = len(commentedLine.split())

    columns = [[] for i in range(numberOfColumns)]

    for i in range(numberOfColumns):
        columns[i].append(float(commentedLine.split()[i]))

    for numericLine in file:
        for i in range(numberOfColumns):
            columns[i].append(float(numericLine.split()[i]))
    file.close()
    
    array = np.array(columns).T
    return array


def read_dat_file(name):
    array = []
            
    with open(name, "r") as f:
        for line in f.readlines():
            array.append(float(line.split()[0]))
    
    return np.array(array)


def plot_tibercad_file(name, couplesOfVariables, returnDF = False, legend = False):
    
    file = open(name, 'r')
    commentedLine = file.readline()

    while 0==0:
        if commentedLine[0]=='#':
            titleLine = commentedLine
            commentedLine = file.readline()
        else:
            break

    titleLine = titleLine.split()[1:]

    numberOfColumns = len(commentedLine.split())

    columns = [[] for i in range(numberOfColumns)]

    for i in range(numberOfColumns):
        columns[i].append(float(commentedLine.split()[i]))

    for numericLine in file:
        for i in range(numberOfColumns):
            columns[i].append(float(numericLine.split()[i]))
    file.close()
    
    for i in range(len(couplesOfVariables)):
        s = '%s vs %s' %(titleLine[couplesOfVariables[i][0]],titleLine[couplesOfVariables[i][1]])
        plt.plot(columns[couplesOfVariables[i][0]],columns[couplesOfVariables[i][1]], label = s)
    if legend:
        plt.legend()
    
    if returnDF:
        df = pd.DataFrame(data = np.array(columns).T, columns = titleLine)
        return df
    else:
        return