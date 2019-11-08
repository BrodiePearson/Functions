import os

def createFolder(dirname):
    try:
        if not os.path.exists(dirname):
            os.makedirs(dirname)
    except OSError:
        print ('Error: Creating directory. ' +  dirname)
        

# Example
createFolder('./data/')
# Creates a folder in the current directory called data
