from os import listdir
from os.path import isfile, join
import os
import re
import time

def getfilenames(dir_path, filenames):

    for f in listdir(dir_path):

        if isfile(join(dir_path, f)):

            filenames.append(join(dir_path, f))
        
        else:

            getfilenames( join(dir_path, f), filenames )


def runscript():
    
    dir_path = os.path.dirname(os.path.realpath(__file__))

    filenames = []
    getfilenames( dir_path + "/Slim", filenames )
    # getfilenames(dir_path + "/SM", filenames)
    # getfilenames(dir_path + "/Plain", filenames)

    r = re.compile(".*\.cu")
    newlist = list(filter(r.match, filenames)) # Read Note below
    # print(newlist)

    compile_command_prefix = "nvcc -O3 -o "
    compile_command_mid = " atb_main.cu "

    for filename in newlist:

        execfilename = filename[:-3]
        run_cmd = execfilename
        # print(filename)

        compile_command = compile_command_prefix + execfilename + compile_command_mid + filename
        print(compile_command)
        os.system( compile_command )

        time.sleep(3)
        os.system( run_cmd )


    # print(filenames)


if __name__ == "__main__":
    runscript()

    # filenames = []
    # getfilenames( os.path.dirname(os.path.realpath(__file__)), filenames )

    # print(filenames)