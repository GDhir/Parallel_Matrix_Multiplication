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
    # getfilenames( dir_path, filenames )
    getfilenames(dir_path + "/ijtile", filenames)
    getfilenames(dir_path + "/ktile", filenames)
    # getfilenames(dir_path + "/jktile", filenames)
    # getfilenames(dir_path + "/jtile", filenames)
    getfilenames(dir_path + "/ikjtile", filenames)
    # getfilenames(dir_path + "/itile", filenames)

    getfilenames(dir_path + "/plainkunroll", filenames)

    r = re.compile(".*\.c")
    newlist = list(filter(r.match, filenames)) # Read Note below
    # newlist.remove( "mmtt_main.c")
    # newlist.remove( "mmtt_seq.c" )
    # print(newlist)

    compile_command_prefix = "gcc -O3 -fopenmp -o "
    compile_command_mid = " mmtt_main.c mmtt_seq.c "

    for filename in newlist:

        execfilename = filename[:-2]
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