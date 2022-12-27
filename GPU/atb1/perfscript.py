from os import listdir
from os.path import isfile, join
import os
import re
import time

def getfilenames(dir_path, filenames, mainfiles):

    filenames[dir_path] = []

    for f in listdir(dir_path):

        if isfile(join(dir_path, f)):

            if( f.endswith( "main.cu" ) ):
                mainfiles[dir_path] = join( dir_path, f )
            else:
                filenames[dir_path].append( join(dir_path, f) )
        
        else:

            getfilenames( join(dir_path, f), filenames, mainfiles )


def runscript():
    
    dir_path = os.path.dirname(os.path.realpath(__file__))

    filenames = dict()
    mainfiles = dict()

    Ni = "1024"
    Nj = "1024"
    Nk = "1024"

    getfilenames( dir_path, filenames, mainfiles )
    # getfilenames(dir_path + "/SM", filenames)
    # getfilenames(dir_path + "/Plain", filenames)

    for dirval, filevals in filenames.items():

        r = re.compile(".*\.cu")
        newlist = list(filter(r.match, filevals)) # Read Note below
        # print(newlist)

        compile_command_prefix = "nvcc -O3 -o "
        compile_command_mid = " " + mainfiles[dirval] + " "

        for filename in newlist:

            execfilename = filename[:-3]
            run_cmd = execfilename + " " + Ni + " " + Nj + " " + Nk
            print(run_cmd)

            compile_command = compile_command_prefix + 
            execfilename + compile_command_mid + filename
            # print(compile_command)
            os.system( compile_command )

            time.sleep(3)
            os.system( run_cmd )



def runscript1():
    
    dir_path = os.path.dirname(os.path.realpath(__file__))

    filenames = dict()
    mainfiles = dict()

    Ni = "128"
    Nj = ""
    Nk = "128"

    getfilenames(dir_path, filenames, mainfiles)
    # getfilenames( dir_path + "/Plain_ijkunroll", filenames, mainfiles )
    # getfilenames(dir_path + "/SM", filenames)
    # getfilenames(dir_path + "/Plain", filenames)

    vals = [128, 256, 512, 1024, 2048, 4096]
    # vals = [8192, 8192*4, 8192*16, 8192*32]

    for val in vals:

        Nj = str(val)

        for dirval, filevals in filenames.items():

            r = re.compile(".*\.cu")
            newlist = list(filter(r.match, filevals)) # Read Note below
            # print(newlist)

            compile_command_prefix = "nvcc -O3 -o "
            compile_command_mid = " " + mainfiles[dirval] + " "

            for filename in newlist:

                execfilename = filename[:-3]
                run_cmd = execfilename + " " + Ni + " " + Nj + " " + Nk
                print(run_cmd)

                compile_command = compile_command_prefix + \
                execfilename + compile_command_mid + filename
                # print(compile_command)
                os.system( compile_command )

                time.sleep(3)
                os.system( run_cmd )


    # print(filenames)



if __name__ == "__main__":
    runscript1()