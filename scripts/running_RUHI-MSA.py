
import subprocess

files_to_run = ['JOHNPP_all', 'JOHNTP_all', 'TP_all' ] #files to analyse
filename_1 = 'JOHNPP_l1'

cmd = "main_pedro_saved.cpp"

#proc = subprocess.Popen(["./a.out", filename])
#proc.wait()


for filename in files_to_run:
    #print (filename)
    proc = subprocess.Popen([".././a.out", filename]) # run file a.out with filename as argument
    proc.wait()
