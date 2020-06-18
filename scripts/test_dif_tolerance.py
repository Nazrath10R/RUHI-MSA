import subprocess

tolerance_to_test = ['20', '15', '10', '5']

cmd = "main_pedro_saved.cpp"
filename = "JOHNPP_l1"

for tolerance in tolerance_to_test:
    #print (filename)
    proc = subprocess.Popen(["./a.out", filename, tolerance]) # run file a.out with filename as argument
    proc.wait()
