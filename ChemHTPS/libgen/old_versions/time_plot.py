import matplotlib.pyplot as plt
infile=open('ttest_libgen.out')
time_list=[]
speedup_list=[]
efficiency_list=[]
procs_list=[96,80,72,60,48,36,24,12,10,8,6,4,2,1]
#line_no=0
for i,line in enumerate(infile):
    
    if 'Total' in line:
        line_parts=line.split('time_taken')
        time_taken=float(line_parts[1][:-1])
        
        time_list.append(time_taken)
seq_time=time_list[len(time_list)-1]

for i,time_taken in enumerate(time_list):
    speedup=seq_time/time_taken
    speedup_list.append(speedup)
    efficiency_list.append(speedup/procs_list[i])

            
plt.plot(procs_list,time_list,linewidth=2.0)
plt.scatter(procs_list,time_list,marker='o')
plt.title("Library generation of 93170 molecules",fontsize=22)
plt.xlabel("Number of processors used", fontsize=18)
plt.ylabel("Time taken (seconds)", fontsize=18)

plt.show()

plt.plot(procs_list,speedup_list,linewidth=2.0)
plt.scatter(procs_list,speedup_list,marker='o')
plt.title("Library generation of 93170 molecules",fontsize=22)
plt.xlabel("Number of processors used", fontsize=18)
plt.ylabel("Parallel speedup", fontsize=18)

plt.show()

plt.plot(procs_list,efficiency_list,linewidth=2.0)
plt.scatter(procs_list,efficiency_list,marker='o')
plt.title("Library generation of 93170 molecules",fontsize=22)
plt.xlabel("Number of processors used", fontsize=18)
plt.ylabel("Parallel efficiency", fontsize=18)

plt.show()
