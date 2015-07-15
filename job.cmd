# @ shell=/bin/bash
#
# Sample script for LoadLeveler
#
# @ error = job.$(jobid).err
# @ output = job.$(jobid).out
# @ job_type = parallel
# @ node_usage= not_shared
# @ node = 1
# @ tasks_per_node = 1
# @ resources = ConsumableCpus(4)
# @ node_resources = ConsumableMemory(1gb)
# @ network.MPI = sn_all,not_shared,us
# @ wall_clock_limit = 24:00:00
# @ notification = complete
# @ notify_user = $(user)@rzg.mpg.de
# @ queue

echo "#### begin env"
hostname -f
free
grep -m2 "^model" /proc/cpuinfo
lscpu
env | sort
echo "#### end env"

# here we go

#     n ConsumableMemory
# 11250              1gb
# 15909              2gb
# 22500              4gb
# 31819              8gb
# 45000             16gb
# 63639             32gb
# 90000             64gb

n=11250

cd /ptmp/${USER}/
outfile=dsyev_ex_${n}_.${LOADL_STEP_ID}
/usr/bin/time -v /u/${USER}/devel/testev/dsyev_ex "${n}" 2>${outfile}.err >${outfile}.out
