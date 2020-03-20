#!/usr/bin/env python
import subprocess
import sys

jobid = sys.argv[1]

output = str(subprocess.check_output("sacct -j %s --format State --noheader | head -1 | awk '{print $1}'" % jobid, shell=True).strip())

# running_status=["PENDING", "CONFIGURING", "COMPLETING", "RUNNING", "SUSPENDED"]

# if "COMPLETED" in output:
#   print("success")
# elif any(r in output for r in running_status):
#   print("running")
# else:
#   print("failed")

fail_status = ["BOOT_FAIL", "OUT_OF_MEMORY", "CANCELLED", "DEADLINE", "FAILED", "NODE_FAIL", "PREEMPTED", "TIMEOUT", "SUSPENDED"]

if "COMPLETED" in output:
  print("success")
elif any(r in output for r in fail_status):
  print("failed")
else:
  print("running")
  