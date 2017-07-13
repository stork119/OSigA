for pid in $(ps -ef | grep "parallel_computing_run" | awk '{print $2}'); do kill -9 $pid; done

