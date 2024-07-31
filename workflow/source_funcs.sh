function started_job_message() {
    job_message=$1
    echo "Started: $job_message"

    # print end time message
    start_time=$(date "+%Y.%m.%d.%H.%M")
    echo "Start time: $start_time"
}

function ended_job_message() {
    job_message=$1
    echo "Ended: $job_message"

    # print end time message
    end_time=$(date "+%Y.%m.%d.%H.%M")
    echo "End time: $end_time"
}
