#!/bin/bash

mkdir "$TMPDIR"/GMM_ord/

cd "$HOME"/GMM_ord

sbatch -a 1-100 submit_jobs.sh
