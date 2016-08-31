#!/bin/sh
echo "Executing gorillas.R script in the background"
Rscript --vanilla gorillas.R > gorillas.log 2>&1 &
echo "Check the progress with command 'tail -f gorillas.log'"
echo "Check the processor usage with command 'top'"
## End of script