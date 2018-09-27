WinSize=10000
StepSize=5000
MinSites=10

Rscript slidewindowFST.R ORT $WinSize $StepSize $MinSites &
Rscript slidewindowFST.R RAC $WinSize $StepSize $MinSites &

wait
