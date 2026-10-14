#!/bin/bash
# Stops MaxQuant before the WSL disk, the Windows C: drive (which hosts the WSL vhdx) or RAM runs out.
LOG=/home/sysoft/mq-hye/watchdog.log
while pgrep -f MaxQuantCmd.dll >/dev/null; do
  root=$(df -BG --output=avail / | tail -1 | tr -dc 0-9)
  cdrv=$(df -BG --output=avail /mnt/c | tail -1 | tr -dc 0-9)
  mem=$(awk '/MemAvailable/{print int($2/1048576)}' /proc/meminfo)
  swp=$(awk '/SwapFree/{print int($2/1048576)}' /proc/meminfo)
  echo "$(date '+%F %T') root=${root}G C=${cdrv}G memAvail=${mem}G swapFree=${swp}G" >> $LOG
  if [ "$root" -lt 12 ] || [ "$cdrv" -lt 20 ] || { [ "$mem" -lt 3 ] && [ "$swp" -lt 4 ]; }; then
    echo "$(date '+%F %T') ALERT stopping MaxQuant (root=${root} C=${cdrv} mem=${mem} swap=${swp})" >> $LOG
    pkill -TERM -f MaxQuantCmd.dll; sleep 30; pkill -KILL -f MaxQuantCmd.dll
    exit 1
  fi
  sleep 60
done
echo "$(date '+%F %T') MaxQuant no longer running — watchdog exit" >> $LOG
