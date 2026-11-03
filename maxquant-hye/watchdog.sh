#!/bin/bash
# Stops MaxQuant before the WSL disk, the Windows C: drive (which hosts the WSL vhdx) or RAM runs out.
# Start it right after launching MaxQuant; it waits up to 5 min for the process to appear.
NAME=MaxQuant
PAT=MaxQuantCmd.dll
LOG=/home/sysoft/mq-hye/watchdog.log
for i in $(seq 30); do pgrep -f "$PAT" >/dev/null && break; sleep 10; done
pgrep -f "$PAT" >/dev/null || { echo "$(date '+%F %T') $NAME not found (pattern $PAT) — watchdog exit" >> "$LOG"; exit 1; }
while pgrep -f "$PAT" >/dev/null; do
  root=$(df -BG --output=avail / | tail -1 | tr -dc 0-9)
  cdrv=$(df -BG --output=avail /mnt/c | tail -1 | tr -dc 0-9)
  mem=$(awk '/MemAvailable/{print int($2/1048576)}' /proc/meminfo)
  swp=$(awk '/SwapFree/{print int($2/1048576)}' /proc/meminfo)
  echo "$(date '+%F %T') root=${root}G C=${cdrv}G memAvail=${mem}G swapFree=${swp}G" >> "$LOG"
  # an empty reading would make [ -lt ] fail silently and switch that guard off: log it and skip that check this round
  for v in root cdrv mem swp; do
    [ -n "${!v}" ] || { echo "$(date '+%F %T') WARNING could not read $v" >> "$LOG"; printf -v "$v" 999999; }
  done
  if [ "$root" -lt 12 ] || [ "$cdrv" -lt 20 ] || { [ "$mem" -lt 3 ] && [ "$swp" -lt 4 ]; }; then
    echo "$(date '+%F %T') ALERT stopping $NAME (root=${root} C=${cdrv} mem=${mem} swap=${swp})" >> "$LOG"
    pkill -TERM -f "$PAT"; sleep 30; pkill -KILL -f "$PAT"
    exit 1
  fi
  sleep 60
done
echo "$(date '+%F %T') $NAME no longer running — watchdog exit" >> "$LOG"
