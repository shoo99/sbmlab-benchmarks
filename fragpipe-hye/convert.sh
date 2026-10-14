#!/bin/bash
# raw -> indexed mzML (vendor peak picking), 3 at a time
cd /home/sysoft/fp-hye
ls /home/sysoft/mq-hye/raw/*.raw | xargs -P 3 -I{} sh -c '/home/sysoft/.dotnet/dotnet ext/trfp/ThermoRawFileParser.dll -i "{}" -o mzml -f 2 > "mzml/$(basename {} .raw).log" 2>&1; echo "done {} $?"'
echo CONVERT_ALL_DONE
