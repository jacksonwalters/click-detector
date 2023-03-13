import csv
import re

with open('1669365039000_channelABCD_2022-11-25_10-30-39-ecolocation.flac.csv', 'r', encoding='utf-8-sig') as f:
    reader = csv.reader(f)
    rows = list(reader)

from datetime import datetime
from datetime import timedelta

prev_timestamp = None
for row in rows:
    timestamp_str = row[0] # assume the timestamp is in the first column
    try:
        parsed_time = float(timestamp_str)
    except:
        print('b')
    print(parsed_time)
    if prev_timestamp is not None and parsed_time - prev_timestamp > 10:
        row[0].append(['', '', '', '', ''])  # add a new empty row
    
    try:
        prev_timestamp = float(parsed_time)
    except:
        print("a")



with open('modified.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerows(rows) 
